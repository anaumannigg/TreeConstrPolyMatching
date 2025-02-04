#include "../include/linear_program.h"

//creates and solves an integer linear program w.r.t the nodes and edge weights within the given Graph g
void LinearProgram::solveILP(const Graph& g, int num_polys1, int num_polys2, Solution* solution) {



    try {
	//gurobi check
	GRBEnv env = GRBEnv(true);

    //deactivate console logging
    env.set("LogToConsole", "0");

    //env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables per edge in G (variable should be 1 when chosen and 0 when not chosen
    std::vector<GRBVar> var;
    std::vector<double> weight;

    //remember adjacent vertices per edge
    std::vector<std::vector<GRBVar>> var_adj1_v(g.vertex_set().size()), var_adj2_v(g.vertex_set().size());
    //remember order of edges
    std::vector<int> sources,targets;


    Graph::vertex_iterator v, vend;
    for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
        //to not count edges twice, only insert new variables once per outgoing edge of the left side of the graph
        if (!g[*v].referenced_map) {
            for (auto e = out_edges(*v, g); e.first != e.second; ++e.first) {
                GRBVar edge_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "");
                var.push_back(edge_var);
                double w = g[*e.first].weight;

                weight.push_back(w);

                //remember edge for atkis and polys1
                var_adj1_v[*v].push_back(edge_var);
                var_adj2_v[e.first->m_target].push_back(edge_var);

                //remember edge for solution retreival
                sources.push_back(*v);

                targets.push_back(e.first->m_target);
            }
        }
    }

    //build target function, which is the sum over all weighted edges
    GRBLinExpr target_func;
    target_func.addTerms(&weight[0], &var[0], var.size());

    //set the objective to maximizing the target function
    model.setObjective(target_func, GRB_MAXIMIZE);

    //add legality constraints
    for (bool map_switch : {false, true}) {
        int num_polys = !map_switch ? num_polys1 : num_polys2;
        for (int i = 0; i < num_polys; i++) {
            //collect all vertices, which refer to polygon i
            std::vector<int> referring_vertices;
            for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
                if (g[*v].referenced_map == map_switch && std::find(g[*v].referenced_polys.begin(), g[*v].referenced_polys.end(), i) != g[*v].referenced_polys.end()) {
                    //map fits and i is one of the referenced polygons, remember vertex
                    referring_vertices.push_back(*v);
                }
            }

            //if there are referring vertices in the current map, add expression
            if (referring_vertices.size() > 0) {
                //sum of variables referring to the respective nodes should not be greater than 1, s.t. every polygon gets picked at most once
                GRBLinExpr constr_expr;
                std::vector<std::vector<GRBVar>> adjacent_edges = !map_switch ? var_adj1_v : var_adj2_v;
                for (const auto& rv : referring_vertices) {
                    double coeff = 1;
                    for (const auto& adj_edge : adjacent_edges[rv]) {
                        constr_expr.addTerms(&coeff, &adj_edge, 1);
                    }
                }
                //add constraint to model, that the polygon may not be chosen twice
                model.addConstr(constr_expr, GRB_LESS_EQUAL, 1.0, "");
            }

        }
    }

    // Optimize model
    model.optimize();

    int matches_added = 0;
    //retreive solution: for each selected edge, put the referred polygons of source and target into the solution
    for (int i = 0; i < var.size(); i++) {
        if (var[i].get(GRB_DoubleAttr_X) == 1.0) {
            //edge is selected, put pair of source and target into solution including the weight of the connecting edge
            double weight = 0.0; bool weight_found = false;
            for (auto e = out_edges(sources[i], g); e.first != e.second; ++e.first) {
                if ((int)e.first->m_target == targets[i]) {
                    weight = g[*e.first].weight;
                    weight_found = true;
                    break;
                }
            }
            assert(weight_found && "did not find weight when iterating through edges!");

            solution->addMatch(g[sources[i]].referenced_polys, g[targets[i]].referenced_polys, weight);

            matches_added++;
        }
    }

    }
    catch (GRBException e) {
        cout << "Gurobi crashed with error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }

}

//creates and solves an integer linear program w.r.t the nodes and edge weights within the given Graph g
//takes tree constraints into account
void LinearProgram::solveILP_trees(TreeConstrainedCandidateGraph& cg_tree, int num_polys1, int num_polys2, Solution* solution) {
    Graph& g = cg_tree.get_graph();

    //make sure the graph is directed upwards from the leafs to ensure an efficient collection of constraints
    cg_tree.setTreeDirection(true);

    try {
        //gurobi check
        GRBEnv env = GRBEnv(true);

        // Restrict Gurobi to 1 thread per process
        env.set(GRB_IntParam_Threads, 1);

        //deactivate console logging
        env.set("LogToConsole", "0");

        //env.set("LogFile", "mip1.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create variables per edge in G (variable should be 1 when chosen and 0 when not chosen
        std::vector<GRBVar> var;
        std::vector<double> weight;

        //remember adjacent vertices per edge
        std::vector<std::vector<GRBVar>> var_adj1_v(g.vertex_set().size()), var_adj2_v(g.vertex_set().size());
        //remember order of edges
        std::vector<int> sources,targets;

        Graph::vertex_iterator v, vend;
        for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
            //to not count edges twice, only insert new variables once per outgoing edge of the left side of the graph
            if (!g[*v].referenced_map) {
                for  (auto e = out_edges(*v, g); e.first != e.second; ++e.first) {
                    GRBVar edge_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "");
                    var.push_back(edge_var);
                    double w = g[*e.first].weight;

                    weight.push_back(w);

                    //remember edge for atkis and polys1
                    var_adj1_v[*v].push_back(edge_var);
                    var_adj2_v[e.first->m_target].push_back(edge_var);

                    //remember edge for solution retreival
                    sources.push_back(*v);

                    targets.push_back(e.first->m_target);
                }
            }
        }

        //build target function, which is the sum over all weighted edges
        GRBLinExpr target_func;
        target_func.addTerms(&weight[0], &var[0], var.size());

        //set the objective to maximizing the target function
        model.setObjective(target_func, GRB_MAXIMIZE);

        //add legality constraints
        for (bool map_switch : {false, true}) {
            int num_polys = !map_switch ? num_polys1 : num_polys2;
            for (int i = 0; i < num_polys; i++) {

                //collect all vertices, which refer to polygon i
                std::vector<Vertex> referring_vertices = cg_tree.getPathToRoot(map_switch,i,PathtoRootMode::START_AT_POLYGON_ID);

                //if there are referring vertices in the current map, add expression
                if (referring_vertices.size() > 0) {
                    //sum of variables referring to the respective nodes should not be greater than 1, s.t. every polygon gets picked at most once
                    GRBLinExpr constr_expr;
                    std::vector<std::vector<GRBVar>> adjacent_edges = !map_switch ? var_adj1_v : var_adj2_v;
                    for (const auto& rv : referring_vertices) {
                        double coeff = 1;
                        for (const auto& adj_edge : adjacent_edges[rv]) {
                            constr_expr.addTerms(&coeff, &adj_edge, 1);
                        }
                    }
                    //add constraint to model, that the polygon may not be chosen twice
                    model.addConstr(constr_expr, GRB_LESS_EQUAL, 1.0, "");
                }

            }
        }

        // Optimize model
        model.optimize();

        int matches_added = 0;
        //retreive solution: for each selected edge, put the referred polygons of source and target into the solution
        for (int i = 0; i < var.size(); i++) {
            if (var[i].get(GRB_DoubleAttr_X) == 1.0) {
                //edge is selected, put pair of source and target into solution including the weight of the connecting edge
                double weight = 0.0; bool weight_found = false;
                for (auto e = out_edges(sources[i], g); e.first != e.second; ++e.first) {
                    if ((int)e.first->m_target == targets[i]) {
                        weight = g[*e.first].weight;
                        weight_found = true;
                        break;
                    }
                }
                assert(weight_found && "did not find weight when iterating through edges!");

                solution->addMatch(g[sources[i]].referenced_polys, g[targets[i]].referenced_polys, weight);

                matches_added++;
            }
        }


    }
    catch (GRBException e) {
        cout << "Gurobi crashed with error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }

}

//solves the linear program defined by the trees fractionally and returns the solution as a vector
//note that due to deletions in the Canzar Recursion, edges might be deleted. The returned vector thus has a constant
//length of the edge number of the initial graph. Edges that have been deleted have fractional value -1.0
std::vector<double> solveLP_trees_fractionally(TreeConstrainedCandidateGraph& cg_tree, int num_polys1, int num_polys2) {
    Graph& g = cg_tree.get_graph();

    //make sure the graph is directed upwards from the leafs to ensure an efficient collection of constraints
    cg_tree.setTreeDirection(true);

    std::vector<double> fractional_values(cg_tree.getMaxEdgeID(),-1);

    try {
        //gurobi check
        GRBEnv env = GRBEnv(true);

        // Restrict Gurobi to 1 thread per process
        env.set(GRB_IntParam_Threads, 1);

        //deactivate console logging
        env.set("LogToConsole", "0");

        //env.set("LogFile", "mip1.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create variables per edge in G (variable should be 1 when chosen and 0 when not chosen
        std::vector<GRBVar> var;
        std::vector<double> weight;

        //remember adjacent vertices per edge
        std::vector<std::vector<GRBVar>> var_adj1_v(g.vertex_set().size()), var_adj2_v(g.vertex_set().size());
        //remember order of edges
        std::vector<int> sources,targets;

        //remember edge ids for solution retrieval
        std::vector<int> edge_ids;

        for (auto ei = edges(g); ei.first != ei.second; ++ei.first) {
            auto e = *ei.first;

            edge_ids.push_back(g[e].id);

            GRBVar edge_var = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "");
            var.push_back(edge_var);
            double w = g[e].weight;

            weight.push_back(w);

            //make sure source is 1 and target is 2
            Vertex source = !g[e.m_source].referenced_map ? e.m_source : e.m_target;
            Vertex target = g[e.m_source].referenced_map ? e.m_source : e.m_target;

            //remember edge for atkis and polys1
            var_adj1_v[source].push_back(edge_var);
            var_adj2_v[target].push_back(edge_var);

            //remember edge for solution retrieval
            sources.push_back(source);

            targets.push_back(target);
        }

        //build target function, which is the sum over all weighted edges
        GRBLinExpr target_func;
        target_func.addTerms(&weight[0], &var[0], var.size());

        //set the objective to maximizing the target function
        model.setObjective(target_func, GRB_MAXIMIZE);

        //add legality constraints
        for (bool map_switch : {false, true}) {
            int num_polys = !map_switch ? num_polys1 : num_polys2;
            for (int i = 0; i < num_polys; i++) {
                //collect all vertices, which refer to polygon i
                std::vector<Vertex> referring_vertices = cg_tree.getPathToRoot(map_switch,i,PathtoRootMode::START_AT_POLYGON_ID);

                //if there are referring vertices in the current map, add expression
                if (referring_vertices.size() > 0) {
                    //sum of variables referring to the respective nodes should not be greater than 1, s.t. every polygon gets picked at most once
                    GRBLinExpr constr_expr;
                    std::vector<std::vector<GRBVar>> adjacent_edges = !map_switch ? var_adj1_v : var_adj2_v;
                    for (const auto& rv : referring_vertices) {
                        double coeff = 1;
                        for (const auto& adj_edge : adjacent_edges[rv]) {
                            constr_expr.addTerms(&coeff, &adj_edge, 1);
                        }
                    }
                    //add constraint to model, that the polygon may not be chosen twice
                    model.addConstr(constr_expr, GRB_LESS_EQUAL, 1.0, "");
                }

            }
        }

        // Optimize model
        model.optimize();

        //retreive solution
        for (int i = 0; i < var.size(); i++) {
            fractional_values[edge_ids[i]] = var[i].get(GRB_DoubleAttr_X);
        }

    }
    catch (GRBException e) {
        cout << "Gurobi crashed with error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }

    return fractional_values;

}

//the matching algorithm proposed in 'On Tree Constrained Matchings and Generalizations' (Algorithm 1)
std::vector<int> CanzarMatching(TreeConstrainedCandidateGraph cg_tree, int num_polys1, int num_polys2) {
    //get optimal fractional solution
    auto X = solveLP_trees_fractionally(cg_tree, num_polys1, num_polys2);

    //collect set of edges with fractional value 0
    std::vector<int> F_0;
    for(int x_id = 0; x_id < X.size(); x_id++) if(X[x_id] == 0.0) F_0.push_back(x_id);

    //check if set of 0-edges is not empty
    if(!F_0.empty()){
        //try to solve the matching problem again without the 0-weighted edges
        //create a deep copy of the tree
        TreeConstrainedCandidateGraph cg_tree_reduced(cg_tree);

        //collect sources and targets of all 0-edges
        std::vector<int> sources,targets;

        for (auto ei = edges(cg_tree_reduced.get_graph()); ei.first != ei.second; ++ei.first) {
            auto e = *ei.first;
            int e_id = cg_tree_reduced.get_graph()[e].id;

            if(std::binary_search(F_0.begin(),F_0.end(),e_id)) {
                //edge should be remembered for deletion
                sources.push_back(e.m_source);
                targets.push_back(e.m_target);
            }
        }

        //remove all 0-edges from the graph
        for(int e = 0; e<sources.size();e++) {
            cg_tree_reduced.delete_edge(sources[e],targets[e]);
        }

        //cout << "removed edges and recursed on graph with " << boost::num_edges(cg_tree_reduced.get_graph()) << " edges" << endl;
        //recurse
        auto M = CanzarMatching(cg_tree_reduced, num_polys1,num_polys2);
        return M;

    }

    //check if there exists an edge, such that all its conflicting edges are in sum < alpha = 3 in fractional value
    //conflicts are always the paths to the root
    int candidate_id = -1;
    double w_e = 0.0;
    std::vector<Graph::edge_descriptor> conflicting_edges;
    std::vector<int> conflicting_edge_ids;

    int dbg_edge_id = -1;
    for (auto ei = edges(cg_tree.get_graph()); ei.first != ei.second; ++ei.first) {
        dbg_edge_id++;
        auto e = *ei.first;

        //skip edges with decision variable value 1 as those have no conflict and thus no recursion is needed
        if(X[cg_tree.get_graph()[e].id] == 1.0) continue;

        //make sure source is 1 and target is 2
        Vertex source = !cg_tree.get_graph()[e.m_source].referenced_map ? e.m_source : e.m_target;
        Vertex target = cg_tree.get_graph()[e.m_source].referenced_map ? e.m_source : e.m_target;

        //collect paths to root
        std::vector<Vertex> root_path1 = cg_tree.getPathToRoot(0,source,PathtoRootMode::START_AT_VERTEX);
        std::vector<Vertex> subtree1 = cg_tree.getSubtree(0,source);
        std::vector<Vertex> root_path2 = cg_tree.getPathToRoot(1,target,PathtoRootMode::START_AT_VERTEX);
        std::vector<Vertex> subtree2 = cg_tree.getSubtree(1,target);

        //also collect all adjacent vertices of source and target as those pose conflicts as well
        auto adj_source = boost::adjacent_vertices(source,cg_tree.get_graph());
        std::vector<Vertex> adj2; for (auto ai = adj_source.first; ai != adj_source.second; ai++) {adj2.push_back(*ai);}
        auto adj_target = boost::adjacent_vertices(target,cg_tree.get_graph());
        std::vector<Vertex> adj1; for (auto ai = adj_target.first; ai != adj_target.second; ai++) {adj1.push_back(*ai);}



        std::vector<Vertex> vertices1 = root_path1;
        vertices1.insert(vertices1.end(),subtree1.begin(),subtree1.end());
        vertices1.insert(vertices1.end(),adj1.begin(),adj1.end());
        std::vector<Vertex> vertices2 = root_path2;
        vertices2.insert(vertices2.end(),subtree2.begin(),subtree2.end());
        vertices2.insert(vertices2.end(),adj2.begin(),adj2.end());

        //sort and make unique
        std::sort(vertices1.begin(),vertices1.end());
        vertices1.erase(std::unique(vertices1.begin(),vertices1.end()), vertices1.end());
        std::sort(vertices2.begin(),vertices2.end());
        vertices2.erase(std::unique(vertices2.begin(),vertices2.end()), vertices2.end());


        double conflict_sum = 0.0;

        //get all conflicts, meaning the edges adjacent to an 1 and/or 2 vertices out of the compiled set
        for (const auto& o : vertices1) {
            auto adj_edges = boost::out_edges(o,cg_tree.get_graph());
            for (auto ei = adj_edges.first; ei != adj_edges.second; ei++) {
                auto e_n = *ei;
                //skip the edge itself
                if (e_n == e) continue;
                conflict_sum += cg_tree.get_graph()[e_n].weight;
                conflicting_edges.push_back(e_n);
                conflicting_edge_ids.push_back(cg_tree.get_graph()[e_n].id);

            }
        }

        //note that we now only want to consider the edges that are not already considered (which may happen for edges incident to one of the vertices1)
        for (const auto& a : vertices2) {
            auto adj_edges = boost::out_edges(a,cg_tree.get_graph());
            for (auto ei = adj_edges.first; ei != adj_edges.second; ei++) {
                auto e_n = *ei;
                //skip the edge itself
                if (e_n == e) continue;
                //only add conflict if not considered yet
                if (std::find(conflicting_edges.begin(),conflicting_edges.end(),e_n) == conflicting_edges.end()) {
                    conflict_sum += cg_tree.get_graph()[e_n].weight;
                    conflicting_edges.push_back(e_n);
                    conflicting_edge_ids.push_back(cg_tree.get_graph()[e_n].id);
                }

            }
        }

        //if conflict sum is less than alpha = 3, a candidate is found
        if(conflict_sum < 3.0) {
            //cout << "setting candidate: " << e << endl;
            candidate_id = cg_tree.get_graph()[e].id;
            w_e = cg_tree.get_graph()[e].weight;
            break;
        }


    }

    //check if candidate has been found
    if(candidate_id > -1) {

        //modify weights of all edges in conflict with e
        for(const auto& n_e : conflicting_edges) {
            cg_tree.get_graph()[n_e].weight -= w_e;
        }

        //recurse
        auto M = CanzarMatching(cg_tree,num_polys1,num_polys2);

        //cout << "after recursion: " << boost::num_edges(cg_tree.get_graph()) << " edges in graph" << endl;

        //check if none of the edges in the conflicts are selected anymore
        std::sort(conflicting_edge_ids.begin(),conflicting_edge_ids.end());
        std::vector<int> edges_still_selected;
        std::set_intersection(conflicting_edge_ids.begin(),conflicting_edge_ids.end(),M.begin(),M.end(),std::back_inserter(edges_still_selected));
        if(edges_still_selected.empty()) {
            M.push_back(candidate_id);
            std::sort(M.begin(),M.end());
            return M;
        }

    }else {
        //arrived in base case
        //here, an empty matching should be returned
        //however, this point can also be reached, when the Linear Program returns with a fully integer solution
        //in this case, return the integer solution
        bool is_integer=true;
        for(const auto& x : X) {
            if(x!= 0.0 && x!=1.0 && x != -1.0) {
                is_integer = false;
                break;
            }
        }
        if(is_integer) {
            std::vector<int> M;
            for(int e=0; e < X.size(); e++) {
                if(X[e] == 1.0) M.push_back(e);
            }
            return M;
        }
        else return {};
    }
    return {};
}

void LinearProgram::solveViaCanzar(TreeConstrainedCandidateGraph& cg_tree, int num_polys1, int num_polys2, Solution& solution) {
    Graph& g = cg_tree.get_graph();
    //make sure the edges in the tree point upwards to the root
    cg_tree.setTreeDirection(true);

    //approximate Matching using Canzar's Algorithm
    auto M = CanzarMatching(cg_tree,num_polys1,num_polys2);

    //collect all found matches and add them to the solution
    // if no edges were selected, return
    if(M.size() == 0) return;
    int m_id = 0;
    //iterate over all edges in G
    for (auto ei = edges(cg_tree.get_graph()); ei.first != ei.second; ++ei.first) {
        auto e = *ei.first;
        int e_id = g[e].id;
        //M is sorted
        while(m_id < M.size()-1 && M[m_id] < e_id) m_id++;

        //check if edge is in matching
        if(e_id == M[m_id]) {
            //collect sets of polys to add to the solution
            bool source_map = g[e.m_source].referenced_map;
            std::vector<int> polys1 = !source_map? g[e.m_source].referenced_polys : g[e.m_target].referenced_polys;
            std::vector<int> polys2 = source_map? g[e.m_source].referenced_polys : g[e.m_target].referenced_polys;

            //add to solution
            solution.addMatch(polys1,polys2,g[e].weight);
        }
    }
}