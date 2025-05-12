#include "../include/tree_computations.h"

//CLASS TREECONSTRAINEDCANDIDATEGRAPH

void TreeConstrainedCandidateGraph::add_edge(int source, int target, std::optional<double> weight){
    auto e = boost::add_edge(source, target, this->g);
    this->g[e.first].id = this->max_edge_id++;
    if(weight.has_value()) {
        this->g[e.first].weight = weight.value();
    }
    else {
        this->g[e.first].weight = 0.0;
    }
}

void TreeConstrainedCandidateGraph::add_tree_edge(bool map, int source, int target) {
    boost::add_edge(this->corresponding_tree_vertex[source],this->corresponding_tree_vertex[target],this->trees[map]);
}

void TreeConstrainedCandidateGraph::delete_edge(int source, int target){
    auto source_vertex = vertex(source, this->g);
    auto target_vertex = vertex(target, this->g);

    remove_edge(source_vertex, target_vertex, this->g);
    remove_edge(target_vertex, source_vertex, this->g);
}

void TreeConstrainedCandidateGraph::computeBipartiteWeightedEdges(const ObjectiveInfo& obj_info, const MapOverlay& mo, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2) {
    //start out with left root
    //the function will fix a vertex in the left tree and scan through the right tree to find all vertices representing
    //intersecting polys
    assert(this->roots[0] != -1 && this->roots[1] != -1 && "could not compute bipartite edges due to invalid roots!");

    //display warning if combined measure is set but polygons are not provided
    if (obj_info.objective == OBJECTIVE::JACCARD_HAUSDORFF && (polys1.empty() || polys2.empty())) {
        std::cerr << "Cannot compute hausdorff distance without polygon references in edge weight computation!" << endl;
    }

    //in the case of a combined measure, the hausdorff distance needs to be normalized in the end, keep track
    //of this
    std::vector<std::pair<detail::edge_desc_impl<undirected_tag,unsigned long>,double>> hausdorff_distances_per_edge;
    double max_hausdorff_dist = std::numeric_limits<double>::lowest();
    double min_hausdorff_dist = std::numeric_limits<double>::max();


    for(const auto& left_v : this->trees[0].vertex_set()) {
        std::vector<int> poly1_ids = this->g[this->corresponding_graph_vertex[0][left_v]].referenced_polys;
        std::queue<int> right_v_queue;
        right_v_queue.push(this->roots[1]);
        while(!right_v_queue.empty()) {

            int right_v = right_v_queue.front();
            right_v_queue.pop();

            //check if represented polys intersect
            std::vector<int> poly2_ids = this->g[this->corresponding_graph_vertex[1][right_v]].referenced_polys;

            if(mo.doOverlap(poly1_ids,poly2_ids,0.0)) {
                double IoU = mo.getIoU(poly1_ids,poly2_ids);
                //add edge to graph
                auto e = boost::add_edge(this->corresponding_graph_vertex[0][left_v],this->corresponding_graph_vertex[1][right_v],this->g);

                //next, the edge weight needs to be computed, this depends on the set option for the objective
                if (obj_info.objective == OBJECTIVE::JACCARD) {
                    this->g[e.first].weight = IoU-obj_info.lambda;
                } else if (obj_info.objective == OBJECTIVE::JACCARD_HAUSDORFF) {
                    //set the weight to the weighted Jaccard index first
                    this->g[e.first].weight = obj_info.weights.first * IoU;

                    //here, the Hausdorff-distance needs to be computed additionally
                    //collect relevant polygons
                    std::vector<const Polygon_wh*> polys1_of_set,polys2_of_set;
                    for (const auto& id : poly1_ids) polys1_of_set.push_back(&polys1[id]);
                    for (const auto& id : poly2_ids) polys2_of_set.push_back(&polys2[id]);

                    //we need the inverted hausdorff-distance to actually represent higher <-> better
                    double hausdorff =  1 / Polygon_wh::hausdorff_distance(polys1_of_set,polys2_of_set);
                    //remember max and min value for normalization later
                    if (hausdorff > max_hausdorff_dist) max_hausdorff_dist = hausdorff;
                    else if (hausdorff < min_hausdorff_dist) min_hausdorff_dist = hausdorff;
                    //remember the Hausdorff-distance to normalize it in the end
                    hausdorff_distances_per_edge.emplace_back(e.first,hausdorff);
                }



                this->g[e.first].id = this->max_edge_id++;

                //push incident vertices to queue
                auto adjacent_v = boost::adjacent_vertices(right_v,this->trees[1]);
                for (auto av = adjacent_v.first; av != adjacent_v.second; ++av) {
                    right_v_queue.push(*av);
                }
            }
        }

    }

    //if combined measure is chosen, renormalize and add hausdorff distance to weights
    if (obj_info.objective == OBJECTIVE::JACCARD_HAUSDORFF) {
        for (const auto& edge : hausdorff_distances_per_edge) {
            double normalized_hausdorff = (edge.second - min_hausdorff_dist) / max_hausdorff_dist;
            this->g[edge.first].weight += obj_info.weights.second * normalized_hausdorff;

            //subtract lambda
            this->g[edge.first].weight -= obj_info.lambda;
        }
    }


}

int TreeConstrainedCandidateGraph::add_vertex(bool ref_map, std::vector<int> ref_polys){
    //add vertex in every graph
    int new_vertex = boost::add_vertex(this->g);
    this->g[new_vertex].referenced_map = ref_map;
    this->g[new_vertex].referenced_polys = std::move(ref_polys);

    int new_tree_vertex = boost::add_vertex(this->trees[ref_map]);
    this->corresponding_graph_vertex[ref_map].push_back(new_vertex);
    this->corresponding_tree_vertex.push_back(new_tree_vertex);

    return new_vertex;
}

void TreeConstrainedCandidateGraph::delete_vertex(int v_id){
    //get the referred tree vertex
    auto t = this->getTree(v_id);
    int v_tree = this->corresponding_tree_vertex[v_id];

    //remove vertex from bipartite graph
    clear_vertex(v_id, this->g);
    remove_vertex(v_id, this->g);

    //remove vertex from the corresponding tree as well
    clear_vertex(v_tree, this->trees[t]);
    remove_vertex(v_tree, this->trees[t]);

    //adapt internal references
    this->corresponding_tree_vertex.erase(this->corresponding_tree_vertex.begin() + v_id);
    this->corresponding_graph_vertex[t].erase(this->corresponding_graph_vertex[t].begin() + v_tree);

    //every vertex, that had a higher ID than v_id, now has an id of one less
    //this must be considered in the 'corresponding...'-arrays
    for(auto& ctv : this->corresponding_tree_vertex) {
        if(ctv > v_tree) ctv--;
    }
    for(auto& cgv : this->corresponding_graph_vertex[t]) {
        if(cgv > v_id) cgv--;
    }

    if(this->roots[t] > v_tree) this->roots[t]--;


}

VertexProps TreeConstrainedCandidateGraph::get_vertex(int v_id){
    return this->g[v_id];
}

Vertex TreeConstrainedCandidateGraph::get_root(bool map){
    return this->roots[map];
}

void TreeConstrainedCandidateGraph::set_root(bool map, int root) {
    this->roots[map] = this->corresponding_tree_vertex[root];
}

void TreeConstrainedCandidateGraph::set_roots(std::vector<int> roots) {
    this->roots = roots;
}

bool TreeConstrainedCandidateGraph::is_root(int v_id) {
    return this->corresponding_tree_vertex[v_id] == this->roots[this->getTree(v_id)];
}

Graph&  TreeConstrainedCandidateGraph::get_graph() {
    return this->g;
}

void TreeConstrainedCandidateGraph::printVertexOverview() {
    cout << "\n---TREE CONSTRAINED CANDIDATE GRAPH VERTEX OVERVIEW---\n\n";
    //iterate over all vertices
    for (int v = 0; v < this->g.vertex_set().size(); v++) {
        std::string map_name = !this->g[v].referenced_map ? "1" : "2";
        cout << "vertex " << v << " represents "<<map_name << " polygons: {";
        for (const auto& ref_poly : this->g[v].referenced_polys) cout << ref_poly << ",";
        cout << "}" << endl;
    }

    cout << "\n---END TREE CONSTRAINED CANDIDATE GRAPH VERTEX OVERVIEW---\n\n";
}

std::string TreeConstrainedCandidateGraph::StringVertexOverview() {
    std::string overview_str;
    overview_str.append("\n---TREE CONSTRAINED CANDIDATE GRAPH VERTEX OVERVIEW---\n\n");
    //iterate over all vertices
    for (int v = 0; v < this->g.vertex_set().size(); v++) {
        std::string map_name = !this->g[v].referenced_map ? "1" : "2";
        overview_str.append("vertex " + std::to_string(v) + " represents " + map_name + " polygons: {");
        for (const auto& ref_poly : this->g[v].referenced_polys) overview_str.append(std::to_string(ref_poly) + ",");
        overview_str.append("}\n");
    }

    overview_str.append("\n---END TREE CONSTRAINED CANDIDATE GRAPH VERTEX OVERVIEW---\n\n");
    return overview_str;
}

void TreeConstrainedCandidateGraph::printTrees() {
    cout << "\n--- TREES ---\n\n";
    cout << "TREE 0: root " << this->roots[0] << endl;
    cout << "represented vertices in bigraph:" << endl;
    for(int v=0; v< this->corresponding_graph_vertex[0].size(); v++) cout << v << " - " << this->corresponding_graph_vertex[0][v] << endl;
    write_graphviz(cout,this->trees[0]);
    cout << "\nTREE 1: root " << this->roots[1] << endl;
    cout << "represented vertices in bigraph:" << endl;
    for(int v=0; v< this->corresponding_graph_vertex[1].size(); v++) cout << v << " - " << this->corresponding_graph_vertex[1][v] << endl;
    write_graphviz(cout, this->trees[1]);
    cout << "\n---END TREES --- \n\n";
}

std::string TreeConstrainedCandidateGraph::StringTrees() {
    std::string trees_str;
    trees_str.append("\n--- TREES ---\n\n");
    trees_str.append("TREE 0: root " + std::to_string(this->roots[0]) + "\n");
    write_graphviz(cout,this->trees[0]);
    trees_str.append("\nTREE 1: root " + std::to_string(this->roots[1]) + "\n");
    write_graphviz(cout, this->trees[1]);
    trees_str.append("\n---END TREES --- \n\n");

    return trees_str;
}

void TreeConstrainedCandidateGraph::setRootsforSingleVertexTrees(){
    for(int map=0; map <2; map++) {
        if(num_vertices(this->trees[map]) == 1) {
            this->roots[map] = 0;
        }
    }
}

void TreeConstrainedCandidateGraph::updateLeafs() {
    for(int tree = 0; tree < 2; tree++) {
        std::vector<Vertex> t_leafs;
        auto vertices_g = vertices(this->g);
        for (auto it = vertices_g.first; it != vertices_g.second; ++it) {
            if(this->g[this->corresponding_graph_vertex[tree][*it]].referenced_polys.size() == 1)
                t_leafs.push_back(*it);
        }

        this->leafs[tree] = t_leafs;
    }
}

int TreeConstrainedCandidateGraph::getTree(Vertex v) {
    return this->g[v].referenced_map;
}

Vertex TreeConstrainedCandidateGraph::getParent(Vertex v) {
    this->setTreeDirection(true);
    auto t = this->getTree(v);

    auto v_tree = this->corresponding_tree_vertex[v];

    assert(boost::out_degree(v_tree,this->trees[t]) == 1 && " vertex has not just one parent, tree condition hurt!");

    auto upwards_edge = boost::out_edges(v_tree,this->trees[t]);

    //since the graph is a tree, there will always be exactly one edge from a child to its parent node
    auto root  = upwards_edge.first->m_target;

    return this->corresponding_graph_vertex[t][root];

}

void TreeConstrainedCandidateGraph::setTreeDirection(bool dir) {
    //first make sure a tree with edges exists within the object
    bool tree = false;
    if(num_edges(this->trees[0])>0) tree = false;
    else if(num_edges(this->trees[1])>0) tree = true;
    else return; //no tree has edges, method can return


    bool invert = false;

    //check if invertion is necessary
    if(!dir) {
        //direction should be
        if(out_degree(this->roots[tree],this->trees[tree])==0) invert=true;
    }else {
        if(out_degree(this->roots[tree],this->trees[tree])>0) invert=true;
    }

    if(invert) {
        for (const auto &tree: {0, 1}) {
            // Create a new graph to store the inverted edges
            DiGraph t_inverted(num_vertices(this->trees[tree]));

            // Iterate over the edges of the original graph
            DiGraph::edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = edges(this->trees[tree]); ei != ei_end; ++ei) {
                auto u = source(*ei, g);
                auto v = target(*ei, g);
                boost::add_edge(v, u, t_inverted); // Add the inverted edge (v, u) to the new graph
            }

            this->trees[tree] = t_inverted;
        }
    }
}

std::vector<Vertex> TreeConstrainedCandidateGraph::getPathToRoot(bool tree, int source, PathtoRootMode mode) {
    this->setTreeDirection(true);

    //assume mode is start at vertex
    Vertex cur = 0;
    if(source < this->g.vertex_set().size()) cur = this->corresponding_tree_vertex[source];

    //if mode is start at poly, then change current vertex
    if(mode == PathtoRootMode::START_AT_POLYGON_ID) {
        //look for a leaf vertex representing the polygon
        bool v_found = false;

        for (const auto &l: this->leafs[tree]) {
            auto it = std::find(this->g[this->corresponding_graph_vertex[tree][l]].referenced_polys.begin(),
                                this->g[this->corresponding_graph_vertex[tree][l]].referenced_polys.end(),
                                source);
            if (it != this->g[this->corresponding_graph_vertex[tree][l]].referenced_polys.end()) {
                v_found = true;
                cur = l;
            }
        }

        //if poly is not represented by any vertex in the graph, return empty vector
        if(!v_found) return {};
    }
    else if(mode == PathtoRootMode::START_AT_VERTEX) assert(cur>=0 && cur<num_vertices(this->trees[tree]));


    std::vector<Vertex> path;
    path.push_back(this->corresponding_graph_vertex[tree][cur]);
    while(cur != this->roots[tree]) {
        auto av = adjacent_vertices(cur,this->trees[tree]);
        assert(av.first+1 == av.second && "Invalid number of adjacent vertices in path collection, check tree direction!");
        cur = *(av.first);
        path.push_back(this->corresponding_graph_vertex[tree][cur]);
    }
    return path;
}

int TreeConstrainedCandidateGraph::getMaxEdgeID() {
    return this->max_edge_id;
}

std::vector<std::pair<Vertex,Vertex>> TreeConstrainedCandidateGraph::getLeafPairs() {
    //make sure the tree direction is downwards
    this->setTreeDirection(false);

    std::vector<std::pair<Vertex,Vertex>> leaf_pairs;

    for(int t=0; t<2; t++) {
        for(auto v  : this->trees[t].vertex_set()){
            auto out_edges = boost::out_edges(v, this->trees[t]);
            std::vector<Vertex> adjacent_vertices;
            for(auto e_it = out_edges.first; e_it != out_edges.second; e_it++) {
                adjacent_vertices.push_back(e_it->m_target);
            }

            if(adjacent_vertices.size() == 2) {
                if(boost::out_degree(adjacent_vertices[0],this->trees[t]) == 0 && boost::out_degree(adjacent_vertices[1],this->trees[t]) == 0) {
                    //found pairs
                    leaf_pairs.push_back(std::make_pair(this->corresponding_graph_vertex[t][adjacent_vertices[0]],
                                                        this->corresponding_graph_vertex[t][adjacent_vertices[1]]));
                }
            }
        }

    }

    return leaf_pairs;
}

std::vector<Vertex> TreeConstrainedCandidateGraph::getLeafs(int tree) {
    //make sure the trees are in downwards direction
    this->setTreeDirection(false);

    std::vector<Vertex> leafs;

    for(auto v : this->trees[tree].vertex_set()) {
        if(boost::out_degree(v,this->trees[tree]) == 0) leafs.push_back(v);
    }

    return leafs;
}

std::vector<Vertex> TreeConstrainedCandidateGraph::getSubtree(bool tree, int source) {
    std::vector<Vertex> subtree;

    //make sure the trees are in downwards direction
    this->setTreeDirection(false);

    std::stack<Vertex> st;
    //push children of source
    auto adj_v = boost::adjacent_vertices(this->corresponding_tree_vertex[source], this->trees[tree]);
    for(auto v_it = adj_v.first; v_it != adj_v.second; v_it++) {
        st.push(*v_it);
    }
    while(!st.empty()) {
        Vertex v = st.top();
        st.pop();
        subtree.push_back(this->corresponding_graph_vertex[tree][v]);
        //collect all children
        auto adj_v = boost::adjacent_vertices(v, this->trees[tree]);
        for(auto v_it = adj_v.first; v_it != adj_v.second; v_it++) {
            st.push(*v_it);
        }
    }

    return subtree;
}

//END CLASS TREECONSTRAINEDCANDIDATEGRAPH


BiGraph buildBiGraphFromIntersectionGraph(CandidateGraph& cg, const MapOverlay& mo) {
    BiGraph g_bi;

    //first insert all vertices, as those are equal
    Graph g_inter = *cg.get_graph();
    //remember vertices of left and right side for edge computation
    std::vector<int> left_vs,right_vs;
    for(int v_id = 0; v_id < cg.num_vertices(); v_id++) {
        //add vertex and copy attributes
        BiVertex v_new = add_vertex(g_bi);
        g_bi[v_new].id = v_new;
        g_bi[v_new].represented_polys = std::vector<int>();
        for(const auto& rp : cg.get_vertex(v_id).referenced_polys) g_bi[v_new].represented_polys.push_back(rp);
        g_bi[v_new].represented_map = cg.get_vertex(v_id).referenced_map;

        if(!g_bi[v_new].represented_map) left_vs.push_back(v_new);
        else right_vs.push_back(v_new);
    }

    //check if the intersection graph has edges (if it is called with an incomplete intersection graph,
    // the edges should be computed from scratch using the MapOverlay)
    if(num_edges(g_inter) > 0) {
        //bidirected only need to be inserted where the intersection graph has undirected edges
        for (auto ei = edges(g_inter); ei.first != ei.second; ++ei.first) {
            auto e = *ei.first;

            int u = boost::source(e, g_inter);
            int v = boost::target(e, g_inter);

            //compute edge weights of both directions of the bidirectional edge and insert them into the graph
            for (int dir = 0; dir < 2; dir++) {
                int source = dir == 0 ? u : v;
                int target = dir == 0 ? v : u;

                double inter_area = !g_inter[source].referenced_map ?
                                    mo.getIntersectionArea(g_inter[source].referenced_polys,
                                                           g_inter[target].referenced_polys)
                                                                    : mo.getIntersectionArea(
                                g_inter[target].referenced_polys, g_inter[source].referenced_polys);
                double source_area = mo.getArea(g_inter[source].referenced_map, g_inter[source].referenced_polys);

                //insert edge and set its weight
                double edge_weight = inter_area / source_area;
                auto e_new = add_edge(source, target, g_bi);
                if (e_new.second) {
                    g_bi[e_new.first].weight = edge_weight;
                }
            }
        }
    }
    else {
        //no edges in intersection graph, compute the edges doing parwise checks
        for(const auto& u: left_vs) {
            for(const auto& v : right_vs) {
                if(mo.doOverlap(g_inter[u].referenced_polys,g_inter[v].referenced_polys,0.0)) {
                    //compute edge weights of both directions of the bidirectional edge and insert them into the graph
                    for (int dir = 0; dir < 2; dir++) {
                        int source = dir == 0 ? u : v;
                        int target = dir == 0 ? v : u;
                        double inter_area = !g_inter[source].referenced_map ?
                                            mo.getIntersectionArea(g_inter[source].referenced_polys,
                                                                   g_inter[target].referenced_polys)
                                                                            : mo.getIntersectionArea(
                                        g_inter[target].referenced_polys, g_inter[source].referenced_polys);

                        double source_area = mo.getArea(g_inter[source].referenced_map,
                                                        g_inter[source].referenced_polys);

                        //insert edge and set its weight
                        double edge_weight = inter_area / source_area;
                        auto e_new = add_edge(source, target, g_bi);
                        if (e_new.second) {
                            g_bi[e_new.first].weight = edge_weight;
                        }
                    }
                }
            }
        }

    }
    return g_bi;
}

// Custom vertex writer to print all vertex properties
class custom_vertex_writer {
public:
    custom_vertex_writer(const BiGraph& g) : g(g) {}

    template <class Vertex>
    void operator()(std::ostream& out, const Vertex& v) const {
        const auto& vp = g[v];
        out << "[label=\"inclusion_score: " << vp.inclusion_score
            << "\\nrepresented_map: " << (vp.represented_map ? "1" : "0")
            << "\\nrepresented_polys: ";

        for (size_t i = 0; i < vp.represented_polys.size(); ++i) {
            out << vp.represented_polys[i];
            if (i < vp.represented_polys.size() - 1) {
                out << ", ";
            }
        }
        out << "\"]";
    }

private:
    const BiGraph& g;
};
void write_graphviz_with_properties(std::ostream& out, const BiGraph& g) {
    // Create a property map for vertex labels
    custom_vertex_writer vw(g);

    // Create a property map for edge weights
    auto edge_writer = boost::make_label_writer(boost::get(&BiGraphEdgeProperties::weight, g));

    // Write the graph in Graphviz DOT format
    boost::write_graphviz(out, g,
                          vw,
                          edge_writer);
}

//TODO: improve in the case that the graph has been modified, then only a subset of scores has to be recomputed
//computes the inclusion scores of the vertices in the graph, defined as the sum of the max 2 incoming edge weights
//returns incoming edges per vertex for later reference
std::vector<std::vector<BiGraph::edge_descriptor>> computeVertexScores(BiGraph& g_bi) {
    Graph::vertex_iterator v, vend;

    //first collect all the incoming edges per vertex in O(|E|)
    std::vector<std::vector<BiGraph::edge_descriptor>> incoming_edges(num_vertices(g_bi));
    for(auto ei = edges(g_bi); ei.first != ei.second; ++ei.first) {
        auto e = *ei.first;
        incoming_edges[e.m_target].push_back(e);
    }

    for(int v_id = 0; v_id < num_vertices(g_bi); v_id++) {
        double score = 0.0;
        std::vector<double> weights;
        for(const auto& ie : incoming_edges[v_id]) weights.push_back(g_bi[ie].weight);

        //sort in reverse order
        std::sort(weights.begin(),weights.end(),std::greater<>());

        //set score to sum of two biggest weights, if no two edges exist set to 0
        if(weights.size() > 1) {
            score = weights[0] + weights[1];
        }
        g_bi[v_id].inclusion_score = score;
    }

    return incoming_edges;
}

TreeConstrainedCandidateGraph buildTreesFromBiGraph(BiGraph g_bi, CandidateGraph& cg, const MapOverlay& mo, double lambda) {
    //assert that the graph consists of one connected component
    std::vector<int> components(num_vertices(g_bi));
    int num_components = connected_components(g_bi,&components[0]);
    assert(num_components == 1 && "cannot build trees from a bipartite graph with > 1 components!");

    TreeConstrainedCandidateGraph cg_tree(cg);

    //perform vertex merges as long as possible
    while(num_vertices(g_bi) > 2) {
        //cout << "loop with graph size " << num_vertices(g_bi) << endl;
        //write_graphviz_with_properties(std::cout,g_bi);

        //compute the inclusion scores of the vertices
        auto incoming_edges = computeVertexScores(g_bi);

        //find maximum score (while skipping vertices with indegree < 2)
        double max_score = -1.0;
        Vertex max_vertex;
        for(const auto& v : g_bi.vertex_set()) {
            if(incoming_edges[v].size() > 1 && g_bi[v].inclusion_score > max_score) {
                max_score = g_bi[v].inclusion_score;
                max_vertex = v;
            }
        }


        //make sure max vertex has been found
        assert(max_score != -1.0 && "could not find max vertex!");

        // Use the adjacent_vertices function to get a pair of iterators
        auto adjVertices = boost::adjacent_vertices(max_vertex, g_bi);
        //get the two adjacent vertices with the highest edge weight
        std::vector<std::pair<Vertex,double>> adjVertices_top2;
        for (auto it = adjVertices.first; it != adjVertices.second; ++it) {

            double weight_to_max_vertex = g_bi[boost::edge(*it,max_vertex,g_bi).first].weight;

            adjVertices_top2.push_back(std::make_pair(*it, weight_to_max_vertex));
        }
        std::sort(adjVertices_top2.begin(),adjVertices_top2.end(),
                  [](const std::pair<Vertex,double>& a, std::pair<Vertex,double>& b) {
            return a.second > b.second;
        });
        if(adjVertices_top2.size() > 2) adjVertices_top2.resize(2);

        //check if there are adjacent vertices to be merged
        if(adjVertices_top2.size() > 1) {

            //set map
            bool map = g_bi[adjVertices_top2[0].first].represented_map;

            std::vector<int> merged_vertex_set = {(int) adjVertices_top2[0].first};
            std::vector<int> merged_poly_set = g_bi[adjVertices_top2[0].first].represented_polys;
            merged_poly_set.insert(merged_poly_set.end(),
                                   g_bi[adjVertices_top2[1].first].represented_polys.begin(),
                                   g_bi[adjVertices_top2[1].first].represented_polys.end());
            std::sort(merged_poly_set.begin(),merged_poly_set.end());

            //insert new vertex into the tree constrained candidate graph
            int v_new = cg_tree.add_vertex(map, merged_poly_set);

            //set new vertex as root (in the last iteration, the correct root will be set)
            cg_tree.set_root(map,v_new);

            //add tree edges to children
            for(const auto &av: adjVertices_top2) {
                cg_tree.add_tree_edge(map, v_new,g_bi[av.first].id);
            }


            //modify Bidirectional Graph accordingly, such that merging can be continued
            //this can be done in place, as one can add represented polys to the first vertex and remove the other one
            g_bi[adjVertices_top2[0].first].id = v_new;
            g_bi[adjVertices_top2[0].first].represented_polys.insert(
                    g_bi[adjVertices_top2[0].first].represented_polys.end(),
                    g_bi[adjVertices_top2[1].first].represented_polys.begin(),
                    g_bi[adjVertices_top2[1].first].represented_polys.end());
            std::sort(g_bi[adjVertices_top2[0].first].represented_polys.begin(),g_bi[adjVertices_top2[0].first].represented_polys.end());

            //add edges to newly created vertex (all vertices adjacent to the delete vertices should not also be adjacent to the merged vertex)
            for(const auto& ie : incoming_edges[adjVertices_top2[1].first]) {
                int u = ie.m_source;
                int v = adjVertices_top2[0].first;

                //if edges do not exist yet, compute them
                if(!boost::edge(u,v,g_bi).second) {

                    //compute edge weights of both directions of the bidirectional edge and insert them into the graph
                    for (int dir = 0; dir < 2; dir++) {
                        int source = dir == 0 ? u : v;
                        int target = dir == 0 ? v : u;

                        double inter_area = !g_bi[source].represented_map ?
                                            mo.getIntersectionArea(g_bi[source].represented_polys,
                                                                   g_bi[target].represented_polys)
                                                                            : mo.getIntersectionArea(
                                        g_bi[target].represented_polys, g_bi[source].represented_polys);
                        double source_area = mo.getArea(g_bi[source].represented_map,
                                                        g_bi[source].represented_polys);

                        //insert edge and set its weight
                        double edge_weight = inter_area / source_area;
                        auto e_new = add_edge(source, target, g_bi);
                        if (e_new.second) {
                            g_bi[e_new.first].weight = edge_weight;
                        }
                    }
                }
            }
            //delete all edges incident to vertex
            clear_vertex(adjVertices_top2[1].first, g_bi);
            //remove_vertex
            remove_vertex(adjVertices_top2[1].first, g_bi);
        }
    }
    //in some cases, one tree might only have one vertex, which should be set to be the root
    cg_tree.setRootsforSingleVertexTrees();

    return cg_tree;
}

//helper function to get the centroid of a polygon
Point getCentroid(const Polygon& p) {
    Point c(0,0);
    for (const auto& v : p.vertices()) {
        c.x() += v.x();
        c.y() += v.y();
    }
    c.x() /= (double)p.vertices().size();
    c.y() /= (double)p.vertices().size();

    return c;
}

//TODO: implement this function for benching against the merge method
TreeConstrainedCandidateGraph buildTreesViaKruskal(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, CandidateGraph& cg, MapOverlay& mo, double lambda) {
    //init return graph
    TreeConstrainedCandidateGraph cg_tree(cg);

    //init sets of polygons per map
    std::vector<std::vector<int>> poly_ids1,poly_ids2;
    std::vector<int> cg_id_lookup1,cg_id_lookup2;
    for (int v=0; v<cg.num_vertices(); v++) {
        if (!cg.get_vertex(v).referenced_map) {
            poly_ids1.push_back(cg_tree.get_vertex(v).referenced_polys);
            cg_id_lookup1.push_back(v);
        }
        else {
            poly_ids2.push_back(cg_tree.get_vertex(v).referenced_polys);
            cg_id_lookup2.push_back(v);
        }
    }

    //we will treat each dataset separately
    bool map = false;
    for (const auto& polys_and_ids : {make_tuple(polys1,poly_ids1,cg_id_lookup1),make_tuple(polys2,poly_ids2,cg_id_lookup2)}) {
        //cout << "MAP: " << map << endl;
        auto polys = get<0>(polys_and_ids);
        auto poly_sets = get<1>(polys_and_ids);
        auto id_lookup = get<2>(polys_and_ids);
        //build arrangement of polygons

        //collect segments of outer boundaries and holes
        std::vector<Segment> segments;
        for (const auto& poly : polys) {
            for (const auto& e : poly.outer_boundary().edges()) {
                segments.push_back(e);
            }
            for (const auto& h : poly.holes()) {
                for (const auto& e : h.edges()) {
                    segments.push_back(e);
                }
            }
        }

        //build arrangement
        Arrangement arr;
        CGAL::insert(arr, segments.begin(),segments.end());

        //init represented poly set per face
        for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
            f->set_data(-1);
        }

        //find the faces that correspond to sets of polygons (= vertices in the candidate graph)
        std::vector<Arrangement::Face_const_handle> face_of_poly(poly_sets.size());
        std::vector<bool> face_of_poly_set(poly_sets.size(), false);

        //set up vector of bounding boxes for quick first inclusion check
        std::vector<CGAL::Bbox_2> bboxes;
        for (const auto& p : polys) bboxes.push_back(p.bbox());

        for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
            if (f->has_outer_ccb()) {
                bool set_face = false;
                //essentially, we only need a sample point in the face
                auto edge = f->outer_ccb();
                Vector dir(edge->source()->point(),edge->target()->point());
                Point sample = edge->source()->point() + 0.5 * dir + 1e-8 * rotate(dir,PI/2);

                auto e_loop = edge;
                Polygon p_sanity;
                do {
                    p_sanity.push_back(edge->source()->point());
                }while (++edge != e_loop);
                edge++;
                while (!p_sanity.has_on_bounded_side(sample)) {
                    dir = Vector(edge->source()->point(),edge->target()->point());
                    sample = edge->source()->point() + 0.5 * dir + 1e-8 * rotate(dir,PI/2);
                    if (++edge == e_loop) {
                        std::cerr << "could not find sample in kruskal-like tree construction!" << endl;
                        std::cerr << "POLYGON((";
                        for (const auto& v : p_sanity) std::cerr << std::fixed << std::setprecision(5) << v << ", ";
                        std::cerr << "))" << endl;
                        break;
                    }
                }

                //cout << std::fixed << std::setprecision(5) << "DIR length: " << to_double(dir.squared_length()) << endl;
                //cout << "distance sample and source: "<< std::setprecision(8) << to_double(CGAL::squared_distance(sample,edge->source()->point())) << endl;

                //now we need to find the polygon that actually includes this sample
                for (int p_set_id = 0; p_set_id < poly_sets.size(); p_set_id++) {
                    for (const auto& p_id : poly_sets[p_set_id]) {
                        //rough check with bbox for speedup
                        if (sample.x() > bboxes[p_id].xmin() && sample.x() < bboxes[p_id].xmax()
                            && sample.y() > bboxes[p_id].ymin() && sample.y() < bboxes[p_id].ymax()) {
                            //lies in bounding box, now check for inclusion
                            if (polys[p_id].outer_boundary().has_on_bounded_side(sample)) {
                                //now we know the sample is within the outer boundary, exclude that it lies within a hole
                                bool lies_in_hole = false;
                                for (const auto& h : polys[p_id].holes()) {
                                    if (h.has_on_bounded_side(sample)) {
                                        lies_in_hole = true;
                                        break;
                                    }
                                }
                                if (!lies_in_hole) {
                                    //remember which set of polygons the face belongs to
                                    f->set_data(p_set_id);
                                    face_of_poly[p_set_id] = f;
                                    face_of_poly_set[p_set_id] = true;
                                    set_face = true;
                                    //other polygos of the set do not need to be considered anymore
                                    break;
                                }
                            }

                            }
                    }
                }
            }
        }

        //check if all pointers have been set
        int f_poly_id = 0;
        for (const auto& f : face_of_poly_set) {
            if (!f) {
                std::cerr << " face of poly set not set in kruskal-like tree construction!" << endl;
                for (const auto& p_ids : poly_sets[f_poly_id]) {
                    cout << "POLYGON((";
                    for (const auto& v : polys[p_ids].outer_boundary().vertices()) cout << std::fixed << std::setprecision(3) << v << ", " << endl;
                        cout << "))" << endl;
                }

            }
            f_poly_id++;
        }

        //cout << "collected faces, building graph of size " << poly_sets.size() << endl;
        //now we need to create the weighted neighborhood graph g for the tree building process
        Graph g(poly_sets.size());

        //init represented polygons
        auto [v_begin,v_end] = boost::vertices(g);
        int vid = 0;
        for (auto vit = v_begin; vit != v_end; vit++) {
            g[*vit].referenced_polys = poly_sets[*vit];
            g[*vit].id = id_lookup[vid++];
        }

        //precompute the centroids to avoid redundances
        std::vector<Point> centroids;
        for (const auto& p : polys) centroids.push_back(getCentroid(p.outer_boundary()));

        for (int p_set_id = 0; p_set_id < poly_sets.size()-1; p_set_id++) {
            //cout << "poly " << p_set_id << " of " << poly_sets.size() << endl;
            //do avoid redundant edges we will always only consider edges to polygons with higher ID
            //remember all already considered polys
            std::vector<bool> considered_polys(poly_sets.size()-p_set_id-1,false);

            //first collect polygons that share an edge with the current polygon
            auto face = face_of_poly[p_set_id];

            std::vector<int> incident_polys;
            std::vector<double> common_edge_lengths;

            if (face->has_outer_ccb()) {

                auto e = face->outer_ccb();
                auto e_loop = e;
                do {
                    int neighboring_face_poly_set_id = e->twin()->face()->data();
                    if (neighboring_face_poly_set_id != -1 && neighboring_face_poly_set_id != p_set_id && neighboring_face_poly_set_id > p_set_id) {
                        //neighboring face is assigned to polygon, remember common edge length for weighting in graph
                        //check if already found
                        int neighbor_index=-1;
                        for (int n=0; n<incident_polys.size(); n++) {
                            if (incident_polys[n] == neighboring_face_poly_set_id) {
                                neighbor_index = n;
                                break;
                            }
                        }
                        if (neighbor_index == -1) {
                            //neighboring poly not found yet, add
                            incident_polys.push_back(neighboring_face_poly_set_id);
                            common_edge_lengths.push_back(0.0);
                            neighbor_index = incident_polys.size()-1;
                        }
                        common_edge_lengths[neighbor_index] += sqrt(to_double(CGAL::squared_distance(e->source()->point(),e->target()->point())));
                    }
                }while (++e != e_loop);
            }

            // shift common edge lengths via f(x) = -1 / x , since we want the edges representing common edge lengths to be the most favorable (more than distances)
            // and the algorithm in each step takes the minimum weighted edge
            for (auto& cel : common_edge_lengths) {cel = - 1 / cel;}

            //add edge
            for (int ip = 0; ip < incident_polys.size(); ip++) {
                auto e = boost::add_edge(p_set_id,incident_polys[ip],g);
                g[e.first].weight = common_edge_lengths[ip];

                //mark as considered
                considered_polys[ip] = true;
            }

            //the remaining edges should represent the distances of the centroids of the polygons
            for (int ap = p_set_id+1; ap < poly_sets.size(); ap++) {
                if (!considered_polys[ap]) {
                    //polygon has not been considered yet, take centroid distance
                    //get minimum possible weight between two centroids of the two sets
                    double weight = std::numeric_limits<double>::max();
                    for (const auto& p1 : poly_sets[p_set_id]) {
                        for (const auto& p : poly_sets[ap]) {
                            double dist = sqrt(to_double(CGAL::squared_distance(centroids[p_set_id],centroids[p])));
                            if (dist < weight) weight = dist;
                        }
                    }
                    auto e = boost::add_edge(p_set_id, ap,g);
                    g[e.first].weight = weight;
                }
            }



        }

        //graph is complete now, build the trees via iteratively merging the shortest edges until the graph only contains one vertex
        while (boost::num_vertices(g) > 1) {
            //we are now looking for the minimum weighted edge
            boost::graph_traits<Graph>::edge_descriptor min_edge;
            double min_weight = std::numeric_limits<double>::max();
            auto es = boost::edges(g);
            for (auto e = es.first; e != es.second; ++e) {
                if (g[*e].weight < min_weight) {
                    min_weight = g[*e].weight;
                    min_edge = *e;
                }
            }
            //minimum weighted edge is found, corresponding vertices should be merged
            //we will keep the source vertex and update its represented polygons while deleting the target vertex
            // Iterate over adjacent vertices of the target
            auto target_adj = boost::adjacent_vertices(min_edge.m_target, g);
            for (auto v = target_adj.first; v != target_adj.second; ++v) {
                if (*v == min_edge.m_source) continue; // Skip the source itself
                auto[existing_edge, edge_exists] = boost::edge(min_edge.m_source, *v, g);

                if (edge_exists) {
                    // Update weight to the minimum of the two edges
                    g[existing_edge].weight = std::min(g[existing_edge].weight, g[boost::edge(min_edge.m_target, *v, g).first].weight);
                } else {
                    // Add a new edge with the weight from the target to this vertex
                    auto new_edge = add_edge(min_edge.m_source, *v, g);
                    g[new_edge.first].weight = g[boost::edge(min_edge.m_target, *v, g).first].weight;
                }
            }

            //update represented polygons of source
            g[min_edge.m_source].referenced_polys.insert(g[min_edge.m_source].referenced_polys.end(),
                g[min_edge.m_target].referenced_polys.begin(),g[min_edge.m_target].referenced_polys.end());


            //insert new vertex into the tree constrained candidate graph
            int v_new = cg_tree.add_vertex(map, g[min_edge.m_source].referenced_polys);
            //set new vertex as root (in the last iteration, the correct root will be set)
            cg_tree.set_root(map,v_new);

            //add tree edges to children
            for(const auto &ae: {min_edge.m_source,min_edge.m_target}) {
                cg_tree.add_tree_edge(map, v_new,g[ae].id);
            }
            //update source id
            g[min_edge.m_source].id = v_new;

            // Remove the target vertex
            clear_vertex(min_edge.m_target, g);
            remove_vertex(min_edge.m_target, g);

        }
        //set map fort second set of polgons
        map = true;
    }

    //we will use the centroids of the polygons to determine their closeness
     return cg_tree;
}
