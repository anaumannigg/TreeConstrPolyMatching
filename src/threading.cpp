#include "../include/threading.h"

// class THREADING

void Threading::solvingWorkerTHREAD(
                  std::queue<ConnectedSetTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  GRBEnv& env,
                  std::atomic<int>& processed_counter) {
    while (true) {
        ConnectedSetTask task;

        // Critical section for pulling a task
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            cv.wait(lock, [&] { return !tasks.empty() || done_flag; });

            if (tasks.empty()) break;

            task = tasks.front();
            tasks.pop();
        }

        // Run your solver function (adapted for one task)
        solveConnectedSet(task, env, processed_counter);

    }
}

void Threading::decompositionWorkerTHREAD(std::queue<DecompositionTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::vector<Polygon_wh>& merged_polys,
                  const Localization& rtree1, const Localization& rtree2,
                  std::atomic<int>& processed_counter) {
    while (true) {
        DecompositionTask task;

        // Critical section for pulling a task
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            cv.wait(lock, [&] { return !tasks.empty() || done_flag; });

            if (tasks.empty()) break;

            task = tasks.front();
            tasks.pop();
        }

        // Run your solver function (adapted for one task)
        decomposeIntoConnectedComponents(task, polys1,polys2,merged_polys,rtree1,rtree2,processed_counter);

    }
}

//a thread for printing out the status of computation
void Threading::statusTHREAD(std::atomic<int>& processed_counter, int num_sets) {
    while (true) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        int processed = processed_counter.load();
        // Move the cursor up one line and clear it
        std::cout << "\x1b[1A\x1b[2K";
        std::cout << "Processed: " << processed << " / " << num_sets << " sets." << std::endl;
        if (processed >= num_sets) break; // Exit when all sets are processed
    }
}

//run n : m map matching algorithm on the maps provided by file1 and file2 for each of the provided epsilon values
//returns Solution to append to global Solution
void Threading::solveConnectedSet(ConnectedSetTask task, GRBEnv& env, std::atomic<int>& processed_counter) {
    //retreive config info
    double lambda = task.lambda;
    std::string data_name = task.data_name;
    int batch_id = task.task_id;
    bool exploit_opt_props = task.options.exploit_opt_props;
    TREE_BUILD_MODE tree_mode = task.options.tree_mode;
    SOLUTION_MODE sol_mode = task.options.mode;

    //prepare local storage for thread solutions
    std::vector<int> set_ids;
    std::vector<Solution> sols;
    std::vector<std::vector<double>> execution_times;
    std::vector<std::pair<int,int>> set_sizes;

    //create object that defines the behavior for edge weight computations
    ObjectiveInfo obj_info_for_edge_weight_comp{task.options.objective,task.options.objective_weights,lambda};


    //get path of binary file for this thread
    std::string binary_path = "../input/" + data_name + "/tmp/data_batch" + std::to_string(batch_id) + ".bin";

    //create binary reader object, which will parse the bin file line by line
    BinaryPolygonFileReader bin_reader(binary_path);

    //init set id field
    size_t set_id = 0;
    std::vector<Polygon_wh> polys1,polys2;

    //this thread should handle the sets given in the range of indices provided
    while(bin_reader.readNextSet(set_id, polys1,polys2)) {
        //cout << "SOLVING SET " << set_id << endl;

        //remember set ID
        set_ids.push_back(set_id);

        //remember if set sizes matter (only for analysis in opt exploit)
        int set_sizes_before = set_sizes.size();

        //measure time for analysis output
        std::chrono::steady_clock::time_point set_begin = std::chrono::steady_clock::now();
        //initialize timings
        double setupTime = 0.0, exploitTime = 0.0, solvingTime = 0.0, treesTime = 0.0;

        //initialize sol_local storage
        Solution sol_local(polys1.size(), polys2.size());
        try {

            //collect and output all matched unmatched polygons from both sets
            std::vector <Polygon_wh> unconsidered1, considered1;
            std::vector <Polygon_wh> unconsidered2, considered2;


            //init r-trees of both datasets for spatial queries
            Localization rtree1(polys1);
            Localization rtree2(polys2);

            //precompute all areas within the arrangement of lines formed by the polygons in order to save computation time
            MapOverlay mo = MapOverlay(polys1, polys2);
            mo.assignFaces(polys1, polys2, rtree1, rtree2);



            //initialize vector of graphs for sol_local on threads
            std::vector <CandidateGraph> g_vec;

            //init lists to remember, which polygons were already considered
            std::vector<bool> poly1_visited;
            poly1_visited.resize(polys1.size(), false);
            std::vector<bool> poly2_visited;
            poly2_visited.resize(polys2.size(), false);

            int max_graph_size = 0;
            int number_of_cons_vertices = 0;

            //iterate over polygons in set 1 and find neighbors
            int index1 = 0;
            for (const auto &p1: polys1) {

                //if polygon is already visited_continue
                if (poly1_visited[index1]) {
                    index1++;
                    continue;
                } else poly1_visited[index1] = true;

                //build candidate graph of the connected component including the init polygon
                CandidateGraph cg;

                //add first node to graph
                std::vector<int> v_polys = {index1};
                int vertex_zero = cg.add_vertex(false, v_polys);

                //while new neighbors are found, iteratively collect and add neighbors in both maps
                //map switch variable switches between each loop (0, 1)
                bool map_switch = true;

                //remember neighbors from last loop
                std::vector<int> query_input;
                query_input.push_back(index1);

                //begin exploring the dataset and building the intersection graph
                do {
                    //vector to store all found neighbors
                    std::vector<int> neighbors;
                    //vector to store all neighbors, which are fully included in the same opposing polygon and thus should be cumulated already
                    std::vector <std::vector<int>> neighbors_grouped;

                    //decide, in which map to queue
                    Localization *rtree = !map_switch ? &rtree1 : &rtree2;
                    std::vector <Polygon_wh> polys = map_switch ? polys1 : polys2;

                    //queue each polygon in query input
                    for (const auto &q: query_input) {

                        //find all neighbors of query input q
                        std::vector<int> neighbors_all;
                        rtree->get_neighbors(polys[q], &neighbors_all);

                        //find all fully included neighbors, which should be introduced as a cumulative vertex right away
                        std::vector<int> neighbors_fully_included;
                        rtree->get_neighbors_fully_included(polys[q], &neighbors_fully_included);

                        //do not consider the fully included neighbors anymore
                        std::sort(neighbors_all.begin(), neighbors_all.end());
                        std::sort(neighbors_fully_included.begin(), neighbors_fully_included.end());

                        //if more than one neighbor is fully included, delete them from the all neighbors set as they should be a cumulated node right away
                        if (neighbors_fully_included.size() > 0) {
                            for (const auto &nfi: neighbors_fully_included)
                                neighbors_all.erase(std::remove(neighbors_all.begin(), neighbors_all.end(), nfi),
                                                    neighbors_all.end());
                        }

                        //add to neighbors and cumulative neighbor sets
                        neighbors.insert(neighbors.end(), neighbors_all.begin(), neighbors_all.end());
                        //EXPERIMENT:: also consider single fully included vertices as group
                        if (neighbors_fully_included.size() > 0) {
                            neighbors_grouped.push_back(neighbors_fully_included);
                        }


                    }

                    //make neighbor vector unique
                    std::sort(neighbors.begin(), neighbors.end());
                    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
                    std::sort(neighbors_grouped.begin(), neighbors_grouped.end());
                    neighbors_grouped.erase(std::unique(neighbors_grouped.begin(), neighbors_grouped.end()),
                                            neighbors_grouped.end());

                    //delete all grouped polygons from the single neighbors list, duplicates may happen due to exploration of multiple polygons beforehand
                    for (const auto &n_group: neighbors_grouped) {
                        for (const auto &n: n_group)
                            neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), n), neighbors.end());
                    }

                    //empty query input list for next loop
                    query_input.clear();

                    //remember, if any new nodes were added to the graph in this iteration
                    bool added_new_vertex = false;
                    //filter neighbors to only those, which were not found yet, also remember their indices within the corresponding list
                    for (const auto &n: neighbors) {
                        std::vector<bool> *search_list = !map_switch ? &poly1_visited : &poly2_visited;
                        if (!(*search_list)[n]) {
                            //mark as visited
                            (*search_list)[n] = true;
                            //save vertex info
                            std::vector<int> v_polys = {n};
                            //neighbor polygon has not been considered yet, create a node for it in the graph and set vertex info
                            int new_vertex = cg.add_vertex(map_switch, v_polys);

                            //remember query input for next loop
                            query_input.push_back(n);

                            added_new_vertex = true;
                        }
                    }

                    //filter grouped neighbors to only those, which were not found yet, also remember their indices within the corresponding list
                    for (const auto &n_group: neighbors_grouped) {
                        std::vector<bool> *search_list = !map_switch ? &poly1_visited : &poly2_visited;

                        //check if at least one of the polygons in the set is not visited yet
                        bool completely_visited = true;
                        for (const auto &n: n_group) {
                            if (!(*search_list)[n]) {
                                completely_visited = false;
                                break;
                            }
                        }

                        if (!completely_visited) {
                            //mark as visited
                            for (const auto &n: n_group) (*search_list)[n] = true;
                            //save vertex info
                            std::vector<int> v_polys = n_group;
                            //neighbor polygon has not been considered yet, create a node for it in the graph and set vertex info
                            int new_vertex = cg.add_vertex(map_switch, v_polys);

                            //remember query input for next loop
                            for (const auto &n: n_group) {
                                //EXPERIMENT: add another property: only further explore polygons that have another polygon that is at least 50% included in them
                                Localization *rtree_opp = map_switch ? &rtree1 : &rtree2;
                                std::vector <Polygon_wh> polys_opp = !map_switch ? polys1 : polys2;
                                std::vector<int> neighbors;
                                rtree_opp->get_neighbors_majorly_included(polys_opp[n], &neighbors, 0.5);
                                //only consider polygon for further exploration if there is a neighbor fulfilling the property
                                if (neighbors.size() > 1) {
                                    query_input.push_back(n);
                                }
                            }

                            added_new_vertex = true;
                        }

                    }

                    //break, once no more neighbors were found
                    if (!added_new_vertex) break;

                    map_switch = !map_switch;
                } while (true);


                Graph *g = cg.get_graph();
                if (g->vertex_set().size() > max_graph_size) max_graph_size = g->vertex_set().size();
                number_of_cons_vertices += g->vertex_set().size();

                //next step: build cumulated graph with nodes representing sets of polygons
                //if graph is small anyway, do not get into cumulated graph computation
                if (g->vertex_set().size() <= 2) {
                    //found isolated vertex, insert it into the sol_local as isolated
                    std::vector<int> o;
                    o.push_back(index1);
                    std::vector<int> a;

                    if (g->vertex_set().size() == 1) {
                        //do not add a match in this case, as singular polygon matches should be assigned negative indices via "completeMatching" in the end
                        //remember node1 as isolated node for output
                        unconsidered1.push_back(polys1[(*g)[0].referenced_polys[0]]);
                    } else {
                        //graph has size 2, poly2 should also be pushed back to sol_local
                        //poly index has to be vertex 1, as vertex 0 is always the initial poly1
                        for (const auto &index2: (*g)[1].referenced_polys) a.push_back(index2);

                        double q = mo.getIoU(o, a) - lambda;

                        //only form a match, if it leads to a better sol_local quality
                        if (q > 0) sol_local.addMatch(o, a, q);

                    }
                } else {
                    //graph has size >= 3, remember for thread traversal
                    g_vec.push_back(cg);

                }

                index1++;
            }

            //here ends the first big block of preparation for the final matching, take time
            std::chrono::steady_clock::time_point set_end = std::chrono::steady_clock::now();
            setupTime += std::chrono::duration_cast<std::chrono::nanoseconds>(set_end - set_begin).count() / 1e9;

            //only need to compute edge weights, cumulative vertices and solve ILP, if there are graphs with size > 3
            //those are now stored within g_vec
            if (g_vec.size() > 0) {
                //decide if preprocessing w.r.t. properties of optimal solutions should be done
                set_begin = std::chrono::steady_clock::now();
                if (exploit_opt_props) {
                    //first compute all inclusion properties between polygons
                    setInclusionProperties(g_vec, rtree1, rtree2, mo, lambda, polys1, polys2);

                    //next perform the pregrouping, where polygons of the same map, that are always in the same match,
                    //get grouped into a single cumulative vertex
                    int n_vertices_after_pregrouping = 0;
                    for (auto &g: g_vec) {
                        pregroupIncludedPolygons(&g, polys1, polys2, rtree1, rtree2, mo, lambda);
                        n_vertices_after_pregrouping += g.num_vertices();
                    }

                    //remember total vertex count after pregrouping for analysis
                    int set_size_before = n_vertices_after_pregrouping;


                    //compute edges of each graph
                    computeEdges(g_vec, rtree1, rtree2, mo, lambda,  polys1, polys2);

                    int n_vertices_after_simple_matches = 0;

                    //start precomputation of simple 1:1 matches, which can be excluded from the graph
                    for (auto &g: g_vec) {
                        precomputeSimpleMatches(&g, mo, lambda, &sol_local);
                        n_vertices_after_simple_matches += g.num_vertices();
                    }


                    //remember total amount of vertices after simple match precomputation for analysis
                    int set_size_after = n_vertices_after_simple_matches;
                    set_sizes.emplace_back(set_size_before, set_size_after);

                    //precomputation might have led to disconnected graphs, which can be deconstructed into connected components again
                    //create temporary memory for new g_threads vector
                    std::vector <CandidateGraph> g_vec_decomposed;
                    for (int g_id = 0; g_id < g_vec.size(); g_id++) {
                        Graph g = *g_vec[g_id].get_graph();

                        //skip for empty graphs
                        if (num_vertices(g) == 0) continue;

                        //store component per vertex
                        std::vector<int> comps(num_vertices(g));

                        //compute connected components
                        int num_components = connected_components(g, &comps[0]);

                        //if graph can be decomposed into more than one component, it should be split up
                        if (num_components > 1) {
                            //store copies of graph
                            std::vector <CandidateGraph> g_comps;
                            for (int c = 0; c < num_components; c++) g_comps.push_back(g_vec[g_id].copy());

                            //loop down to not change vertex numbers
                            for (int v = comps.size() - 1; v >= 0; v--) {
                                //the graphs in g_comps should represent the connected components, delete vertex from every graph except the one representing the component it is included in
                                for (int c = 0; c < num_components; c++) {
                                    if (comps[v] != c) {
                                        g_comps[c].delete_vertex(v);
                                    }

                                }
                            }

                            //add new graphs to vector
                            for (int c = 0; c < num_components; c++) {
                                if (g_comps[c].num_vertices() > 1) g_vec_decomposed.push_back(g_comps[c]);
                            }


                        } else {
                            //graph stays the same, insert to new list while filtering out empty graphs or those with only one vertex that might occur due to prematching simple matches
                            if (g_vec[g_id].num_vertices() > 1) g_vec_decomposed.push_back(g_vec[g_id]);
                        }
                    }

                    //reset g_threads pointer
                    g_vec.clear();
                    g_vec = g_vec_decomposed;
                }

                set_end = std::chrono::steady_clock::now();
                exploitTime += std::chrono::duration_cast<std::chrono::nanoseconds>(set_end - set_begin).count() / 1e9;

                //perform tree computation as well as solving the tree constrained matching problem per connected component
                for(auto& g : g_vec) {
                    //measure tree building time
                    set_begin = std::chrono::steady_clock::now();

                    TreeConstrainedCandidateGraph g_tree;

                    if (tree_mode == TREE_BUILD_MODE::INFORMED) {
                        //build Bidirectional Graph from the Candidate Graph
                        BiGraph g_bi = buildBiGraphFromIntersectionGraph(g,mo);
                        g_tree = buildTreesFromBiGraph(g_bi,g,mo,lambda);
                    }
                    else if (tree_mode == TREE_BUILD_MODE::KRUSKAL) {
                        g_tree = buildTreesViaKruskal(polys1,polys2,g, mo,lambda);
                    }

                    //compute weighted edges of desired objective function
                    g_tree.computeBipartiteWeightedEdges(obj_info_for_edge_weight_comp,mo,polys1,polys2);

                    set_end = std::chrono::steady_clock::now();
                    treesTime += std::chrono::duration_cast<std::chrono::nanoseconds>(set_end - set_begin).count() / 1e9;

                    //g_tree.printVertexOverview();
                    //g_tree.printTrees();
                    //cout << "graph has: " << boost::num_edges(g_tree.get_graph()) << " edges." << endl;

                    //measure actual solution time
                    set_begin = std::chrono::steady_clock::now();

                    //solve according to chosen option
                    if(sol_mode == SOLUTION_MODE::OPT) {
                        //solve ILP on tree constrained graph
                        LinearProgram::solveILP_trees(env, g_tree, polys1.size(), polys2.size(), &sol_local);
                    }
                    else if(sol_mode == SOLUTION_MODE::CANZAR3APPROX) {
                        //solve approximation algorithm
                        LinearProgram::solveViaCanzar(env, g_tree,polys1.size(),polys2.size(),sol_local);
                    }

                    set_end = std::chrono::steady_clock::now();
                    solvingTime += std::chrono::duration_cast<std::chrono::nanoseconds>(set_end - set_begin).count() / 1e9;


                }


            }
        }
        catch (std::exception &e) {
            std:cerr << "Unexpected error in solving set " << set_id << endl;
        }

        //complete the matching in the sol_local via assigning unmatched polygons negative indices
        sol_local.completeMatching();

        //set references in local solution for later retrieval in global solution
        sol_local.setReferences(polys1,polys2);

        sols.push_back(sol_local);

        std::chrono::steady_clock::time_point set_end = std::chrono::steady_clock::now();
        execution_times.push_back({setupTime,exploitTime,treesTime,solvingTime});

        //if set sizes did not differ, fill up with 0-info
        if (set_sizes_before == set_sizes.size()) set_sizes.emplace_back(0,0);

        //reset polys
        polys1.clear();
        polys2.clear();

        processed_counter++;
    }

    //write batch of solved sets to file
    std::string dump_path = "../input/" + data_name + "/tmp/sol_batch" + std::to_string(batch_id) + ".bin";
    write_thread_results_to_file(dump_path,set_ids,set_sizes,sols,execution_times);

}

void Threading::write_thread_results_to_file(const std::string& filename,
                                  const std::vector<int>& set_ids,
                                  const std::vector<std::pair<int,int>>& set_sizes,
                                  const std::vector<Solution>& sols,
                                  const std::vector<std::vector<double>>& execution_times) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }

    //cout << "writing: " << set_sizes.size() << " solutions into " << filename << " | sols: " << sols.size() << " | sizes: " << set_ids.size() << " | times: " << execution_times.size()  << endl;

    int size_count = sols.size();
    out.write(reinterpret_cast<const char*>(&size_count), sizeof(int));
    for (int i = 0; i < size_count; ++i) {
        out.write(reinterpret_cast<const char*>(&set_ids[i]), sizeof(int));  // index
        out.write(reinterpret_cast<const char*>(&set_sizes[i].first), sizeof(int));
        out.write(reinterpret_cast<const char*>(&set_sizes[i].second), sizeof(int));
        writeSolutionToBinaryFile(out, sols[i]);

        // write execution times
        const auto& times = execution_times[i];
        writeDoubleVecToBinaryFile(out, times);
    }

    out.close();
}

void Threading::read_thread_results_from_file(const std::string& filename,
                                   std::vector<int>& set_ids,
                                   std::vector<std::pair<int,int>>& set_sizes,
                                   std::vector<Solution>& sols,
                                   std::vector<std::vector<double>>& execution_times) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }

    int size_count;
    in.read(reinterpret_cast<char*>(&size_count), sizeof(int));

    for (int i = 0; i < size_count; ++i) {
        int set_id;
        int first, second;

        in.read(reinterpret_cast<char*>(&set_id), sizeof(int));
        in.read(reinterpret_cast<char*>(&first), sizeof(int));
        in.read(reinterpret_cast<char*>(&second), sizeof(int));

        set_ids.push_back(set_id);
        set_sizes.emplace_back(first, second);
        sols.push_back(readSolutionFromBinaryFile(in));
        execution_times.push_back(readDoubleVecFromBinaryFile(in));

    }

    in.close();
}

void Threading::decomposeIntoConnectedComponents(DecompositionTask task, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::vector<Polygon_wh>& merged_polys,
                            const Localization& rtree1, const Localization& rtree2, std::atomic<int>& processed_counter) {

    int batch_begin = task.interval.first;
    int batch_end = batch_begin + task.batch_size;


    while (true) {
        //break if all batches have been handled
        if (batch_begin > task.interval.second) break;

        //if batch end exceeds remaining polygon number, reduce
        if (batch_end > task.interval.second) batch_end = task.interval.second;

        //get batch ID
        int batch_id = batch_begin/task.batch_size;

        //set output path for batch
        std::string binary_path = "../input/" + task.data_name + "/tmp/data_batch" + to_string(batch_id) + ".bin";

        //ensure file does not exist yet to not append to already existing file
        std::ifstream file(binary_path);
        if (file.good()) {
            file.close();
            std::remove(binary_path.c_str());
        }

        for (int set = batch_begin; set <= batch_end; ++set) {
            try {
                //initiate the polygon for the collision detection as the outer boundary, as we will only use this to determine
                //the groups, this is sufficient for the speed-up and saves time on checking hole intersections
                const Polygon_wh m_poly = merged_polys[set];

                //collect all intersecting polygons in both maps
                std::vector<int> group_indices1, group_indices2;
                rtree1.get_neighbors(m_poly, &group_indices1);
                rtree2.get_neighbors(m_poly, &group_indices2);

                std::vector<Polygon_wh> group1, group2;
                for (const auto& i : group_indices1) group1.push_back(polys1[i]);
                for (const auto& i : group_indices2) group2.push_back(polys2[i]);

                //cout << "set " << p << ": ";
                //cout << "{"; for (const auto& i : group_indices1) cout << i <<", "; cout << "},";
                //cout << "{"; for (const auto& i : group_indices2) cout << i <<", "; cout << "}\n";

                //store size of set
                //set_sizes[p] = group1.size() + group2.size();

                //write to binary_file
                writePolysToBinaryFile(group1,group2,binary_path,set);

                //memory cleanup
                group1.clear();
                group2.clear();
                group_indices1.clear();
                group_indices2.clear();
            } catch (const std::exception& e) {
                std::cerr << "Unexpected error in decomposition of merged poly " << set << endl;
                //set_sizes[p] = 0;
                writePolysToBinaryFile(std::vector<Polygon_wh>(), std::vector<Polygon_wh>(), binary_path,set);
            }

            processed_counter++;
        }

        //increase batch ids
        batch_begin += task.batch_size;
        batch_end += task.batch_size;


    }


}

Solution Threading::combineThreadSolutions(const std::string& data_name, int num_batches, int num_polys1, int num_polys2, std::vector<int>& set_ids,
                                        std::vector<std::pair<int,int>>& set_sizes,
                                        std::vector<std::vector<double>>& execution_times) {
    Solution sol(num_polys1, num_polys2);

    //id offset
    int offset = set_ids.size();

    //read from all batches
    for (int b=0; b < num_batches; ++b) {
        //prepare memory
        std::vector<Solution> b_sols;

        //set path
        std::string filename = "../input/" + data_name + "/tmp/sol_batch" + to_string(b) + ".bin";

        //read sols, append to other fields
        read_thread_results_from_file(filename,set_ids,set_sizes,b_sols,execution_times);

        for (int c = 0; c < b_sols.size(); ++c) {
            sol.insert(b_sols[c]);
        }

        //update ID offset
        offset += b_sols.size();

    }

    return sol;
}


void Threading::cleanup_temp_folder(const std::string& folder_path) {
    namespace fs = std::filesystem;

    try {
        if (fs::exists(folder_path) && fs::is_directory(folder_path)) {
            fs::remove_all(folder_path);  // deletes contents + folder itself
        } else {
            std::cout << "Temp Folder doesn't exist: " << folder_path << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error during cleanup: " << e.what() << std::endl;
    }
}


// end class THREADING