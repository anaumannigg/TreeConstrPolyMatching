#include "../include/cgal_includes.h"
#include "../include/shapefile_io_operations.h"
#include "../include/localization.h"
#include "../include/graph_computations.h"
#include "../include/linear_program.h"
#include "../include/command_line_parser.h"
#include "../include/tree_computations.h"
#include "../include/logger.h"
#include "../include/binary_io.h"

#include <exception>
#include <filesystem>
namespace fs = std::filesystem;



//a thread for printing out the status of computation
void statusTHREAD(std::atomic<int>& processed_counter, int num_sets) {
    while (true) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        int processed = processed_counter.load();
        // Move the cursor up one line and clear it
        std::cout << "\x1b[1A\x1b[2K";
        std::cout << "Processed: " << processed << " / " << num_sets << " sets." << std::endl;
        if (processed >= num_sets) break; // Exit when all sets are processed
    }
}

void setDecompositionTHREAD(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::vector<Polygon_wh>& merged_polys,
                            std::pair<int,int> interval, int thread_id, const Localization& rtree1, const Localization& rtree2,
                            std::vector<std::vector<int>>& poly1_id_lookup, std::vector<std::vector<int>>& poly2_id_lookup, std::vector<int>& set_sizes,
                            std::string folderMERGED, std::string data_name, std::atomic<int>& processed_counter) {

    //specifiy path for binary file
    //since the number of decomposition threads is equal to the number of worker threads later,
    //each decomposition thread creates exactly one binary file for a worker thread to read later on
    std::string binary_path = folderMERGED + "data_thread" + to_string(thread_id);

    //ensure file does not exist yet to not append to already existing file
    std::ifstream file(binary_path);
    if (file.good()) {
        file.close();
        std::remove(binary_path.c_str());
    }

    for (int p = interval.first; p <= interval.second; p++) {
        try {
            //initate the polygon for the collision detection as the outer boundary, as we will only use this to determine
            //the groups, this is sufficient for the speed-up and saves time on checking hole intersections
            const Polygon_wh m_poly = merged_polys[p];

            //collect all intersecting polygons in both maps
            std::vector<int> group_indices1, group_indices2;
            rtree1.get_neighbors(m_poly, &group_indices1);
            rtree2.get_neighbors(m_poly, &group_indices2);

            //remember indices for lookup
            poly1_id_lookup[p] = group_indices1;
            poly2_id_lookup[p]= group_indices2;

            std::vector<Polygon_wh> group1, group2;
            for (const auto& i : group_indices1) group1.push_back(polys1[i]);
            for (const auto& i : group_indices2) group2.push_back(polys2[i]);

            //store size of set
            set_sizes[p] = group1.size() + group2.size();

            //write to binary_file
            writePolysToBinaryFile(group1,group2,binary_path,p);
        } catch (const std::exception& e) {
            std::cerr << "Unexpected error in decomposition of merged poly " << p << endl;
            set_sizes[p] = 0;
            writePolysToBinaryFile(std::vector<Polygon_wh>(), std::vector<Polygon_wh>(), binary_path,p);
        }

        processed_counter++;
    }
}

//run n : m map matching algorithm on the maps provided by file1 and file2 for each of the provided epsilon values
//returns Solution to append to global Solution
void solveConnectedSetTHREAD(std::string data_name, double lambda, SOLUTION_MODE sol_mode, TREE_BUILD_MODE tree_mode, int thread_id,
                             bool exploit_opt_props, std::vector<std::pair<int,int>>& set_sizes,
                             std::vector<Solution>& sols, std::vector<std::vector<double>>& execution_times, std::atomic<int>& processed_counter) {

    //get path of binary file for this thread
    std::string binary_path = "../input/" + data_name + "/merged/data_thread" + std::to_string(thread_id);

    //create binary reader object, which will parse the bin file line by line
    BinaryPolygonFileReader bin_reader(binary_path);



    //init set id field
    size_t set_id = 0;
    std::vector<Polygon_wh> polys1,polys2;

    //this thread should handle the sets given in the range of indices provided
    while(bin_reader.readNextSet(set_id, polys1,polys2)) {
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
                    set_sizes[set_id].first = n_vertices_after_pregrouping;


                    //compute edges of each graph
                    computeEdges(g_vec, rtree1, rtree2, mo, lambda,  polys1, polys2);

                    int n_vertices_after_simple_matches = 0;

                    //start precomputation of simple 1:1 matches, which can be excluded from the graph
                    for (auto &g: g_vec) {
                        precomputeSimpleMatches(&g, mo, lambda, &sol_local);
                        n_vertices_after_simple_matches += g.num_vertices();
                    }


                    //remember total amount of vertices after simple match precomputation for analysis
                    set_sizes[set_id].second = n_vertices_after_simple_matches;

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
                        LinearProgram::solveILP_trees(g_tree, polys1.size(), polys2.size(), &sol_local);
                    }
                    else if(sol_mode == SOLUTION_MODE::CANZAR3APPROX) {
                        //solve approximation algorithm
                        LinearProgram::solveViaCanzar(g_tree,polys1.size(),polys2.size(),sol_local);
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


        sols[set_id] = sol_local;

        std::chrono::steady_clock::time_point set_end = std::chrono::steady_clock::now();
        execution_times[set_id] = {setupTime,exploitTime,treesTime,solvingTime};

        //reset polys
        polys1.clear();
        polys2.clear();

        processed_counter++;
    }

}

//takes two maps and a third map that contains the merged version of both maps with a union of all polygons with a common boundary
//decomposes the maps w.r.t. the merged map into separated subsets and saves them in a subfolder 'tmp' of the current folder
void run_NtoM_with_PreDecomposition(const CommandLineOptions& options) {
    //extract options
    std::string data_name = options.dataset_name;
    int num_threads = options.num_threads;
    SOLUTION_MODE sol_mode = options.mode;
    TREE_BUILD_MODE tree_mode = options.tree_mode;


    //retreive paths
    std::string file1 = "../input/" + data_name + "/" + data_name + "_1";
    std::string file2 = "../input/" + data_name + "/" + data_name + "_2";




    std::string folderMERGED = "../input/" + data_name + "/merged/";
    //create folder of it does not exist yet
    if (!fs::exists(folderMERGED)) fs::create_directory(folderMERGED);



    //overall timekeeping
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //read polygons from shapefile
    std::vector <Polygon_wh> polys1, polys2, merged_polys;
    int num_invalid1 = ReadShapeFile(file1, polys1, false);
    int num_invalid2 = ReadShapeFile(file2, polys2, false);

    cout << "read " << polys1.size() << " polygons from file 1. (" << num_invalid1 << " ignored due to invalid geometries)" << endl;
    cout << "read " << polys2.size() << " polygons from file 2. (" << num_invalid2 << " ignored due to invalid geometries)" << endl;


    //if a file with the merged information is already provided, load it, else compute the union and save it for
    //later use
    std::string fileMERGED = "../input/" + data_name + "/" + data_name + "_merged";
    if (fs::exists(fileMERGED + ".shp")){
        int num_invalid_merged = ReadShapeFile(fileMERGED, merged_polys, false);
        cout << "read " << merged_polys.size() << " polygons from merged file. (" << num_invalid_merged << " ignored due to invalid geometries)"<< endl;

    }
    else {
        merged_polys = merge(polys1,polys2);
        writeToShapeFile(merged_polys, fileMERGED);
        cout << "computed " << merged_polys.size() << " components via merging." << endl;
    }



    //check if any polygon within the sets is non-simple, as this will lead to issues in the code
    //warn the user if a non-simple polygon is found
    for(auto* polyset: {&polys1,&polys2,&merged_polys}) {
        for(auto& p : *polyset)  {
            if(!p.outer_boundary().is_simple()) {
                //try to fix with usua
                auto fixed_boundary = PolygonOperations::fixIfNotSimple(p.outer_boundary());

                if (!fixed_boundary.is_simple()) {
                    cout << "WARNING: Found non-simple polygon in the input, that will be handled as its convex hull: ";
                    cout << "POLYGON((";
                    for(const auto& v : p.outer_boundary().vertices()) cout << std::fixed << v << ", ";
                    cout<<"))\n";

                    //replace polygon by the convex hull of
                    std::vector<Point> ch;
                    CGAL::convex_hull_2(p.outer_boundary().vertices_begin(),p.outer_boundary().vertices_end(),std::back_inserter(ch));

                    p.clear();
                    Polygon p_outer_b;
                    for(const auto& ch_p : ch) p_outer_b.push_back(ch_p);

                    p = Polygon_wh(p_outer_b);
                }
                else {
                    p = Polygon_wh(fixed_boundary);
                }



                assert(p.outer_boundary().is_simple() && "Could not simplify Polygon by replacing it with its convex hull!");
            }
        }
    }

    //take time of entire decomposition
    std::chrono::steady_clock::time_point decomp_begin = std::chrono::steady_clock::now();

	//init r-trees of both datasets in order for quicker neighbor finding
	Localization rtree1(polys1);
	Localization rtree2(polys2);

	//init original vertex ID lookup, as IDs change on grouping
    int num_components = merged_polys.size();
	std::vector<std::vector<int>> poly_id_lookup1(num_components), poly_id_lookup2(num_components);

    //set up intervals for the computation on the connected sets using multiple threads
    //the intervals will be used for the decomposition as well as the solving
    num_threads = std::min(num_threads, num_components);
    int sets_per_thread = (int)(num_components / num_threads);

    std::vector<std::pair<int,int>> thread_set_intervals;
    for(int t=0; t<num_threads-1;t++) {
        thread_set_intervals.emplace_back(t * sets_per_thread, (t + 1) * sets_per_thread - 1);
    }
    thread_set_intervals.emplace_back((num_threads-1) * sets_per_thread, num_components - 1);


    //loop over all merged polygons and create sets (<-> connected components in the intersection graph)
    std::vector<int> set_sizes(num_components);

    cout << "Starting Decomposition into connected sets.\n\n";

    //start status thread
    std::atomic<int> processed_counter(0);
    std::thread status_thread(statusTHREAD, std::ref(processed_counter),num_components);

    //start worker threads
    std::vector<std::thread> decomp_threads(num_threads);
    for(int t=0; t < num_threads; t++) {
        decomp_threads[t] = std::thread(setDecompositionTHREAD,
                                     std::ref(polys1), std::ref(polys2),std::ref(merged_polys),
                                     thread_set_intervals[t], t, std::ref(rtree1), std::ref(rtree2),
                                     std::ref(poly_id_lookup1), std::ref(poly_id_lookup2), std::ref(set_sizes),
                                     folderMERGED, data_name, std::ref(processed_counter));
    }

    //join threads
    for(int t=0; t < num_threads; t++) {
        decomp_threads[t].join();
    }
    //ensure the status thread finished
    processed_counter = num_components;
    status_thread.join();

	cout << "Completed inital decomposition into connected subsets." << endl;
    std::chrono::steady_clock::time_point decomp_end = std::chrono::steady_clock::now();
    double timing_decomp = std::chrono::duration_cast<std::chrono::nanoseconds>(decomp_end - decomp_begin).count() / 1e9;


    //execute for all lambdas
    for (const auto& lambda : options.lambdas) {
        //initialize logger
        Logger logger("../export/" + data_name + "/" + options.log_name + ".csv");
        logger.setInstanceName(data_name);
        logger.setLambda(lambda);
        logger.setExploitMode(options.exploit_opt_props);
        if (options.tree_mode == TREE_BUILD_MODE::INFORMED) logger.setTreeMode("merge");
        else logger.setTreeMode("kruskal");
        if (options.mode == SOLUTION_MODE::OPT) logger.setSolMode("opt");
        else if (options.mode == SOLUTION_MODE::CANZAR3APPROX) logger.setSolMode("canzar3approx");
        logger.setInputPolygons(polys1.size(),polys2.size());
        logger.setConnectedComponents(merged_polys.size());



        //make sure the target folder exists
        std::string target_folder = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100));
        if(!fs::exists(target_folder)) fs::create_directories(target_folder);

        //measure time of global execution
        std::chrono::steady_clock::time_point exec_begin = std::chrono::steady_clock::now();

        //initialize global solution
        Solution sol(polys1.size(),polys2.size());

        //remember for each instance, if the exploration could be performed within the time limits
        std::vector<bool> completed_exploration; completed_exploration.resize(num_components, true);
        std::vector<std::vector<double>> execution_times; execution_times.resize(num_components, std::vector<double>());




        std::vector<Solution> thread_solutions(num_components);
        std::vector<bool> thread_explored(num_components);
        std::vector<std::pair<int,int>> set_sizes_after_precomputations(num_components);

        cout << "Starting to solve for lambda: " << lambda << ".\n\n";
        //measure time
        std::chrono::steady_clock::time_point solving_begin = std::chrono::steady_clock::now();

        //start status thread
        processed_counter = 0;
        status_thread = std::thread(statusTHREAD, std::ref(processed_counter), num_components);

        //start worker threads
        std::vector<std::thread> set_threads(num_threads);
        for(int t=0; t < num_threads; t++) {
            set_threads[t] = std::thread(solveConnectedSetTHREAD,
                                         data_name, lambda, sol_mode, tree_mode,  t,
                                         options.exploit_opt_props, std::ref(set_sizes_after_precomputations),
                                         std::ref(thread_solutions), std::ref(execution_times), std::ref(processed_counter));
        }

        //join threads
        for(int t=0; t < num_threads; t++) {
            set_threads[t].join();
        }
        // Ensure the status thread finishes
        processed_counter = num_components;
        status_thread.join();

        //insert all part solutions
        for(int s=0; s < num_components; s++) {
            sol.insert(thread_solutions[s], poly_id_lookup1[s], poly_id_lookup2[s]);
        }


        //complete global matching (setting negative match IDs for unmatched polygons)
        sol.completeMatching();

        //this is where the pure computation ends, take time and write to logger
        std::chrono::steady_clock::time_point exec_end = std::chrono::steady_clock::now();
        cout << "Logging..." << endl;
        std::vector<double> timings; timings.resize(6,0.0);
        timings[0] = timing_decomp;
        timings[5] = std::chrono::duration_cast<std::chrono::nanoseconds>(exec_end - exec_begin).count() / 1e9;
        for (int i=1; i<5; i++) {
            for (int j=0; j< execution_times.size(); j++) {
                timings[i] += execution_times[j][i-1];
            }
        }
        logger.setTimings(timings);
        logger.setObjective(sol.getTargetValue());
        logger.setTotalNumMatches(sol.getMatchCount());
        logger.setMatchDistribution(sol.getMatchDistribution());
        logger.setObjectiveDistribution(sol.getObjectiveDistribution());
        logger.log();
        std::cout << "\x1b[1A\x1b[2K";
        cout << "Logging completed." << endl;

        //export the global solution
        std::string output_path1 = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "1_matched_NtoM";
        std::string output_path2 = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "2_matched_NtoM";
        std::string csv_output_path = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "_data_NtoM";
        std::string analysis_output_path = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "_analysis";

        //write labeled polygons to shapefiles
        writeToShapeFile(polys1, sol.getMatchIndices(0), sol.getMatchWeights(0), output_path1);
        writeToShapeFile(polys2, sol.getMatchIndices(1), sol.getMatchWeights(1), output_path2);
        //write analysis data to csvs
        writeToCSV(sol, csv_output_path);


        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Computed n:m matching and exported results. (" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s])" << endl;
        cout << "Objective = " << sol.getTargetValue() << endl;
    }

}

int main(int argc, char* argv[]) {
    try {
        CommandLineOptions options = parse_command_line(argc,argv);

        if (options.dataset_name.empty() ||options.lambdas.empty()) {
            std::cerr << "Error: Both -d and -l (or -lr) options are required." << endl;
            print_help();
            return 1;
        }

        run_NtoM_with_PreDecomposition(options);

        cout << "Completed matching " << options.dataset_name << " for lambda(s) ";
        for (const auto& l : options.lambdas) {cout << l << ", ";} cout << endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << endl;
        print_help();
        return 1;
    }


	return 0;
}
