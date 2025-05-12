#include "../include/threading.h"

#include <exception>
#include <filesystem>
namespace fs = std::filesystem;

//measuring the peak memory used by the process
size_t get_peak_memory_kb() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.substr(0, 6) == "VmHWM:") { // VmHWM = High Water Mark (peak resident set size)
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit;
            return std::stoul(value); // in kilobytes
        }
    }
    return 0;
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




    std::string folderTMP = "../input/" + data_name + "/tmp/";
    //create folder of it does not exist yet
    if (!fs::exists(folderTMP)) fs::create_directory(folderTMP);



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

                    int global_id = p.global_id;
                    p = Polygon_wh(p_outer_b);
                    p.global_id = global_id;
                }
                else {
                    int global_id = p.global_id;
                    p = Polygon_wh(fixed_boundary);
                    p.global_id = global_id;
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

	//get number of components
    int num_components = merged_polys.size();

    //set up intervals for the computation on the connected sets using multiple threads
    //the intervals will be used for the decomposition as well as the solving
    num_threads = std::min(num_threads, num_components);
    //set how many sets should be bundled into one task (1 batches)
    int sets_per_task = options.batch_size;
    int num_batches = num_components / sets_per_task;
    //there might be an additional smaller batch to consider all polygons
    bool last_batch_is_smaller = false;
    if (num_components % sets_per_task != 0) {
        num_batches++;
        last_batch_is_smaller = true;
    }

    std::vector<std::pair<int,int>> task_set_intervals;
    for(int t=0; t < num_batches;t++) {
        if (last_batch_is_smaller && t == num_batches - 1) {
            task_set_intervals.emplace_back((num_batches-1) * sets_per_task, num_components - 1);
            break;
        }
        task_set_intervals.emplace_back(t * sets_per_task, (t + 1) * sets_per_task - 1);
    }


    //create tasks to be handled by threads
    std::queue<DecompositionTask> decomp_tasks;
    for (const auto& t : task_set_intervals) {
        decomp_tasks.push((DecompositionTask){t, options.batch_size,data_name});
    }


    cout << "Starting Decomposition into connected sets.\n\n";

    //start status thread
    std::atomic<int> processed_counter(0);
    std::thread status_thread(Threading::statusTHREAD, std::ref(processed_counter),num_components);

    //init mutexes
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::atomic<bool> done_flag(false);

    //start worker threads
    std::vector<std::thread> decomp_threads;
    for(int t=0; t < num_threads; t++) {
        decomp_threads.emplace_back(Threading::decompositionWorkerTHREAD,
                                    std::ref(decomp_tasks),std::ref(queue_mutex),std::ref(cv), std::ref(done_flag),
                                     std::ref(polys1), std::ref(polys2),std::ref(merged_polys),
                                     std::ref(rtree1), std::ref(rtree2),std::ref(processed_counter));
    }

    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        done_flag = true;
    }
    cv.notify_all();

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

    //initialize one global Gurobi Environment to avoid collisions on license check in multiple threads
    GRBEnv env = GRBEnv(true);
    //deactivate console logging
    env.set("LogToConsole", "0");
    //start environment
    env.start();


    //execute for all lambdas
    for (const auto& lambda : options.lambdas) {
        //initialize logger
        Logger logger("../export/" + data_name + "/" + options.log_name + ".csv");
        logger.setInstanceName(data_name);
        logger.setLambda(lambda);
        logger.setExploitMode(options.exploit_opt_props);
        if (options.tree_mode == TREE_BUILD_MODE::INFORMED) logger.setTreeMode("informed");
        else logger.setTreeMode("kruskal");
        if (options.mode == SOLUTION_MODE::OPT) logger.setSolMode("opt");
        else if (options.mode == SOLUTION_MODE::CANZAR3APPROX) logger.setSolMode("canzar3approx");
        if (options.objective == OBJECTIVE::JACCARD) logger.setObjectiveMode("jaccard");
        else if (options.objective == OBJECTIVE::JACCARD_HAUSDORFF)
            logger.setObjectiveMode("jaccard_hausdorff_"+to_string(options.objective_weights.first) + "_" + to_string(options.objective_weights.second));
        logger.setInputPolygons(polys1.size(),polys2.size());
        logger.setConnectedComponents(merged_polys.size());



        //make sure the target folder exists
        //define subfolder export name depending on objective
        std::string target_subfolder_obj = "lambda" + to_string((int)(lambda * 100));
        if (options.objective == OBJECTIVE::JACCARD_HAUSDORFF) {
            target_subfolder_obj = "hd-jac-" + to_string((int)(options.objective_weights.first * 100))
                                    + "-" + to_string((int)(options.objective_weights.second * 100));
        }


        std::string target_folder = "../export/" + data_name + "/" + target_subfolder_obj;
        //if the jaccard hausdorff objective is chosen, change output folder

        if(!fs::exists(target_folder)) fs::create_directories(target_folder);

        //measure time of global execution
        std::chrono::steady_clock::time_point exec_begin = std::chrono::steady_clock::now();

        //prepare thread task pool
        std::queue<ConnectedSetTask> solving_tasks;
        for (int i=0; i< num_batches; i++) {
            solving_tasks.push((ConnectedSetTask){data_name, lambda, options,i});
        }


        cout << "Starting to solve for lambda: " << lambda << ".\n\n";
        //measure time
        std::chrono::steady_clock::time_point solving_begin = std::chrono::steady_clock::now();

        //start status thread
        processed_counter = 0;
        status_thread = std::thread(Threading::statusTHREAD, std::ref(processed_counter), num_components);


        //start worker threads
        done_flag = false;
        std::vector<std::thread> set_threads(num_threads);
        for(int t=0; t < num_threads; t++) {
            set_threads[t] = std::thread(Threading::solvingWorkerTHREAD,
                                         std::ref(solving_tasks),std::ref(queue_mutex),std::ref(cv),std::ref(done_flag),
                                         std::ref(env), std::ref(processed_counter));
        }

        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            done_flag = true;
        }
        cv.notify_all();

        //join threads
        for(int t=0; t < num_threads; t++) {
            set_threads[t].join();
        }

        // Ensure the status thread finishes
        processed_counter = num_components;
        status_thread.join();

        cout << "done, now combining solutions" << endl;

        //combine all part-solutions produced by the threads
        std::vector<int> set_ids;
        std::vector<std::pair<int,int>> set_sizes;
        std::vector<std::vector<double>> execution_times;
        Solution sol = Threading::combineThreadSolutions(data_name,num_batches,polys1.size(),polys2.size(),set_ids,
                                set_sizes,execution_times);

        cout << "completing matching..." << endl;

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
        logger.setMemoryUsage(get_peak_memory_kb());
        logger.setObjective(sol.getTargetValue());
        logger.setTotalNumMatches(sol.getMatchCount());
        logger.setMatchDistribution(sol.getMatchDistribution());
        logger.setObjectiveDistribution(sol.getObjectiveDistribution());
        logger.log();
        std::cout << "\x1b[1A\x1b[2K";
        cout << "Logging completed." << endl;

        //export the global solution
        std::string output_path1 = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "1_matched_NtoM";
        std::string output_path2 = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "2_matched_NtoM";
        std::string csv_output_path = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "_data_NtoM";
        std::string analysis_output_path = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "_analysis";

        //write labeled polygons to shapefiles
        writeToShapeFile(polys1, sol.getMatchIndices(0), sol.getMatchWeights(0), output_path1);
        writeToShapeFile(polys2, sol.getMatchIndices(1), sol.getMatchWeights(1), output_path2);
        //write analysis data to csvs
        writeToCSV(sol, csv_output_path);


        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Computed n:m matching and exported results. (" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s])" << endl;
        cout << "Objective = " << sol.getTargetValue() << endl;
    }

    //cleanup tmp storage
    Threading::cleanup_temp_folder("../input/" + data_name + "/tmp/");

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
