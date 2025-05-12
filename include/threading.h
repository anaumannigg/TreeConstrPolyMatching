#ifndef THREADING_H
#define THREADING_H

#include "command_line_parser.h"
#include "binary_io.h"
#include "cgal_includes.h"
#include "polygon_wh.h"
#include "shapefile_io_operations.h"
#include "localization.h"
#include "graph_computations.h"
#include "linear_program.h"

#include <string>
#include <condition_variable>


//define a struct for a decomposition task
struct DecompositionTask {
    std::pair<int, int> interval;
    int batch_size;
    std::string data_name;
};

//define a struct for a solving task
struct ConnectedSetTask {
    std::string data_name;
    double lambda;
    CommandLineOptions options;
    int task_id;
};

class Threading {
public:
    static void solvingWorkerTHREAD(
                  std::queue<ConnectedSetTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  GRBEnv& env,
                  std::atomic<int>& processed_counter);

    static void decompositionWorkerTHREAD(std::queue<DecompositionTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::vector<Polygon_wh>& merged_polys,
                  const Localization& rtree1, const Localization& rtree2, std::atomic<int>& processed_counter);

    static void statusTHREAD(std::atomic<int>& processed_counter, int num_sets);

    static void solveConnectedSet(ConnectedSetTask task, GRBEnv& env,std::atomic<int>& processed_counter);

    static void decomposeIntoConnectedComponents(DecompositionTask task, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::vector<Polygon_wh>& merged_polys,
                            const Localization& rtree1, const Localization& rtree2,
                            std::atomic<int>& processed_counter);

    static Solution combineThreadSolutions(const std::string& data_name, int num_batches, int num_polys1, int num_polys2, std::vector<int>& set_ids,
                                        std::vector<std::pair<int,int>>& set_sizes,
                                        std::vector<std::vector<double>>& execution_times);

    static void cleanup_temp_folder(const std::string& folder_path);

private:
    //helper function to write thread solution to a file
    static void write_thread_results_to_file(const std::string& filename,
                                  const std::vector<int>& set_ids,
                                  const std::vector<std::pair<int,int>>& set_sizes,
                                  const std::vector<Solution>& sols,
                                  const std::vector<std::vector<double>>& execution_times);

    static void read_thread_results_from_file(const std::string& filename,
                                   std::vector<int>& set_ids,
                                   std::vector<std::pair<int,int>>& set_sizes,
                                   std::vector<Solution>& sols,
                                   std::vector<std::vector<double>>& execution_times);
};



#endif //THREADING_H
