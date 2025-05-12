#ifndef POLYGONMATCHING_COMMAND_LINE_PARSER_H
#define POLYGONMATCHING_COMMAND_LINE_PARSER_H

#include <string>
#include <vector>
#include <optional>
#include <stdexcept>
#include <algorithm>

enum class TREE_BUILD_MODE {
    INFORMED,
    KRUSKAL
};

enum class SOLUTION_MODE {
    OPT,
    CANZAR3APPROX
};

enum class OBJECTIVE {
    JACCARD,
    JACCARD_HAUSDORFF,
};



//struct to store set command line options
struct CommandLineOptions {
    std::string dataset_name;
    std::string log_name = "log";
    std::vector<double> lambdas;
    int num_threads = 1; //set default to 1 thread
    bool exploit_opt_props = false;
    TREE_BUILD_MODE tree_mode = TREE_BUILD_MODE::INFORMED; // mode for how the trees are being built
    SOLUTION_MODE mode = SOLUTION_MODE::CANZAR3APPROX; //mode for setting the algorithm to either compute the optimal tree constraint matching or approximate it
    OBJECTIVE objective = OBJECTIVE::JACCARD;
    // weights for Jaccard-index and Hausdorff-distance, should the combined objective be chosen
    std::pair<double, double> objective_weights = std::make_pair(0.0, 0.0);
    int batch_size = 1000; //number of connected sets to be processed until dumping to binary file (higher -> more RAM)
};

void print_help();
CommandLineOptions parse_command_line(int argc, char* argv[]);

#endif //POLYGONMATCHING_COMMAND_LINE_PARSER_H
