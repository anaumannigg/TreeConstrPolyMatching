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



//struct to store set command line options
struct CommandLineOptions {
    std::string dataset_name;
    std::string log_name = "log";
    std::vector<double> lambdas;
    int num_threads = 1; //set default to 1 thread
    bool exploit_opt_props = false;
    TREE_BUILD_MODE tree_mode = TREE_BUILD_MODE::INFORMED; // mode for how the trees are being built
    SOLUTION_MODE mode = SOLUTION_MODE::CANZAR3APPROX; //mode for setting the algorithm to either compute the optimal tree constraint matching or approximate it

};

void print_help();
CommandLineOptions parse_command_line(int argc, char* argv[]);

#endif //POLYGONMATCHING_COMMAND_LINE_PARSER_H
