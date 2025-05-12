#include "../include/command_line_parser.h"

#include <iostream>
#include <sstream>

void print_help() {
    std::cout << "Usage: polygonmatching [-d dataset_name] [-l lambda]\n"
              << "Options:\n"
              << "  -d dataset_name      Provide the name of the dataset. There should be two Shapefiles\n"
              << "                       of SINGLE(!) Polygons with appendices _1 and _2.\n"
              << "                       The code expects it to be located in 'input/dataset_name'.\n"
              << "  -l lambda            Specify the parameter for lambda.\n"
              << "  XOR -lr lambda_range Specify multiple values for lambda as START_VALUE STEP_SIZE END_VALUE.\n"
              << "  -o                   (optional, default is jaccard) Specify the objective used for the matching (jaccard,jaccard-hd with two respective weights)."
              << "  -s                   (optional, default is off) Activates making use of the Lemmas of optimal solutions w.r.t. Jaccard Index as Preprocessing.\n"
              << "  -t threads           (optional, default 1) Specify the number of threads that will be used for computations on connected components.\n"
              << "  -r tRee_mode         (optional, default informed) Set the mode describing, how the trees are built (informed,kruskal).\n"
              << "  -m solution_mode     (optional, default is OPTIMAL) Specify, if the optimal solution (opt) or an approximation (3approx) should be computed.\n"
              << "  -e log_name          (optional, default is 'log') Specify the name of the log file with information about the solution that will be written.\n"
              << "  -h, --help           Display this help message.\n";
}

CommandLineOptions parse_command_line(int argc, char* argv[]) {
    CommandLineOptions options;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            exit(0);
        } else if (arg == "-d" && i + 1 < argc) {
            options.dataset_name = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            options.lambdas.push_back(atof(argv[++i]));
        } else if (arg == "-lr" && i + 3 < argc) {
            //collect all lambdas
            double l_min = atof(argv[++i]);
            double step_size = atof(argv[++i]);
            double l_max = atof(argv[++i]);
            //make up for rounding errors
            for (double l = l_min; l <= l_max+0.01; l += step_size) {
                options.lambdas.push_back(l);
            }
        } else if (arg == "-t" & i + 1 < argc) {
            options.num_threads = atoi(argv[++i]);
        } else if (arg == "-e" & i + 1 < argc) {
            options.log_name = argv[++i];
        } else if (arg == "-s") {
            options.exploit_opt_props = true;
        } else if (arg == "-m" & i + 1 < argc) {
            std::string mode_str = argv[++i];
            //transform to lower case for robustness
            std::transform(mode_str.begin(),mode_str.end(), mode_str.begin(),
                           [](unsigned char c){return std::tolower(c);});
            if(mode_str == "opt") {
                options.mode = SOLUTION_MODE::OPT;
            }
            else if(mode_str == "3approx") {
                options.mode = SOLUTION_MODE::CANZAR3APPROX;
            }
            else {
                throw std::invalid_argument("Unknown mode: " + mode_str + "! Use opt, 2approx, 3approx.");
            }
        } else if (arg == "-o" & i+1 < argc) {
            std::string obj_str = argv[++i];
            //transform to lower case for robustness
            std::transform(obj_str.begin(),obj_str.end(), obj_str.begin(),
                           [](unsigned char c){return std::tolower(c);});
            if (obj_str == "jaccard") {
                options.objective = OBJECTIVE::JACCARD;

            } else if (obj_str == "jaccard-hd") {
                options.objective = OBJECTIVE::JACCARD_HAUSDORFF;
                if (i + 2 < argc) {
                    double weight_jac = atof(argv[++i]);
                    double weight_hd = atof(argv[++i]);
                    options.objective_weights = std::make_pair(weight_jac, weight_hd);
                } else {
                    throw std::invalid_argument("No weights provided for multicriterial objective!");
                }
            }
            else {
                throw std::invalid_argument("Unknown objective: " + obj_str + "! Use jaccard, jaccard-hd (weight_jac weight_hd)");
            }

        } else if (arg == "-r" & i+1 < argc) {
            //transform to lower case for robustness
            std::string mode_str = argv[++i];
            std::transform(mode_str.begin(),mode_str.end(), mode_str.begin(),
                           [](unsigned char c){return std::tolower(c);});

            if (mode_str == "informed") {
                options.tree_mode = TREE_BUILD_MODE::INFORMED;
            }
            else if (mode_str == "kruskal") {
                options.tree_mode = TREE_BUILD_MODE::KRUSKAL;
            }
            else {
                throw std::invalid_argument("Unknown tree build mode: " + arg + "! Use informed or kruskal.");
            }
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }
    return options;
}