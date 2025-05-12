#ifndef POLYGONMATCHING_LOGGER_H
#define POLYGONMATCHING_LOGGER_H

#include <string>
#include <fstream>
#include <vector>

class Logger {
public:
    // Constructor to initialize the file path
    Logger(const std::string& filePath);

    // Config setter
    void setTreeMode(std::string treeMode);
    void setSolMode(std::string solMode);
    void setObjectiveMode(std::string objMode);
    void setExploitMode(bool exploit_opt_props);

    // Input setters
    void setInputPolygons(int osm, int atkis);
    void setInstanceName(const std::string& name);
    void setLambda(double lambda);
    void setConnectedComponents(int cc);

    // Runtime setters (expects decomposition, setup, exploit, trees, solving)
    void setTimings(std::vector<double> timings);

    // Memory setter
    void setMemoryUsage(int mem_usage_kb);

    // Output setters
    void setTotalNumMatches(int arcs);
    void setMatchDistribution(const std::vector<int>& dist);
    void setObjectiveDistribution(const std::vector<double>& dist);
    void setObjective(double objective);

    // Logging and reset
    void log();
    void reset();

private:
    // File path for the CSV
    std::string filePath;

    //config information
    std::string treeMode;
    std::string solMode;
    std::string objectiveMode;
    bool exploited_opt_props;

    // Input information
    int inputPolygons1;
    int inputPolygons2;
    int connectedComponents;
    std::string instanceName;
    double lambda;

    // Runtime information
    double decompTime;
    double solvingSetupTime;
    double solvingExploitTime;
    double solvingTreeTime;
    double solvingOptTime;
    //stores the actual execution time after decomposition (important to measure multithreading improvement)
    double actualExecTime;
    //store memory usage in kilobyte
    int memoryUsageKB;

    // Output information
    int totalNumMatches;
    int totalNumUnmatched;
    int totalNum1to1;
    int totalNum1toN;
    int totalNumMtoN;
    double obj1to1;
    double obj1toN;
    double objMtoN;
    double objective;

    // Helper to initialize the file with headers if empty
    void initializeFile();
};


#endif //POLYGONMATCHING_LOGGER_H
