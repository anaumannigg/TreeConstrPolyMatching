#include "../include/logger.h"

#include <fstream>
#include <iomanip>
#include <assert.h>

#include <iostream>
#include <filesystem>

// Constructor
Logger::Logger(const std::string& filePath) :
        filePath(filePath), inputPolygons1(0), inputPolygons2(0), connectedComponents(0), lambda(0.0),
        decompTime(0.0),solvingSetupTime(0.0),solvingExploitTime(0.0),solvingOptTime(0.0), actualExecTime(0.0), memoryUsageKB(0),
        totalNumMatches(0), totalNum1to1(0), totalNum1toN(0), totalNumMtoN(0), objective(0.0)
        {
    initializeFile();
}

// Private method to add headers to the file if it is empty
void Logger::initializeFile() {
    // Check if the directory exists, and create it if it doesn't
    std::filesystem::path fileDir = std::filesystem::path(filePath).parent_path();

    if (!std::filesystem::exists(fileDir)) {
        if (!std::filesystem::create_directories(fileDir)) {
            std::cerr << "Error creating directory for logging: " << fileDir << std::endl;
            return;
        }
    }

    if (!std::filesystem::exists(filePath) || std::filesystem::file_size(filePath) == 0) {
        std::ofstream outfile(filePath, std::ios::app);
        if (!outfile) {  // Check if the file was successfully opened
            std::cerr << "Error initializing log file: " << filePath << std::endl;
        }
        outfile << "InstanceName,TreeMode,SolMode,ObjectiveMode,OptPropExploit,polygons1,polygons2,ConnectedComponents,Lambda,"
                << "DecompTime,SetupTime,ExploitTime,TreeBuildTime,OptTime,ExecTime,PeakMemoryUsed,TotalMatches,"
                << "TotalUnmatched,Total1to1,Total1toN,TotalMtoN,Obj1to1,Obj1toN,ObjMtoN,Objective"
                << std::endl;
        outfile.close();
    }
}

void Logger::setTreeMode(std::string treeMode) {
    this->treeMode = treeMode;
}

void Logger::setSolMode(std::string solMode) {
    this->solMode = solMode;
}

void Logger::setObjectiveMode(std::string objMode) {
  this->objectiveMode = objMode;
}

void Logger::setExploitMode(bool exploit_opt_props) {
    this->exploited_opt_props = exploit_opt_props;
}

// Input setters
void Logger::setInputPolygons(int polys1, int polys2) {
    inputPolygons1 = polys1;
    inputPolygons2 = polys2;
}

void Logger::setInstanceName(const std::string& name) {
    instanceName = name;
}

void Logger::setConnectedComponents(int cc) {connectedComponents=cc;}

void Logger::setLambda(double lambdaValue) {
    lambda = lambdaValue;
}

// Runtime setters
void Logger::setTimings(std::vector<double> timings) {
    if (timings.size() != 6) {
        return;
    }
    decompTime = timings[0];
    solvingSetupTime = timings[1];
    solvingExploitTime = timings[2];
    solvingTreeTime = timings[3];
    solvingOptTime = timings[4];
    actualExecTime = timings[5];
}

void Logger::setMemoryUsage(int mem_usage_kb) {
    this->memoryUsageKB = mem_usage_kb;
}

// Output setters
void Logger::setTotalNumMatches(int arcs) {
    totalNumMatches = arcs;
}

void Logger::setMatchDistribution(const std::vector<int>& dist) {
    assert(dist.size() == 4 && "cannot set match distribution with invalid amount of fields!");
    this->totalNumUnmatched = dist[0];
    this->totalNum1to1 = dist[1];
    this->totalNum1toN = dist[2];
    this->totalNumMtoN = dist[3];
}

void Logger::setObjectiveDistribution(const std::vector<double>& dist) {
    assert(dist.size() == 3 && "cannot set objective distribution with invalid amount of fields!");
    this->obj1to1  = dist[0];
    this->obj1toN  = dist[1];
    this->objMtoN  = dist[2];
}

void Logger::setObjective(double objectiveValue) {
    objective = objectiveValue;
}

// Log method to write data to CSV
void Logger::log() {
    std::ofstream outfile(filePath, std::ios::app);
    outfile << instanceName << ","
            << treeMode << ","
            << solMode << ","
            << objectiveMode << ","
            << exploited_opt_props << ","
            << inputPolygons1 << ","
            << inputPolygons2 << ","
            << connectedComponents << ","
            << std::fixed << std::setprecision(1)
            << lambda << ","
            << std::fixed << std::setprecision(3)
            << decompTime << ","
            << solvingSetupTime << ","
            << solvingExploitTime << ","
            << solvingTreeTime << ","
            << solvingOptTime << ","
            << actualExecTime << ","
            << memoryUsageKB << ","
            << totalNumMatches << ","
            << totalNumUnmatched << ","
            << totalNum1to1 << ","
            << totalNum1toN << ","
            << totalNumMtoN << ","
            << std::fixed << std::setprecision(3)
            << obj1to1 << ","
            << obj1toN << ","
            << objMtoN << ","
            << objective
         << std::endl;
    outfile.close();
}

// Reset method to clear data
void Logger::reset() {
    inputPolygons1 = 0;
    inputPolygons2 = 0;
    exploited_opt_props = false;
    instanceName.clear();
    treeMode.clear();
    solMode.clear();
    objectiveMode.clear();
    lambda = 0.0;
    decompTime = 0.0;
    solvingSetupTime = 0.0;
    solvingExploitTime = 0.0;
    solvingTreeTime = 0.0;
    solvingOptTime = 0.0;
    actualExecTime = 0.0;
    memoryUsageKB = 0;
    totalNumMatches = 0;
    totalNumUnmatched = 0;
    totalNum1to1 = 0;
    totalNum1toN = 0;
    totalNumMtoN = 0;
    obj1to1 = 0.0;
    obj1toN = 0.0;
    objMtoN = 0.0;
    objective = 0.0;
}
