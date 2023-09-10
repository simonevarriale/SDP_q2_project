#ifndef UTILS_H
#define UTILS_H

#include "../Graph/Graph.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>
#include <cmath>
#include <thread>
#include <set>
#include <utility>

struct PartitionData {
    std::vector<std::vector<bool>> partitions;
    std::vector<double> executionTimes;
    std::vector<double> balanceFactors;
    double averageBalanceFactor;
    std::vector<int> cutSizes;
    int cutSizePartitions;
    double averageCutSize;
    int totalEdgesWeight;
    int totalNodesWeight;
    struct rusage usage;
    std::string fileName;
    double cpu_percentage;
};

double computeMedian(const Eigen::VectorXd& vector);
double parallel_computeMedian(const Eigen::VectorXd& vector);
double calculateBalanceFactor(Graph& graph, const std::vector<bool>& partitionA);
double calculateBalanceFactorPartitions(Graph& G, const std::vector<std::vector<bool>>& partitions);
double calculateAverageBalanceFactor(const std::vector<double>& balanceFactors);
int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
int calculateCutSizePartitions(Graph& G, const std::vector<std::vector<bool>>& partitions);
double calculateAverageCutSize(const std::vector<int>& cutSizes);
std::vector<int> sortIndices(const Eigen::VectorXd& vec);
void read_input(const std::string& filename, Graph* G);
Graph coarsening(Graph& G);
void savePartitionDataToFile(const PartitionData& partitionData);

#endif