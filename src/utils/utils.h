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

double computeMedian(const Eigen::VectorXd& vector);
double parallel_computeMedian(const Eigen::VectorXd& vector);
double calculateBalanceFactor(Graph& graph, const std::vector<bool>& partitionA);
bool isPartitionBalanced(Graph& graph, const std::vector<bool>& partitionA);
void computeInitialGains(Graph& graph, const std::vector<bool>& partitionA, std::vector<int>& gains);
void computeNetGains(Graph& graph, const std::vector<bool>& partitionA, std::vector<int>& netGains);
int calculateNodeGain(Graph& graph, const std::vector<bool>& partitionA, int node, bool moveToPartitionA);
int calculateTotalWeight(Graph& graph);
int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
int calculateCutSizePartitions(Graph& G, const std::vector<std::vector<bool>>& partitions);
std::vector<int> sortIndices(const Eigen::VectorXd& vec);
void read_input(const std::string& filename, Graph* G);
void read_input2(const std::string& filename, Graph* G);
Graph coarsening(Graph& G);
std::unordered_map<int, std::pair<int, int>> coarsenGraph(Graph& G);
std::vector<bool> uncoarsening(Graph G1, std::vector<bool> partition, int graphSize);
std::vector<bool> uncoarsening2(std::unordered_map<int, std::pair<int, int>> coarse, std::vector<bool> partition);
void savePartitionDataToFile(const std::vector<std::vector<bool>>& partitions, const std::vector<double> execTimes, const std::vector<double>& balanceFactors, const std::vector<int>& cutSizes, const std::string& filename);

#endif