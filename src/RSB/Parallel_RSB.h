#ifndef PARALLEL_RSB_H
#define PARALLEL_RSB_H

#include "../Graph/Graph.h"
#include "../utils/utils.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <thread>
#include <mutex>

std::vector<bool> Parallel_RSB(Graph& G, int numThreads);
std::vector<std::vector<bool>> Parallel_pRSB(Graph& G, int p);

#endif