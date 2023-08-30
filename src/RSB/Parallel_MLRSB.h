#ifndef PARALLEL_MLRSB_H
#define PARALLEL_MLRSB_H

#include "../Graph/Graph.h"
#include "MLRSB.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <thread>
#include <mutex>

std::vector<bool> Parallel_MLRSB(Graph& G, int numThreads);
std::vector<std::vector<bool>> Parallel_pMLRSB(Graph& G, int p, int numThreads);

#endif