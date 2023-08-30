#ifndef RSB_H
#define RSB_H

#include "../utils/utils.h"
#include "../Graph/Graph.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/LU>
#include <cmath>

std::vector<bool> RSB(Graph& G, int p);
std::vector<std::vector<bool>> pRSB(Graph& G, int p);

#endif