#ifndef MLRSB_H
#define MLRSB_H

#include "../Graph/Graph.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/LU>
#include <cmath>

Eigen::VectorXd rqi(Eigen::VectorXd fv, Eigen::MatrixXd L, int sizeNodes);
std::vector<bool> MLRSB(Graph& G);
std::vector<std::vector<bool>> pMLRSB(Graph& G, int p);

#endif
