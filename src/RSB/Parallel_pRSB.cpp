#include <iostream>
#include <thread>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>
#include <mutex>
#include "../Graph/Graph.h" // Include your Graph class header here

std::mutex totalCutSizeMutex;
int totalCutSize = 0;

extern int calculateCutSizeBetweenPartitions(Graph& G, const std::vector<bool>& partitionA, const std::vector<bool>& partitionB);
extern std::vector<int> sortIndices(const Eigen::VectorXd& vec);

void partitionIndices(int threadIndex, int totalThreads, int numIndices, std::vector<int>& sortedIndices,
    std::vector<std::vector<bool>>& partitions, int subPartitionSize) {
    int startIndex = threadIndex * subPartitionSize;
    int endIndex = std::min((threadIndex + 1) * subPartitionSize, numIndices);

    for (int j = startIndex; j < endIndex; ++j) {
        partitions[threadIndex][sortedIndices[j]] = true;
    }
}

std::vector<std::vector<bool>> parallelpRSB(Graph& G, int p) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    // Compute Laplacian matrix
    for (int i = 0; i < sizeNodes; i++) {
        for (int j = 0; j < sizeNodes; j++) {
            L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
        Eigen::VectorXd fiedlerVector = eigenvectors.col(1);

        std::vector<int> sortedIndices = sortIndices(fiedlerVector);

        int subPartitionSize = sizeNodes / p;
        std::vector<std::vector<bool>> partitions(p, std::vector<bool>(sizeNodes, false));
        std::vector<double> partitionWeights(p, 0.0);

        std::vector<std::thread> partitionThreads(p);
        for (int i = 0; i < p; ++i) {
            partitionThreads[i] = std::thread(partitionIndices, i, p, fiedlerVector.size(), std::ref(sortedIndices),
                std::ref(partitions), subPartitionSize);
        }

        for (auto& thread : partitionThreads) {
            thread.join();
        }

        // Proceed with weight calculations and cut size calculations as before
        for (int i = 0; i < p; ++i) {
            // Calculate partition weight for partitions[i]
            // ...

            for (int j = i + 1; j < p; ++j) {
                // Calculate cut size between partitions[i] and partitions[j]
                // ...
            }
        }

        std::cout << "Cut size pRSB: " << totalCutSize << std::endl;

        return partitions;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}