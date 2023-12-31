#include "RSB.h"

std::vector<bool> RSB(Graph& G) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    // compute Laplacian matrix
    for (int i = 0; i < sizeNodes; i++) {
        for (int j = 0; j < sizeNodes; j++) {
            L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real(); // ascending order
        Eigen::VectorXd fiedlerVector = eigenvectors.col(1);
        double medianValue = computeMedian(fiedlerVector);
        std::vector<bool> partition(L.rows());
        for (int i = 0; i < L.rows(); ++i) {
            if (fiedlerVector(i) <= medianValue) {
                partition[i] = 0; // Assign node i to partition 0
            }
            else {
                partition[i] = 1; // Assign node i to partition 1
            }
        }

        double weightA = 0.0;
        double weightB = 0.0;
        for (int i = 0; i < G.num_of_nodes(); i++) {
            if (partition[i]) {
                weightA += G.getNodeWeight(i);
            }
            else {
                weightB += G.getNodeWeight(i);
            }
        }

        return partition;

    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}

std::vector<std::vector<bool>> pRSB(Graph& G, int p) {
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

        // Sort fiedlerVector to get partitioning indices
        std::vector<int> sortedIndices = sortIndices(fiedlerVector);

        // Calculate the size of each sub-partition
        int subPartitionSize = sizeNodes / p;
        std::vector<std::vector<bool>> partitions(p, std::vector<bool>(sizeNodes, false));

        for (int i = 0; i < p; ++i) {
            int startIndex = i * subPartitionSize;
            int endIndex = (i == p - 1) ? sizeNodes : (i + 1) * subPartitionSize;

            for (int j = startIndex; j < endIndex; ++j) {
                partitions[i][sortedIndices[j]] = true;
            }
        }

        std::vector<double> partitionWeights(p, 0.0);
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < G.num_of_nodes(); ++j) {
                if (partitions[i][j]) {
                    partitionWeights[i] += G.getNodeWeight(j);
                }
            }
        }

        return partitions;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}
