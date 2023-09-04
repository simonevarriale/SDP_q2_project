#include "Parallel_RSB.h"

std::mutex partitionMutex;

void processNodePartition(const Eigen::VectorXd& fiedlerVector, const Graph& G, int start, int end, std::vector<bool>& partition, double& weightA, double& weightB) {
    for (int i = start; i <= end; ++i) {
        if (fiedlerVector(i) <= parallel_computeMedian(fiedlerVector)) {
            std::lock_guard<std::mutex> lock(partitionMutex);
            partition[i] = 0;
            weightA += G.getNodeWeight(i);
        }
        else {
            std::lock_guard<std::mutex> lock(partitionMutex);
            partition[i] = 1;
            weightB += G.getNodeWeight(i);
        }
    }
}

void partitionIndices(int threadIndex, int totalThreads, int numIndices, std::vector<int>& sortedIndices,
    std::vector<std::vector<bool>>& partitions, int subPartitionSize) {
    int startIndex = threadIndex * subPartitionSize;
    int endIndex = std::min((threadIndex + 1) * subPartitionSize, numIndices);

    for (int j = startIndex; j < endIndex; ++j) {
        partitions[threadIndex][sortedIndices[j]] = true;
    }
}

void computeLRow(int row, int sizeNodes, const std::vector<std::vector<int>>& matDeg, const std::vector<std::vector<std::vector<int>>>& matAdj, Eigen::MatrixXd& L) {
    for (int j = 0; j < sizeNodes; j++) {
        L(row, j) = matDeg[row][j] - matAdj[row][j][0] * matAdj[row][j][1];
    }
}

std::vector<bool> Parallel_RSB(Graph& G, int numThreads) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);
    std::vector<std::thread> threads;

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    int chunkSize = sizeNodes / numThreads;
    for (int t = 0; t < numThreads; t++) {
        int start = t * chunkSize;
        int end = (t == numThreads - 1) ? sizeNodes : (t + 1) * chunkSize;
        threads.emplace_back([&L, &matDeg, &matAdj, start, end, sizeNodes]() {
            for (int i = start; i < end; i++) {
                for (int j = 0; j < sizeNodes; j++) {
                    L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
                }
            }
            });
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
        Eigen::VectorXd fiedlerVector = eigenvectors.col(1);

        double medianValue = parallel_computeMedian(fiedlerVector);
        std::vector<bool> partition(L.rows());
        double weightA = 0.0;
        double weightB = 0.0;


        int chunkSize = sizeNodes / numThreads;
        int start = 0;
        for (int t = 0; t < numThreads; t++) {
            int end = (t == numThreads - 1) ? sizeNodes - 1 : start + chunkSize - 1;
            threads.emplace_back(processNodePartition, std::ref(fiedlerVector), std::cref(G), start, end, std::ref(partition), std::ref(weightA), std::ref(weightB));
            start = end + 1;
        }

        for (auto& thread : threads) {
            thread.join();
        }

        return partition;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}

std::vector<std::vector<bool>> Parallel_pRSB(Graph& G, int p) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    std::vector<std::thread> rowThreads(sizeNodes);
    for (int i = 0; i < sizeNodes; ++i) {
        rowThreads[i] = std::thread(computeLRow, i, sizeNodes, std::ref(matDeg), std::ref(matAdj), std::ref(L));
    }

    for (auto& thread : rowThreads) {
        thread.join();
    }

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

        return partitions;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}