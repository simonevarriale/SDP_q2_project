#include "../Graph/Graph.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/LU>
#include <cmath>
#include <thread>
#include <mutex>

std::mutex mutex;
std::mutex partitionMutex;
std::mutex weightAMutex;
std::mutex weightBMutex;
bool parallel_hasImproved = true;
int numThreads = 4;

extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
extern Eigen::VectorXd rqi(Eigen::VectorXd fv, Eigen::MatrixXd L, int sizeNodes);
extern Graph coarsening(Graph& G);

double parallel_computeMedian(const Eigen::VectorXd& vector) {
    Eigen::VectorXd sortedVector = vector;
    std::sort(sortedVector.data(), sortedVector.data() + sortedVector.size());
    int size = sortedVector.size();
    double medianValue;

    if (size % 2 == 0) {
        // If the vector size is even, take the average of the middle two values
        int mid = size / 2;
        std::mutex mutex;
        std::vector<double> middleValues(2);

        // Compute the middle two values in parallel
        std::vector<std::thread> threads;
        for (int i = 0; i < 2; i++) {
            threads.emplace_back([&sortedVector, &middleValues, mid, i, &mutex]() {
                std::lock_guard<std::mutex> lock(mutex);
                middleValues[i] = sortedVector(mid - 1 + i);
                });
        }

        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }

        medianValue = 0.5 * (middleValues[0] + middleValues[1]);
    }
    else {
        // If the vector size is odd, return the middle value
        medianValue = sortedVector(size / 2);
    }

    return medianValue;
}

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

std::vector<bool> Parallel_RSB(Graph& G, int p) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);
    int numThreads = p; // Number of threads
    std::vector<std::thread> threads;

    G.computeMatrixDegree();
    // G.computeAdjacencyMatrix();
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

        std::cout << "Partition Balance Factor Parallel RSB: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
        std::cout << "Cut size Parallel RSB: " << calculateCutSize(G, partition) << std::endl;

        return partition;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}

Eigen::VectorXd parallel_interpolate(Eigen::VectorXd fv1, Eigen::MatrixXd L, int sizeNodes) {
    // std::cout << "parallel_interpolate" << std::endl;
    Eigen::VectorXd fv(sizeNodes);

    int sum = 0;
    int num = 0;

    for (int i = 0; i < fv1.size(); i++) {
        fv[i] = fv1[i];
    }

    std::vector<std::thread> threads;
    std::vector<int> indices;

    int numThreads = 4;
    int chunkSize = (sizeNodes - fv1.size() + numThreads - 1) / numThreads; // calculate the chunk size
    int start = fv1.size(); // set the starting index

    for (int t = 0; t < numThreads; t++) {
        int end = std::min(start + chunkSize, sizeNodes); // calculate the ending index for this chunk
        threads.emplace_back([&fv, &L, &sum, &num, start, end]() {
            for (int i = start; i < end; i++) {
                for (int j = 0; j < L.cols(); j++) {
                    if (L(i, j) != 0) {
                        std::lock_guard<std::mutex> lock(mutex);
                        sum += fv[j];
                        num++;
                    }
                }
                std::lock_guard<std::mutex> lock(mutex);
                fv[i] = sum / num;
                sum = 0;
                num = 0;
            }
            });
        start = end;
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    return fv;
}

Eigen::VectorXd parallel_fiedler(Graph& G) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);
    Eigen::VectorXd fv;
    Eigen::VectorXd fv1;

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    std::vector<std::thread> threads;

    // compute Laplacian matrix in parallel
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

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    if (sizeNodes > 50 && parallel_hasImproved) {
        // std::cout << "Coarsening" << std::endl;
        Graph G1 = coarsening(G);
        if (G1.num_of_nodes() == G.num_of_nodes() - 1)
            parallel_hasImproved = false;
        // std::cout << "Fine Coarsening con " << G1.num_of_nodes() << " nodes." << std::endl;
        fv1 = parallel_fiedler(G1);
        fv = parallel_interpolate(fv1, L, sizeNodes);
        fv = rqi(fv, L, sizeNodes);
        // std::cout << "Fine parallel_rqi" << std::endl;
    }
    else {
        // std::cout << "Laplacian matrix:" << std::endl;
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

        if (eigenSolver.info() == Eigen::Success) {
            Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
            fv = eigenvectors.col(1);
        }
    }

    return fv;
}

std::vector<bool> Parallel_MLRSB(Graph& G, int p) {
    int sizeNodes = G.num_of_nodes();
    Eigen::VectorXd fiedlerV;

    fiedlerV = parallel_fiedler(G);

    double medianValue = parallel_computeMedian(fiedlerV);
    std::vector<bool> partition(G.num_of_nodes());

    std::vector<std::thread> threads;
    std::mutex mutex;

    int chunkSize = G.num_of_nodes() / numThreads;
    for (int t = 0; t < numThreads; t++) {
        int start = t * chunkSize;
        int end = (t == numThreads - 1) ? G.num_of_nodes() : (t + 1) * chunkSize;
        threads.emplace_back([&partition, &fiedlerV, &medianValue, start, end, &mutex]() {
            for (int i = start; i < end; ++i) {
                if (fiedlerV(i) <= medianValue) {
                    std::lock_guard<std::mutex> lock(mutex);
                    partition[i] = 0; // Assign node i to partition 0
                }
                else {
                    std::lock_guard<std::mutex> lock(mutex);
                    partition[i] = 1; // Assign node i to partition 1
                }
            }
            });
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
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

    std::cout << "Partition Balance Factor Parallel MLRSB: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
    std::cout << "Cut size Parallel MLRSB: " << calculateCutSize(G, partition) << std::endl;

    return partition;
}