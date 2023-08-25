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
#include <condition_variable>
#include <queue>

class ThreadPool {
public:
    ThreadPool(int numThreads) : stop(false) {
        for (int i = 0; i < numThreads; i++) {
            threads.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(mutex);
                        condition.wait(lock, [this]() { return stop || !tasks.empty(); });
                        if (stop && tasks.empty()) {
                            return;
                        }
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
                });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(mutex);
            stop = true;
        }
        condition.notify_all();
        for (auto& thread : threads) {
            thread.join();
        }
    }

    template<class F, class... Args>
    void enqueue(F&& f, Args&&... args) {
        {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.emplace([=]() { std::forward<F>(f)(std::forward<Args>(args)...); });
        }
        condition.notify_one();
    }

private:
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable condition;
    bool stop;
};

std::mutex mutex;
bool parallel_hasImproved = true;

extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
extern Graph coarsening(Graph& G);

std::mutex partitionMutex;
std::mutex weightAMutex;
std::mutex weightBMutex;

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
    G.computeAdjacencyMatrix();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    if (sizeNodes < std::thread::hardware_concurrency()) {
        numThreads = sizeNodes;
    }
    else {
        numThreads = std::thread::hardware_concurrency();
    }
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

        std::cout << "Partition Balance Factor: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
        std::cout << "Cut size RSB: " << calculateCutSize(G, partition) << std::endl;

        return partition;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}

Eigen::VectorXd parallel_rqi(Eigen::VectorXd fv, Eigen::MatrixXd L, int sizeNodes) {
    std::cout << "RQI" << std::endl;
    float theta = fv.transpose() * L * fv;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sizeNodes, sizeNodes);
    double p;
    std::mutex mutex;

    int maxIterations = 1000;
    int iteration = 0;

    do {
        std::vector<std::thread> threads;
        std::vector<double> results(sizeNodes);

        // Compute x for each i in parallel
        int numThreads;
        if (sizeNodes < std::thread::hardware_concurrency()) {
            numThreads = sizeNodes;
        }
        else {
            numThreads = std::thread::hardware_concurrency();
        }
        int chunkSize = sizeNodes / numThreads;
        for (int t = 0; t < numThreads; t++) {
            int start = t * chunkSize;
            int end = (t == numThreads - 1) ? sizeNodes : (t + 1) * chunkSize;
            threads.emplace_back([&fv, &L, &theta, &I, &results, start, end]() {
                for (int i = start; i < end; i++) {
                    Eigen::VectorXd x = (L - theta * I).lu().solve(fv);
                    results[i] = x[i];
                }
                });
        }

        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }

        // Compute fv and theta using the computed x values
        double norm = 0;
        for (int i = 0; i < sizeNodes; i++) {
            fv[i] = results[i];
            norm += fv[i] * fv[i];
        }
        norm = sqrt(norm);
        fv /= norm;
        theta = fv.transpose() * L * fv;
        p = ((L * fv).transpose() * (L * fv) - theta * theta);
        p = sqrt(p);
        std::cout << "p: " << p << std::endl;

        iteration++;
    } while (p < 0.0000001 && iteration < maxIterations);

    return fv;
}

Eigen::VectorXd parallel_interpolate(Eigen::VectorXd fv1, Eigen::MatrixXd L, int sizeNodes) {
    std::cout << "parallel_interpolate" << std::endl;
    Eigen::VectorXd fv(sizeNodes);

    int sum = 0;
    int num = 0;

    for (int i = 0; i < fv1.size(); i++) {
        fv[i] = fv1[i];
    }

    std::vector<std::thread> threads;
    std::vector<int> indices;

    // Create threads to compute fv for each i
    for (int i = fv1.size(); i < sizeNodes; i++) {
        indices.push_back(i);
        threads.emplace_back([&fv, &L, &sum, &num, i]() {
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
            });
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
    // int numThreads = std::min(sizeNodes, std::thread::hardware_concurrency());
    int numThreads;
    if (sizeNodes < std::thread::hardware_concurrency()) {
        numThreads = sizeNodes;
    }
    else {
        numThreads = std::thread::hardware_concurrency();
    }
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
        std::cout << "Coarsening" << std::endl;
        Graph G1 = coarsening(G);
        if (G1.num_of_nodes() == G.num_of_nodes() - 1)
            parallel_hasImproved = false;
        std::cout << "Fine Coarsening con " << G1.num_of_nodes() << " nodes." << std::endl;
        fv1 = parallel_fiedler(G1);
        fv = parallel_interpolate(fv1, L, sizeNodes);
        fv = parallel_rqi(fv, L, sizeNodes);
        std::cout << "Fine parallel_rqi" << std::endl;
    }
    else {
        std::cout << "Laplacian matrix:" << std::endl;
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

        if (eigenSolver.info() == Eigen::Success) {
            Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
            fv = eigenvectors.col(1);
        }
    }

    return fv;
}

// Eigen::VectorXd parallel_rqi(Eigen::VectorXd fv, Eigen::MatrixXd L, int sizeNodes) {
//     std::cout << "RQI" << std::endl;
//     float theta = fv.transpose() * L * fv;
//     Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sizeNodes, sizeNodes);
//     double p;
//     std::mutex mutex;

//     do {
//         std::vector<std::thread> threads;
//         std::vector<double> results(sizeNodes);

//         // Compute x for each i in parallel
//         for (int i = 0; i < sizeNodes; i++) {
//             threads.emplace_back([&fv, &L, &theta, &I, &results, i]() {
//                 Eigen::VectorXd x = (L - theta * I).lu().solve(fv);
//                 results[i] = x[i];
//                 });
//         }

//         // Wait for all threads to finish
//         for (auto& thread : threads) {
//             thread.join();
//         }

//         // Compute fv and theta using the computed x values
//         double norm = 0;
//         for (int i = 0; i < sizeNodes; i++) {
//             fv[i] = results[i];
//             norm += fv[i] * fv[i];
//         }
//         norm = sqrt(norm);
//         fv /= norm;
//         theta = fv.transpose() * L * fv;
//         p = ((L * fv).transpose() * (L * fv) - theta * theta);
//         p = sqrt(p);
//         std::cout << "p: " << p << std::endl;
//     } while (p < 0.0000001);

//     return fv;
// }

// Eigen::VectorXd parallel_fiedler(Graph& G) {
//     int sizeNodes = G.num_of_nodes();
//     Eigen::MatrixXd L(sizeNodes, sizeNodes);
//     Eigen::VectorXd fv;
//     Eigen::VectorXd fv1;

//     G.computeMatrixDegree();
//     auto matDeg = G.getMatDegree();
//     auto matAdj = G.getMatAdj();

//     std::vector<std::thread> threads;

//     // compute Laplacian matrix in parallel
//     for (int i = 0; i < sizeNodes; i++) {
//         threads.emplace_back([&L, &matDeg, &matAdj, i, sizeNodes]() {
//             for (int j = 0; j < sizeNodes; j++) {
//                 L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
//             }
//             });
//     }

//     // Wait for all threads to finish
//     for (auto& thread : threads) {
//         thread.join();
//     }

//     if (sizeNodes > 50 && parallel_hasImproved) {
//         std::cout << "Coarsening" << std::endl;
//         Graph G1 = coarsening(G);
//         if (G1.num_of_nodes() == G.num_of_nodes() - 1)
//             parallel_hasImproved = false;
//         std::cout << "Fine Coarsening con " << G1.num_of_nodes() << " nodes." << std::endl;
//         fv1 = parallel_fiedler(G1);
//         fv = parallel_interpolate(fv1, L, sizeNodes);
//         fv = parallel_rqi(fv, L, sizeNodes);
//         std::cout << "Fine parallel_rqi" << std::endl;
//     }
//     else {
//         std::cout << "Laplacian matrix:" << std::endl;
//         Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

//         if (eigenSolver.info() == Eigen::Success) {
//             Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
//             fv = eigenvectors.col(1);
//         }
//     }

//     return fv;
// }

// std::vector<bool> Parallel_MLRSB(Graph& G, int p) {
//     Eigen::VectorXd fiedlerV;

//     fiedlerV = parallel_fiedler(G);

//     double medianValue = parallel_computeMedian(fiedlerV);
//     std::vector<bool> partition(G.num_of_nodes());
//     std::cout << "Start partitioning" << std::endl;

//     std::vector<std::thread> threads;
//     std::mutex mutex;

//     // Assign nodes to partitions in parallel
//     for (int i = 0; i < G.num_of_nodes(); ++i) {
//         threads.emplace_back([&partition, &fiedlerV, &medianValue, i, &mutex]() {
//             if (fiedlerV(i) <= medianValue) {
//                 std::lock_guard<std::mutex> lock(mutex);
//                 partition[i] = 0; // Assign node i to partition 0
//             }
//             else {
//                 std::lock_guard<std::mutex> lock(mutex);
//                 partition[i] = 1; // Assign node i to partition 1
//             }
//             });
//     }

//     // Wait for all threads to finish
//     for (auto& thread : threads) {
//         thread.join();
//     }

//     std::cout << "End partitioning" << std::endl;

//     double weightA = 0.0;
//     double weightB = 0.0;
//     for (int i = 0; i < G.num_of_nodes(); i++) {
//         if (partition[i]) {
//             weightA += G.getNodeWeight(i);
//         }
//         else {
//             weightB += G.getNodeWeight(i);
//         }
//     }

//     std::cout << "Partition Balance Factor: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
//     std::cout << "Cut size Parallel MLRSB: " << calculateCutSize(G, partition) << std::endl;

//     return partition;
// }

std::vector<bool> Parallel_MLRSB(Graph& G, int p) {
    int sizeNodes = G.num_of_nodes();
    Eigen::VectorXd fiedlerV;

    fiedlerV = parallel_fiedler(G);

    double medianValue = parallel_computeMedian(fiedlerV);
    std::vector<bool> partition(G.num_of_nodes());
    std::cout << "Start partitioning" << std::endl;

    std::vector<std::thread> threads;
    std::mutex mutex;

    // Assign nodes to partitions in parallel
    // int numThreads = std::min(G.num_of_nodes(), std::thread::hardware_concurrency());
    int numThreads;
    if (sizeNodes < std::thread::hardware_concurrency()) {
        numThreads = sizeNodes;
    }
    else {
        numThreads = std::thread::hardware_concurrency();
    }
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

    std::cout << "End partitioning" << std::endl;

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

    std::cout << "Partition Balance Factor: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
    std::cout << "Cut size Parallel MLRSB: " << calculateCutSize(G, partition) << std::endl;

    return partition;
}