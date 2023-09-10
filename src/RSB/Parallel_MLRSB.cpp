#include "Parallel_MLRSB.h"

std::mutex mutex;
bool parallel_hasImproved = true;

Eigen::VectorXd parallel_interpolate(Eigen::VectorXd fv1, Eigen::MatrixXd L, int sizeNodes, int numThreads) {
    Eigen::VectorXd fv(sizeNodes);

    int sum = 0;
    int num = 0;

    for (int i = 0; i < fv1.size(); i++) {
        fv[i] = fv1[i];
    }

    std::vector<std::thread> threads;
    std::vector<int> indices;

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

Eigen::VectorXd parallel_fiedler(Graph& G, int numThreads) {
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
        fv1 = parallel_fiedler(G1, numThreads);
        fv = parallel_interpolate(fv1, L, sizeNodes, numThreads);
        fv = rqi(fv, L, sizeNodes);
    }
    else {
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

        if (eigenSolver.info() == Eigen::Success) {
            Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
            fv = eigenvectors.col(1);
        }
    }

    return fv;
}

std::vector<bool> Parallel_MLRSB(Graph& G, int numThreads) {
    int sizeNodes = G.num_of_nodes();
    Eigen::VectorXd fiedlerV;

    fiedlerV = parallel_fiedler(G, numThreads);

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

    return partition;
}

std::vector<std::vector<bool>> Parallel_pMLRSB(Graph& G, int p, int numThreads) {
    Eigen::VectorXd fiedlerV = parallel_fiedler(G, numThreads);

    std::vector<double> partitionValues(G.num_of_nodes());

    // Parallelize calculating partitionValues using threads
    std::vector<std::thread> threads;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        threads.emplace_back([i, &partitionValues, &fiedlerV]() {
            partitionValues[i] = fiedlerV(i);
            });
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    std::vector<double> partitionThresholds(p - 1);
    std::vector<std::vector<bool>> partitions(p, std::vector<bool>(G.num_of_nodes(), false));

    // Sort the partition values
    std::sort(partitionValues.begin(), partitionValues.end());

    // Calculate partition thresholds
    for (int i = 1; i < p; ++i) {
        partitionThresholds[i - 1] = partitionValues[i * G.num_of_nodes() / p];
    }

    // Assign nodes to partitions based on thresholds
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        bool assigned = false;
        for (std::size_t j = 0; j < p - 1; ++j) {
            if (fiedlerV(i) <= partitionThresholds[j]) {
                partitions[j][i] = true;
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            partitions[p - 1][i] = true;
        }
    }

    return partitions;
}
