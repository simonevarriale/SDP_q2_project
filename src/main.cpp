#include "./Graph/Graph.h"
#include "utils/utils.h"
#include "./RSB/MLRSB.h"
#include "./RSB/Parallel_MLRSB.h"
#include "./RSB/RSB.h"
#include "./RSB/Parallel_RSB.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <sys/resource.h>

extern std::vector<bool> kernighanLin(Graph& graph, std::vector<bool> partitionA = {});
extern std::vector<bool> multilevel_KL(Graph& graph);
extern std::vector<bool> fiducciaMattheyses2(Graph& graph, int maxIterations, std::vector<bool> partitionA = {});
extern std::vector<bool> kernighanLin1(Graph& G, std::vector<bool> partition);

PartitionData partitionData;

std::vector<std::vector<bool>> algorithmsRunner(Graph G, std::string chosenAlgorithm, int numPartitions, int numThreads) {
    std::vector<std::vector<bool>> partitions;
    partitions.resize(numPartitions);
    struct rusage usage1;
    int i = 0;
    getrusage(RUSAGE_SELF, &usage1);
    auto startTime = std::chrono::high_resolution_clock::now();
    if (chosenAlgorithm == "RSB") {
        partitions[0] = RSB(G, numPartitions);
    }
    else if (chosenAlgorithm == "pRSB") {
        partitions = pRSB(G, numPartitions);
    }
    else if (chosenAlgorithm == "MLRSB") {
        partitions[0] = MLRSB(G);
    }
    else if (chosenAlgorithm == "pMLRSB") {
        partitions = pMLRSB(G, numPartitions);
    }
    else if (chosenAlgorithm == "Parallel_RSB") {
        partitions[0] = Parallel_RSB(G, numThreads);
    }
    else if (chosenAlgorithm == "Parallel_pRSB") {
        partitions = Parallel_pRSB(G, numPartitions);
    }
    else if (chosenAlgorithm == "Parallel_MLRSB") {
        partitions[0] = Parallel_MLRSB(G, numPartitions);
    }
    else if (chosenAlgorithm == "Parallel_pMLRSB") {
        partitions = Parallel_pMLRSB(G, numPartitions, numThreads);
    }
    else if (chosenAlgorithm == "KL") {
        partitions[0] = kernighanLin1(G, partitions[0]);
    }
    else if (chosenAlgorithm == "FM") {
        partitions[0] = fiducciaMattheyses2(G, 10);
    }
    else {
        std::cout << "Invalid algorithm name" << std::endl;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    getrusage(RUSAGE_SELF, &partitionData.usage);
    std::chrono::duration<double> execTimeAlgorithm = endTime - startTime;
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / 1000000.0;
    double cpu = (partitionData.usage.ru_utime.tv_sec - usage1.ru_utime.tv_sec) + (partitionData.usage.ru_utime.tv_usec - usage1.ru_utime.tv_usec) / 1000000.0;
    partitionData.cpu_percentage = 100.0 * cpu / elapsed;
    std::cout << "Execution time " << chosenAlgorithm << ": " << execTimeAlgorithm.count() << " seconds" << std::endl;
    std::cout << "CPU percentage used: " << partitionData.cpu_percentage << "%" << std::endl;
    std::cout << "Memory usage: " << partitionData.usage.ru_maxrss / (1024.0 * 1024.0) << " GBs" << std::endl;
    partitionData.executionTimes.push_back(execTimeAlgorithm.count());

    return partitions;
}

int main(int argc, char** argv) {
    Graph G;

    int i = 0, numPartitions = 1, numThreads = 1;
    std::string algorithmName = "MLRSB", inputGraphFile = "./data/graph_500_1024.txt";

    if (argc > 1) {
        inputGraphFile = argv[1];
    }
    if (argc > 2) {
        algorithmName = argv[2];
    }
    if (argc > 3) {
        numPartitions = std::atoi(argv[3]);
    }
    if (argc > 4) {
        numThreads = std::atoi(argv[4]);
    }

    // Read the graph from the input file
    auto startTime = std::chrono::high_resolution_clock::now();
    read_input(inputGraphFile, &G);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTimeRead = endTime - startTime;
    partitionData.totalEdgesWeight = G.getTotalEdgesWeight();
    std::cout << "Total edges weight: " << partitionData.totalEdgesWeight << std::endl;
    std::cout << "Execution time graph reading: " << execTimeRead.count() << " seconds" << std::endl << std::endl;
    partitionData.executionTimes.push_back(execTimeRead.count());

    partitionData.partitions = algorithmsRunner(G, algorithmName, numPartitions, numThreads);

    // Print the partitioning result
    std::cout << std::endl << "Partitioning result on graph '" << inputGraphFile << "'" << std::endl;
    for (auto v : partitionData.partitions) {
        if (v.size() > 0) {
            i++;
            if (partitionData.partitions.size() > 1) {
                std::cout << "Partition " << i << ": ";
            }
            else {
                std::cout << "Partition: ";
            }
            for (auto j : v) {
                std::cout << j << " ";
            }
            std::cout << std::endl;
            std::cout << "Balance Factor: " << calculateBalanceFactor(G, v) << " | Cut Size: " << calculateCutSize(G, v) << std::endl << std::endl;
            partitionData.balanceFactors.push_back(calculateBalanceFactor(G, v));
            partitionData.cutSizes.push_back(calculateCutSize(G, v));
        }
    }
    partitionData.averageBalanceFactor = calculateAverageBalanceFactor(partitionData.balanceFactors);
    partitionData.averageCutSize = calculateAverageCutSize(partitionData.cutSizes);
    // partitionData.cutSizePartitions = calculateCutSizePartitions(G, partitionData.partitions);
    partitionData.fileName = "./results/graph_" + std::to_string(G.num_of_nodes()) + "_" + std::to_string(G.num_of_edges()) + "_" + algorithmName + "_" + std::to_string(numPartitions) + "_" + std::to_string(numThreads) + ".txt";

    savePartitionDataToFile(partitionData);

    return 0;
}