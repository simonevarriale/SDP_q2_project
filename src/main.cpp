#include "./Graph/Graph.h"
#include "utils/utils.h"
#include "./RSB/MLRSB.h"
#include "./RSB/Parallel_MLRSB.h"
#include "./RSB/RSB.h"
#include "./RSB/Parallel_RSB.h"
#include <iostream>
#include <vector>

#include <chrono>

extern std::vector<bool> kernighanLin(Graph& graph, std::vector<bool> partitionA = {});
extern std::vector<bool> multilevel_KL(Graph& graph);
extern std::vector<bool> fiducciaMattheyses2(Graph& graph, int maxIterations, std::vector<bool> partitionA = {});
extern std::vector<bool> kernighanLin1(Graph& G, std::vector<bool> partition);

std::vector<double> executionTimes;

std::vector<std::vector<bool>> algorithmsRunner(Graph G, std::string chosenAlgorithm, int numPartitions, int numThreads) {
    std::vector<std::vector<bool>> partitions;
    partitions.resize(numPartitions);
    int i = 0;

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
        partitions[0] = kernighanLin1(G,partitions[0]);
    }
    else if (chosenAlgorithm == "FM") {
        partitions[0] = fiducciaMattheyses2(G, 10);
    }
    else {
        std::cout << "Invalid algorithm name" << std::endl;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTimeAlgorithm = endTime - startTime;
    std::cout << "Execution time " << chosenAlgorithm << ": " << execTimeAlgorithm.count() << " seconds" << std::endl;
    executionTimes.push_back(execTimeAlgorithm.count());

    // if (partitions.size() > 1) {
    //     if (partitions[1].size() > 0) {
    //         std::cout << "Cut size between partitions: " << calculateCutSizePartitions(G, partitions) << std::endl;
    //     }
    // }

    return partitions;
}

int main(int argc, char** argv) {
    Graph G;

    int i = 0, numPartitions = 2, numThreads = 2;
    std::string algorithmName = "KL";
    std::string inputGraphFile = "./data/graph_75_150.txt";
    std::vector<double> balanceFactors;
    std::vector<int> cutSizes;

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

    /////////////////////////////////////////////////////GRAPH READ//////////////////////////////////////////////
    auto startTime = std::chrono::high_resolution_clock::now();
    read_input(inputGraphFile, &G);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execTimeRead = endTime - startTime;
    std::cout << "Execution time graph reading: " << execTimeRead.count() << " seconds" << std::endl << std::endl;
    executionTimes.push_back(execTimeRead.count());
    // std::cout << "Num of nodes: " << G.num_of_nodes() << " and Num of Edges: " << G.num_of_edges() << std::endl;
    // G1 = G;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto partitions = algorithmsRunner(G, algorithmName, numPartitions, numThreads);

    // Print the partitioning result
    std::cout << std::endl << "Partitioning result:" << std::endl;
    for (auto v : partitions) {
        if (v.size() > 0) {
            i++;
            if (partitions.size() > 1) {
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
            balanceFactors.push_back(calculateBalanceFactor(G, v));
            cutSizes.push_back(calculateCutSize(G, v));
        }
    }

    savePartitionDataToFile(partitions, executionTimes, balanceFactors, cutSizes, inputGraphFile + "_" + algorithmName + "_" + std::to_string(numPartitions) + "_" + std::to_string(numThreads) + "_results.txt");

    return 0;
}