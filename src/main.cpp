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
        partitions[0] = kernighanLin(G);
    }
    else if (chosenAlgorithm == "FM") {
        partitions[0] = fiducciaMattheyses2(G, 10);
    }
    else {
        std::cout << "Invalid algorithm name" << std::endl;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;
    std::cout << "Execution time " << chosenAlgorithm << ": " << duration.count() << " seconds" << std::endl;
    if (partitions.size() > 1) {
        if (partitions[1].size() > 0) {
            std::cout << "Cut size: " << calculateCutSizePartitions(G, partitions) << std::endl;
        }
        else {
            std::cout << "Cut size: " << calculateCutSize(G, partitions[0]) << std::endl;
        }
    }
    else {
        std::cout << "Cut size: " << calculateCutSize(G, partitions[0]) << std::endl;
    }

    return partitions;
}

int main(int argc, char** argv) {
    Graph G;
    // Graph G1;
    int i = 0, numPartitions = 2, numThreads = 2;
    std::string algorithmName = "MLRSB";
    std::string inputGraphFile = "./data/test_graph.txt";
    // auto startTime = std::chrono::high_resolution_clock::now();
    // auto endTime = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = endTime - startTime;
    // std::vector<bool> partition = {};
    // std::string file2 = "./data/test_graph.txt";
    // std::string file3 = "./data/connected_graph.txt";
    // std::string file4 = "./data/big_graph.txt";
    // std::string file5 = "./data/3elt.graph";
    // std::string file6 = "./data/simple.graph";
    // std::string file7 = "./data/add20.graph";
    // std::string file8 = "./data/graph_20230828_171909.txt";
    // std::string file9 = "./data/graph_20230828_171237.txt";//350 nodes, 600 edges

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
    G.computeAdjacencyMatrix();
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;
    std::cout << "Execution time graph reading: " << duration.count() << " seconds" << std::endl << std::endl;
    // std::cout << "Num of nodes: " << G.num_of_nodes() << " and Num of Edges: " << G.num_of_edges() << std::endl;
    // G1 = G;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto partitions = algorithmsRunner(G, algorithmName, numPartitions, numThreads);

    // Print the partitioning result
    std::cout << std::endl << "Partitioning result:" << std::endl;
    for (auto v : partitions) {
        if (v.size() > 0) {
            i++;
            std::cout << "Partition " << i << ": ";
            for (auto j : v) {
                std::cout << j << " ";
            }
            std::cout << std::endl << std::endl;
        }
    }

    return 0;
}