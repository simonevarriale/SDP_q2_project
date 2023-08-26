#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <climits>
#include <chrono>
#include <thread>
#include <mutex>
#include <future>
#include <barrier>
#include "./Graph/Graph.h"

extern std::vector<bool> RSB(Graph& G, int p);
extern std::vector<bool> kernighanLin(Graph& graph, std::vector<bool> partitionA = {});
extern Graph coarsening(Graph& graph);
extern std::vector<bool> uncoarsening(Graph G1, std::vector<bool> partition, int graphSize);
extern std::vector<bool> multilevel_KL(Graph& graph);
//extern std::vector<bool> fiducciaMattheyses(Graph& graph, int maxIterations);
extern std::vector<bool> fiducciaMattheyses2(Graph& graph, int maxIterations, std::vector<bool> partitionA = {});
extern std::vector<bool> fm(Graph& graph);
extern double calculateBalanceF(Graph& graph, const std::vector<bool>& partitionA);
extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
extern std::vector<bool> MLRSB(Graph& G, int p);
extern std::vector<bool> Parallel_MLRSB(Graph& G, int p);
extern std::vector<bool> Parallel_RSB(Graph& G, int p);

void read_input(const std::string& filename, Graph* G) {

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return;
    }

    int numNodes;
    int numEdges;

    int n;
    int weight;

    // Read the number of nodes and edges
    inputFile >> numNodes >> numEdges;

    G->setSizeNodes(numNodes);
    G->setSizeEdges(numEdges);

    // Read the nodes
    for (int i = 0; i < numNodes; ++i) {
        inputFile >> n >> weight;
        G->setNode(n, weight);
    }

    int source;
    int destination;
    // Read the edges
    for (int i = 0; i < numEdges; ++i) {
        inputFile >> source >> destination >> weight;
        G->setEdge(source, destination, weight);
        G->incrementDegree(source);
        G->incrementDegree(destination);
    }

    G->computeAdjacencyMatrix();

    inputFile.close();
}

void read_input2(const std::string& filename, Graph* G) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return;
    }

    int numNodes, numEdges;
    inputFile >> numNodes >> numEdges;

    G->setSizeNodes(numNodes);
    G->setSizeEdges(numEdges);

    for (int i = 0; i < numNodes; i++) {
        G->setNode(i, 1);
        std::string line;
        std::getline(inputFile, line);
        std::istringstream iss(line);
        int adjacentNodeId;
        while (iss >> adjacentNodeId) {
            G->setEdge(i, adjacentNodeId - 1, 1);
        }
    }

    inputFile.close();
}

void read_input_parallel(const std::string& filename, Graph* G) {

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return;
    }

    int numNodes, numEdges;
    inputFile >> numNodes >> numEdges;

    G->setSizeNodes(numNodes);
    G->setSizeEdges(numEdges);

    std::vector<std::thread> threads;
    const int numThreads = /*std::thread::hardware_concurrency()*/ 2;
    const int chunkSize = numNodes / numThreads;

    std::mutex mtx;

    for (int t = 0; t < numThreads; t++) {
        const int start = t * chunkSize;
        const int end = (t == numThreads - 1) ? numNodes : (t + 1) * chunkSize - 1;

        threads.emplace_back([&, start, end]() {
            // std::ifstream threadInputFile(filename);
            // if (!threadInputFile.is_open()) {
            //     std::cout << "Failed to open the file: " << filename << std::endl;
            //     return;
            // }
            for (int i = start; i <= end; i++) {
                mtx.lock();
                G->setNode(i, 1);
                mtx.unlock();

                std::string line;
                std::getline(inputFile, line);
                std::istringstream iss(line);

                int adjacentNodeId;

                mtx.lock();
                while (iss >> adjacentNodeId) {
                    G->setEdge(i, adjacentNodeId - 1, 1);
                }
                mtx.unlock();
            }
            //threadInputFile.close();
            });
    }


    for (auto& thread : threads) {
        thread.join();

    }

    //G->computeAdjacencyMatrix();
    // G->setAdjacencyMatrix();
    inputFile.close();
}

// void read_input_parallel(const std::string& filename, Graph* G) {
//     std::ifstream inputFile(filename);
//     if (!inputFile.is_open()) {
//         std::cout << "Failed to open the file: " << filename << std::endl;
//         return;
//     }

//     int numNodes, numEdges;
//     inputFile >> numNodes >> numEdges;

//     G->setSizeNodes(numNodes);
//     G->setSizeEdges(numEdges);

//     std::vector<std::thread> threads;
//     std::vector<std::mutex> mutexes(numNodes);

//     for (int t = 0; t < 2; t++) {
//         threads.emplace_back([&, t]() {
//             for (int i = t; i < numNodes; i += 2) {
//                 mutexes[i].lock();
//                 G->setNode(i, 1);
//                 std::string line;
//                 std::getline(inputFile, line);
//                 std::istringstream iss(line);
//                 int adjacentNodeId;
//                 while (iss >> adjacentNodeId) {
//                     G->setEdge(i, adjacentNodeId - 1, 1);
//                     // G->setEdge(adjacentNodeId - 1, i, 1); // add this line to fix the issue
//                 }
//                 mutexes[i].unlock();
//             }
//             });
//     }

//     for (auto& thread : threads) {
//         thread.join();
//     }

//     inputFile.close();
// }

std::vector<bool> readPartition(const std::string& filename, int numNodes) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return {};
    }

    std::vector<int> values(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        inputFile >> values[i];
    }

    // Define a threshold for converting integers to boolean values
    int threshold = 0;

    // Initialize the partition as a vector of bool
    std::vector<bool> partition(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        partition[i] = (values[i] > threshold);
    }

    inputFile.close();

    return partition;
}

std::vector<bool> readPartitionS(const std::string& filename, int numNodes) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return {};
    }

    std::vector<int> values(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        std::string line;
        std::getline(inputFile, line);
        std::istringstream iss(line);
        iss >> values[i];
    }

    // Define a threshold for converting integers to boolean values
    int threshold = 0;

    // Initialize the partition as a vector of bool
    std::vector<bool> partition(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        partition[i] = (values[i] > threshold);
    }

    inputFile.close();

    return partition;
}

int main() {

    auto startTime = std::chrono::high_resolution_clock::now();
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;

    Graph G;
    Graph G1;
    std::string file1 = "./data/simple_graph.txt";
    std::string file2 = "./data/test_graph.txt";
    std::string file3 = "./data/connected_graph.txt";
    std::string file4 = "./data/big_graph.txt";
    std::string file5 = "./data/3elt.graph";
    std::string file6 = "./data/simple.graph";
    std::string file7 = "./data/add20.graph";

    //////////////////////////////////////////////////PARTITIONS/////////////////////////////////////////////////
    // auto ourRSBpartition = readPartition("./data/add20_RSB_res.txt", 2395);
    // std::cout << "Partition Read" << std::endl;

    // auto metisRSBpartition = readPartitionS("./data/add20_metis_res.txt", 2395);
    // std::cout << "Partition Read" << std::endl;
    // for (auto i : metisRSBpartition) {
    //     std::cout << metisRSBpartition[i] << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Partition Balance Factor: " << calculateBalanceF(G, metisRSBpartition) << std::endl;
    // std::cout << "Cut size RSB: " << calculateCutSize(G, metisRSBpartition) << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////SEQUENTIAL GRAPH READ//////////////////////////////////////////////
    startTime = std::chrono::high_resolution_clock::now();
    // read_input(file4, &G);
    read_input2(file7, &G);
    // read_input3(file7, &G)
    G.computeAdjacencyMatrix();
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;
    std::cout << "Execution time sequential graph reading: " << duration.count() << " seconds" << std::endl;
    std::cout << "Num of nodes: " << G.num_of_nodes() << " and Num of Edges: " << G.num_of_edges() << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////PARALLEL GRAPH READ//////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // // read_input(file4, &G);
    // read_input_parallel(file7, &G1);
    // // read_input3(file7, &G);
    // // std::cout << "///////////" << "Graphs G" << "///////////" << std::endl;
    // // G.printGraph();
    // // std::cout << "//////////////////////" << std::endl;
    // // std::cout << "///////////" << "Graphs G1" << "///////////" << std::endl;
    // // G1.printGraph();
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // std::cout << "Execution time parallel graph reading: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Num of nodes: " << G1.num_of_nodes() << " and Num of Edges: " << G1.num_of_edges() << std::endl;
    // G1.computeAdjacencyMatrix();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////RSB////////////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // auto partition = RSB(G, 2);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // std::cout << "Execution time RSB: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Partitioning result:" << std::endl;
    // for (auto v : partition) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////PARALLEL RSB////////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // auto partitionParallelRSB = Parallel_RSB(G, 2);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // std::cout << "Execution time RSB: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Partitioning result:" << std::endl;
    // for (auto v : partitionParallelRSB) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////MLRSB////////////////////////////////////////////////////
    startTime = std::chrono::high_resolution_clock::now();
    auto MLRSBpartition = MLRSB(G, 2);
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;
    std::cout << "Execution time MLRSB: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Partitioning result:" << std::endl;
    // for (auto v : MLRSBpartition) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////PARALLEL MLRSB/////////////////////////////////////////////
    startTime = std::chrono::high_resolution_clock::now();
    auto parallelMRSBpartition = Parallel_MLRSB(G, 2);
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;
    std::cout << "Execution time Parallel MLRSB: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Partitioning result:" << std::endl;
    // for (auto v : parallelMRSBpartition) {
    //     std::cout << v << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////KL/////////////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // auto klPart = kernighanLin(G);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // // Print the execution time
    // std::cout << "Execution time KL: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Final partition KL: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //       std::cout << klPart[i] << " ";
    //   }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////MULTILEVEL KL/////////////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // auto multilevel = multilevel_KL(G);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // std::cout << "Execution time multilevel KL: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Final multilevel partition KL: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << multilevel[i] << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////FM/////////////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // auto secondFm = fiducciaMattheyses(G, 2);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // std::cout << "Execution time FM OG: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Final partition FM OG: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << secondFm[i] << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////FM2/////////////////////////////////////////////////////
    // startTime = std::chrono::high_resolution_clock::now();
    // auto fmPart = fiducciaMattheyses2(G, 10, partitionA);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = endTime - startTime;
    // std::cout << "Execution time FM: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Final partition FM: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << fmPart[i] << " ";
    // }
    // std::cout << std::endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return 0;
}