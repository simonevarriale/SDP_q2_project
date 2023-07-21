#include <iostream>
#include <vector>
#include <random>
#include "../Graph/Graph.h"
#include <climits>

typedef struct  {
    int i;
    int j;
    int gMax;
} G_Max;


// Function to generate a random initial bisection
std::vector<bool> generateRandomBisection(Graph& graph) {
    int sizeNodes = graph.num_of_nodes();
    std::vector<bool> partitionA(sizeNodes, false);

    // Use a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    for (int i = 0; i < sizeNodes; ++i) {
        // Assign each node randomly to either set A (true) or set B (false)
        partitionA[i] = dis(gen);
    }

    return partitionA;
}

// Function to calculate the cut size between two sets A and B
int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA) {
    int cutSize = 0;
    auto matAdj = graph.getMatAdj();

    for (int i = 0; i < graph.num_of_nodes(); ++i) {
        for (int j = i + 1; j < graph.num_of_nodes(); ++j) {
            if (partitionA[i] != partitionA[j]) {
                cutSize += matAdj[i][j][1];
            }
        }
    }
    return cutSize;
}

// Function to compute the net gains for each node in set A
void computeNetGains(Graph& graph, const std::vector<bool>& partitionA, std::vector<int>& netGains) {
    auto matAdj = graph.getMatAdj();
    int sizeNodes = graph.num_of_nodes();

    for (int i = 0; i < sizeNodes; ++i) {
        int gain = 0;
        for (int j = 0; j < sizeNodes; ++j) {
            if (partitionA[i] != partitionA[j]) {
                gain += matAdj[i][j][1];
            } else {
                gain -= matAdj[i][j][1];
            }
        }
        netGains[i] = gain;
    }
}

// std::vector<bool> kernighanLin(Graph& graph) {
    // int sizeNodes = graph.num_of_nodes();
    // int halfSize = sizeNodes / 2;

    // auto partitionA = generateRandomBisection(graph); // Initial partition A
    // std::cout << "Initial partition: " << std::endl;
    // for (int i = 0; i < sizeNodes; ++i) {
    //     std::cout << partitionA[i] << " ";
    // }
    // std::cout << std::endl;

    // std::vector<int> netGains(sizeNodes, 0); // Net gains for each node
    // // Perform iterations until no further improvement can be made
    // bool improved = true;
    // while (improved) {
    //     improved = false;
    //     computeNetGains(graph, partitionA, netGains);

    //     for (int i = 0; i < halfSize; ++i) {
    //         int maxGainIdx = -1;
    //         int maxGain = 0;

    //         // Find the node in set A with the maximum net gain
    //         for (int j = 0; j < sizeNodes; ++j) {
    //             if (!partitionA[j] && netGains[j] > maxGain) {
    //                 maxGainIdx = j;
    //                 maxGain = netGains[j];
    //             }
    //         }

    //         int nodeToMove = maxGainIdx;
    //         partitionA[nodeToMove] = true; 

    //         // Find the corresponding node in set B with the maximum net gain
    //         for (int j = 0; j < sizeNodes; ++j) {
    //             if (partitionA[j] && netGains[j] > maxGain) {
    //                 maxGainIdx = j;
    //                 maxGain = netGains[j];
    //             }
    //         }

    //         // Check if there is an improvement in cut size after the swap
    //         int cutBefore = calculateCutSize(graph, partitionA);
    //         partitionA[maxGainIdx] = false;
    //         int cutAfter = calculateCutSize(graph, partitionA);
    //         if (cutAfter < cutBefore) {
    //             improved = true;
    //         } else {
    //             partitionA[nodeToMove] = false;
    //         }
    //     }
    // }
    // return partitionA;
// }

// Kernighan-Lin algorithm
std::vector<bool> kernighanLin(Graph& graph) {
    // auto partitionA = generateRandomBisection(graph); // Initial partition A
    std::vector<bool> partitionA(graph.num_of_nodes(), false); // Initial partition A
    for (int i = 0; i < graph.num_of_nodes()/2; ++i) {
        partitionA[i] = true;
    }
    std::vector<bool> lock(graph.num_of_nodes(), false); // Locks for each node (true if locked, false if unlocked
    std::vector<int> netGains(graph.num_of_nodes(), 0); // Net gains for each node
    std::vector<std::vector<int>> g(graph.num_of_nodes(), std::vector<int>(graph.num_of_nodes(), 0));
    G_Max gMax;
    gMax.gMax = INT_MIN;
    std::vector<G_Max> vecGMax;
    int maxGain = 0;

    do {
        computeNetGains(graph, partitionA, netGains);
        
        for(int k = 0; k < graph.num_of_nodes()/2; k++) {
            for(int i = 0; i < graph.num_of_nodes(); i++) {
                for (int j = i + 1; j < graph.num_of_nodes(); j++) {
                    g[i][j] = netGains[i] + netGains[j] - 2 * graph.getMatAdj()[i][j][1];
                    if(g[i][j] > gMax.gMax) {
                        gMax.gMax = g[i][j];
                        gMax.i = i;
                        gMax.j = j;
                    }
                }
            }

            vecGMax.push_back(gMax);

            lock[gMax.i] = true;
            lock[gMax.j] = true;

            for(int i = 0; i < graph.num_of_nodes(); i++) {
                if(!lock[i]) {
                    if(!partitionA[i]) {
                        netGains[i] = netGains[i] + 2 * graph.getMatAdj()[i][gMax.i][1] - 2 * graph.getMatAdj()[i][gMax.j][1];
                    } else {
                        netGains[i] = netGains[i] - 2 * graph.getMatAdj()[i][gMax.i][1] + 2 * graph.getMatAdj()[i][gMax.j][1];
                    }
                }
            }
        }

        // Find k, such that Gk = Sum(Pi) gains[i] is maximized
        maxGain = 0;
        int maxGainIdx = 0;
        for (int k = 0; k < vecGMax.size(); ++k) {
            // Calculate the sum of gains[0] to gains[k]
            int sum = 0;
            for (int i = 0; i <= k; ++i) {
                sum += vecGMax[i].gMax;
            }
            if (sum > maxGain) {
                maxGain = sum;
                maxGainIdx = k;
            }
        }

        if(maxGain > 0) {
            for(int i = 0; i <= maxGainIdx; i++) {
                partitionA[vecGMax[i].i] = !partitionA[vecGMax[i].i];
                partitionA[vecGMax[i].j] = !partitionA[vecGMax[i].j];
            }
        }

        // Clear the locks
        for (int i = 0; i < graph.num_of_nodes(); ++i) {
            lock[i] = false;
        }
        vecGMax.clear();

    } while(maxGain <= 0);

    return partitionA;
    
}