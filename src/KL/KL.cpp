#include <iostream>
#include <vector>
#include <random>
#include "../Graph/Graph.h"

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

// Kernighan-Lin algorithm
void kernighanLin(Graph& graph) {
    int sizeNodes = graph.num_of_nodes();
    int halfSize = sizeNodes / 2;

    std::vector<bool> partitionA(sizeNodes, false); // Initial partition A
    for (int i = 0; i < halfSize; ++i) {
        partitionA[i] = true;
    }

    std::vector<int> netGains(sizeNodes, 0); // Net gains for each node

    // Perform iterations until no further improvement can be made
    bool improved = true;
    while (improved) {
        improved = false;
        computeNetGains(graph, partitionA, netGains);

        for (int i = 0; i < halfSize; ++i) {
            int maxGainIdx = -1;
            int maxGain = 0;

            // Find the node in set A with the maximum net gain
            for (int j = 0; j < sizeNodes; ++j) {
                if (!partitionA[j] && netGains[j] > maxGain) {
                    maxGainIdx = j;
                    maxGain = netGains[j];
                }
            }

            int nodeToMove = maxGainIdx;
            partitionA[nodeToMove] = true;

            // Find the corresponding node in set B with the maximum net gain
            for (int j = 0; j < sizeNodes; ++j) {
                if (partitionA[j] && netGains[j] > maxGain) {
                    maxGainIdx = j;
                    maxGain = netGains[j];
                }
            }

            // Check if there is an improvement in cut size after the swap
            int cutBefore = calculateCutSize(graph, partitionA);
            partitionA[maxGainIdx] = false;
            int cutAfter = calculateCutSize(graph, partitionA);

            if (cutAfter < cutBefore) {
                improved = true;
            } else {
                partitionA[nodeToMove] = false;
            }
        }
    }
}