#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include "../Graph/Graph.h"

extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);

// Function to compute the initial gain for each node in the graph
void computeInitialGains(Graph& graph, const std::vector<bool>& partitionA, std::vector<int>& gains) {
    gains.resize(graph.num_of_nodes(), 0);
    for (int i = 0; i < graph.num_of_nodes(); i++) {
        for (int j = 0; j < graph.num_of_nodes(); j++) {
            if (partitionA[i] != partitionA[j]) {
                gains[i] += graph.getMatAdj()[i][j][1];
            }
            else {
                gains[i] -= graph.getMatAdj()[i][j][1];
            }
        }
        gains[i] *= graph.getNodes().at(i).weight; // Multiply by node weight
    }
}

// Function to calculate the gain for a node when moved to the other partition
int calculateNodeGain(Graph& graph, const std::vector<bool>& partitionA, std::vector<int>& gains, int node, bool moveToPartitionA) {
    int gain = 0;

    for (int i = 0; i < graph.num_of_nodes(); i++) {
        if (partitionA[i] != moveToPartitionA) {
            gain += graph.getMatAdj()[node][i][1];
        }
        else {
            gain -= graph.getMatAdj()[node][i][1];
        }
    }
    gain *= graph.getNodes().at(node).weight; // Multiply by node weight
    return gain;
}

// Function to calculate the gain for a node when moved to the other partition
//gain = (sum of weights of edges to nodes in same partition) - (sum of weights of edges to nodes in other partition) + (sum of weights of nodes in same partition) - (sum of weights of nodes in other partition)
// int calculateNodeGain(Graph& graph, std::vector<bool>& partitionA, std::vector<int>& gains, int node, bool moveToPartitionA) {
//     int gain = 0;

//     // Calculate the gain due to edges
//     for (int i = 0; i < graph.num_of_nodes(); i++) {
//         if (i == node) {
//             continue;
//         }

//         int weight = graph.getMatAdj()[node][i][1];
//         if (partitionA[i] == moveToPartitionA) {
//             gain += weight;
//         }
//         else {
//             gain -= weight;
//         }
//     }

//     // Calculate the gain due to nodes
//     if (moveToPartitionA) {
//         gain += graph.getNodes().at(node).weight;
//     }
//     else {
//         gain -= graph.getNodes().at(node).weight;
//     }

//     return gain;
// }

std::vector<bool> fiducciaMattheyses(Graph& graph, int maxIterations) {
    std::vector<bool> partitionA(graph.num_of_nodes(), false); // Initial partition A
    for (int i = 0; i < graph.num_of_nodes() / 2; ++i) {
        partitionA[i] = true;
    }

    // Initialize the gains vector to the correct size
    std::vector<int> gains(graph.num_of_nodes(), 0);

    // Compute the initial gains for each node in the graph
    computeInitialGains(graph, partitionA, gains);

    // Track the best partitioning solution found during the process
    std::vector<bool> bestPartitionA = partitionA;
    int bestCutSize = calculateCutSize(graph, partitionA);

    int numIterations = 0;
    while (numIterations < maxIterations) {
        int maxGain = INT_MIN;
        int maxGainNode = -1;
        bool moveToPartitionA = false;

        for (int i = 0; i < graph.num_of_nodes(); i++) {
            // Calculate gain when moving the node to the other partition
            int gainToA = calculateNodeGain(graph, partitionA, gains, i, true);
            int gainToB = calculateNodeGain(graph, partitionA, gains, i, false);
            int gain = std::abs(gainToA - gainToB);

            if (gain > maxGain) {
                maxGain = gain;
                maxGainNode = i;
                moveToPartitionA = (gainToA >= gainToB);
            }
        }

        // If no improvement is found, stop the algorithm
        if (maxGain <= 0) {
            break;
        }

        // Move the selected node to the other partition and update gains
        partitionA[maxGainNode] = moveToPartitionA;
        for (int i = 0; i < graph.num_of_nodes(); i++) {
            if (partitionA[i] != partitionA[maxGainNode]) {
                gains[i] -= 2 * graph.getMatAdj()[maxGainNode][i][1];
            }
            else {
                gains[i] += 2 * graph.getMatAdj()[maxGainNode][i][1];
            }
        }

        // Update the balance of the partitions

        // Calculate the cut size for the current partitioning
        int currentCutSize = calculateCutSize(graph, partitionA);

        // Track the best partitioning solution found so far
        if (currentCutSize < bestCutSize) {
            bestCutSize = currentCutSize;
            bestPartitionA = partitionA;
        }

        numIterations++;
    }

    return bestPartitionA;
}