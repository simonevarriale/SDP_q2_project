#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <random>
#include <numeric>
#include <chrono>
#include "../Graph/Graph.h"
#include "BucketPQ.h"

extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
extern int calculateTotalWeight(Graph& graph);

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
        //gains[i] *= graph.getNodes().at(i).weight; // Multiply by node weight
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

std::vector<bool> fm(Graph& graph) {
    // Check if the graph has any edges
    if (graph.getEdges().empty()) {
        // If the graph has no edges, return the partition with all nodes in it
        return std::vector<bool>(graph.num_of_nodes(), false);
    }

    std::map<int, Node> nodes = graph.getNodes();

    // Compute the total weight of the nodes and edges
    int totalWeight = 0;
    for (auto& node : nodes) {
        totalWeight += node.second.weight;
    }
    // for (const Edge& edge : graph.getEdges()) {
    //     totalWeight += edge.weight;
    // }

    // Initialize partition A and B based on node weights
    std::vector<bool> partitionA(graph.num_of_nodes(), false);
    std::vector<bool> partitionB(graph.num_of_nodes(), true);
    int currentWeight = 0;
    for (int i = 0; i < graph.num_of_nodes(); i++) {
        if (currentWeight + graph.getNodeWeight(i) <= totalWeight / 2) {
            partitionA[i] = true;
            partitionB[i] = false;
            currentWeight += graph.getNodeWeight(i);
        }
    }



    // Initialize the best cut and the current cut
    int bestCut = totalWeight;
    int currentCut = bestCut;

    // Initialize the random number generator
    std::random_device rd;
    std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());

    // Repeat the algorithm until no improvement is made
    bool improved = true;
    while (improved) {
        improved = false;

        // Shuffle the nodes
        std::vector<int> nodeOrder(graph.num_of_nodes());
        std::iota(nodeOrder.begin(), nodeOrder.end(), 0);
        std::shuffle(nodeOrder.begin(), nodeOrder.end(), gen);


        // Move nodes from partition B to partition A
        for (int node : nodeOrder) {
            if (partitionB[node]) {
                partitionA[node] = true;
                partitionB[node] = false;

                // Update the cut
                for (const Edge& edge : graph.getEdges()) {
                    if (edge.n1 == node) {
                        if (partitionA[edge.n2]) {
                            currentCut += edge.weight;
                        }
                        else {
                            currentCut -= edge.weight;
                        }
                    }
                    else if (edge.n2 == node) {
                        if (partitionA[edge.n1]) {
                            currentCut += edge.weight;
                        }
                        else {
                            currentCut -= edge.weight;
                        }
                    }
                }

                // Check if the cut has improved
                if (currentCut < bestCut) {
                    bestCut = currentCut;
                    improved = true;
                }
                else {
                    partitionA[node] = false;
                    partitionB[node] = true;
                }
            }
        }
    }

    // Return the partition with the smaller cut
    if (bestCut == totalWeight) {
        // If the cut was not improved, return the original partition
        return partitionA;
    }
    else if (currentCut == bestCut) {
        // If the cut was improved, return the partition with the smaller cut
        return partitionA;
    }
    else {
        return partitionB;
    }
}


std::vector<bool> fiducciaMattheyses2(Graph& graph, int maxIterations,std::vector<bool> partitionA = {} ) {
    
    if (partitionA.empty()) {
        partitionA.resize(graph.num_of_nodes(), false); // Initial partition A
        // Initialize partition A and B based on node weights
        int totalWeight = calculateTotalWeight(graph);
        int currentWeight = 0;
        for (int i = 0; i < graph.num_of_nodes(); i++) {
            if (currentWeight + graph.getNodeWeight(i) <= totalWeight / 2) {
                partitionA[i] = true;
                currentWeight += graph.getNodeWeight(i);
            }
        }
    }

    // Initialize the gains vector to the correct size
    std::vector<int> gains(graph.num_of_nodes(), 0);
    std::vector<bool> lock(graph.num_of_nodes(), 0);

    // Compute the initial gains for each node in the graph
    computeInitialGains(graph, partitionA, gains);

    // Track the best partitioning solution found during the process
    std::vector<bool> bestPartitionA = partitionA;
    BucketPriorityQueue bucket1,bucket2, bucket;

    int bestCutSize = calculateCutSize(graph, partitionA);
    std::vector<int> L;

    int numIterations = 0;
    int maxGainIdx = 0;
    int maxGain = 0;
    
    do {
        int maxGain = INT_MIN;
        int maxGainNode = -1;
        bool moveToPartitionA = false;

        for(int i=0; i<graph.num_of_nodes(); i++){
            if(!partitionA[i]){
                bucket1.insert(i, gains[i]);
            }
            else{
                bucket2.insert(i, gains[i]);
            }
        }
        

        int index=0;
        while(std::find(lock.begin(), lock.end(), false) != lock.end()){
            
            if(index%2==0){
                bucket=bucket1;
            }
            else{
                bucket=bucket2;
            }

            int max = bucket.extractMax();
            L.push_back(max);

            lock[max] = true;

            partitionA[max] = !partitionA[max];
            
            
            for (const Edge& edge : graph.getEdges()) {
                int w = (edge.n1 == max) ? edge.n2 : edge.n1;
                
                if (!lock[w]) {
                    int weight = edge.weight;
                    gains[w] += (partitionA[max]) ? weight : -weight;
                    if (partitionA[w]) {
                        bucket1.insert(w, gains[w]);
                        bucket2.insert(w, 0); // To avoid errors when extracting from an empty queue
                    } else {
                        bucket2.insert(w, gains[w]);
                        bucket1.insert(w, 0);
                    }
                }
            }
            
            
            int sum = 0;
            int n = std::min(std::count(partitionA.begin(), partitionA.end(), true), std::count(partitionA.begin(), partitionA.end(), false));
            

            //problem
            for (int k = 0; k < n; ++k) {
                sum += gains[L[k]];
                if (sum > maxGain) {
                    maxGain = sum;
                    maxGainIdx = k;
                }
            }
            

            if (maxGain > 0) {
                // Apply changes to the partition based on the maxGainIdx
                for (int i = 0; i <= maxGainIdx; ++i) {
                    int v = L[i];
                    partitionA[v] = !partitionA[v];
                }
            }
            

        }


        
       
    } while (maxGain <= 0 && numIterations<maxIterations);

    return partitionA;
}