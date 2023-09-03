#include "../utils/utils.h"
#include <iostream>
#include <vector>
#include <random>
#include "../Graph/Graph.h"
#include <climits>
#include <utility>
#include <algorithm>
#include <set>
#include <unordered_set>

#define SEQUENTIAL_THRESHOLD 20

typedef struct {
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

std::vector<bool> kernighanLin(Graph& graph, std::vector<bool> partitionA = {}) {

    /*if (partitionA.empty()) {
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
    }*/
    if (partitionA.empty()) {
        partitionA.resize(graph.num_of_nodes(), false); // Initial partition A
        // Initialize partition A and B based on node weights
        for (int i = 0; i < graph.num_of_nodes() / 2; i++) {
            partitionA[i] = true;
        }
    }


    std::vector<bool> lock(graph.num_of_nodes(), false); // Locks for each node (true if locked, false if unlocked
    std::vector<int> netGains(graph.num_of_nodes(), 0);  // Net gains for each node
    std::vector<std::vector<int>> g(graph.num_of_nodes(), std::vector<int>(graph.num_of_nodes(), 0));
    G_Max gMax;
    std::vector<G_Max> vecGMax;
    int maxGain = 0;

    bool hasImproved = true;
    std::vector<bool> prevPartition = partitionA;
    // std::cout << "Initial partition size: " << prevPartition.size() << ", Graph size: " << graph.num_of_nodes() << ", Cut size: "<< calculateCutSize(graph, partitionA) <<std::endl;


    int cumulativeWeightA = 0;
    int cumulativeWeightB = 0;

    do {
        gMax.gMax = INT_MIN;
        computeNetGains(graph, partitionA, netGains);
        for (int k = 0; k < graph.num_of_nodes() / 2; k++) {
            for (int i = 0; i < graph.num_of_nodes(); i++) {
                for (int j = i + 1; j < graph.num_of_nodes(); j++) {
                    //prova nodo unlocked
                    if (!lock[i] && !lock[j]) {
                        g[i][j] = netGains[i] + netGains[j] - 2 * graph.getMatAdj()[i][j][1];
                        if (g[i][j] > gMax.gMax) {
                            gMax.gMax = g[i][j];
                            gMax.i = i;
                            gMax.j = j;
                        }
                    }
                }
            }

            vecGMax.push_back(gMax);

            lock[gMax.i] = true;
            lock[gMax.j] = true;

            for (int i = 0; i < graph.num_of_nodes(); i++) {
                if (!lock[i]) {
                    if (!partitionA[i]) {
                        netGains[i] = netGains[i] + 2 * graph.getMatAdj()[i][gMax.i][1] - 2 * graph.getMatAdj()[i][gMax.j][1];
                    }
                    else {
                        netGains[i] = netGains[i] - 2 * graph.getMatAdj()[i][gMax.i][1] + 2 * graph.getMatAdj()[i][gMax.j][1];
                    }
                }
            }
        }

        // Find k, such that Gk = Sum(Pi) gains[i] is maximized
        maxGain = INT_MIN;
        int maxGainIdx = 0;
        int sum = 0;
        for (int k = 0; k < vecGMax.size(); ++k) {
            // Calculate the sum of gains[0] to gains[k]

            for (int i = 0; i <= k; ++i) {
                sum += vecGMax[i].gMax;
            }
            if (sum > maxGain) {
                maxGain = sum;
                maxGainIdx = k;
            }
        }
        // Update cumulative weights after the moves
        cumulativeWeightA = 0;
        cumulativeWeightB = 0;
        for (int i = 0; i < graph.num_of_nodes(); ++i) {
            if (!partitionA[i]) {
                cumulativeWeightA += graph.getNodeWeight(i);
            }
            else {
                cumulativeWeightB += graph.getNodeWeight(i);
            }
        }

        if (maxGain > 0) {
            for (int i = 0; i <= maxGainIdx; i++) {

                if (!partitionA[vecGMax[i].i]) {

                    if (
                        (cumulativeWeightA - graph.getNodeWeight(vecGMax[i].i) + graph.getNodeWeight(vecGMax[i].j) <=
                            (cumulativeWeightB + graph.getNodeWeight(vecGMax[i].i) - graph.getNodeWeight(vecGMax[i].j)) * 1.2)
                        &&
                        (cumulativeWeightA - graph.getNodeWeight(vecGMax[i].i) + graph.getNodeWeight(vecGMax[i].j) >=
                            (cumulativeWeightB + graph.getNodeWeight(vecGMax[i].i) - graph.getNodeWeight(vecGMax[i].j)) * 0.8)

                        ) {

                        partitionA[vecGMax[i].i] = !partitionA[vecGMax[i].i];
                        partitionA[vecGMax[i].j] = !partitionA[vecGMax[i].j];
                    }

                }
                else {
                    if (
                        (cumulativeWeightA + graph.getNodeWeight(vecGMax[i].i) - graph.getNodeWeight(vecGMax[i].j) <=
                            (cumulativeWeightB - graph.getNodeWeight(vecGMax[i].i) + graph.getNodeWeight(vecGMax[i].j)) * 1.2)
                        &&
                        (cumulativeWeightA + graph.getNodeWeight(vecGMax[i].i) - graph.getNodeWeight(vecGMax[i].j) >=
                            (cumulativeWeightB - graph.getNodeWeight(vecGMax[i].i) + graph.getNodeWeight(vecGMax[i].j)) * 0.8)

                        ) {
                        partitionA[vecGMax[i].i] = !partitionA[vecGMax[i].i];
                        partitionA[vecGMax[i].j] = !partitionA[vecGMax[i].j];
                    }

                }


            }
        }

        // Clear the locks
        for (int i = 0; i < graph.num_of_nodes(); ++i) {
            lock[i] = false;
        }
        vecGMax.clear();

        // Check for improvement in the partition
        hasImproved = (prevPartition != partitionA);
        prevPartition = partitionA;
    } while (maxGain <= 0 && hasImproved);

    cumulativeWeightA = 0;
    cumulativeWeightB = 0;
    for (int i = 0; i < graph.num_of_nodes(); ++i) {
        if (partitionA[i]) {
            cumulativeWeightA += graph.getNodeWeight(i);
        }
        else {
            cumulativeWeightB += graph.getNodeWeight(i);
        }
    }
    // std::cout << "Cut size final: "<< calculateCutSize(graph, partitionA) <<std::endl;
    // std::cout<<"Weight partition: "<<cumulativeWeightA<<std::endl;
    return partitionA;
}

std::vector<bool> kernighanLin1(Graph& G, std::vector<bool> partition = {}) {
    int numNodes = G.num_of_nodes();
    int halfNumNodes = numNodes / 2;

    //std::vector<bool> partition(numNodes, false); // Partition assignment
    bool hasImproved = true;


    if (partition.empty()) {
        partition.resize(G.num_of_nodes(), false); // Initial partition A
        // Initialize partition A and B based on node weights
        int totalWeight = calculateTotalWeight(G);
        int currentWeight = 0;
        for (int i = 0; i < G.num_of_nodes(); i++) {
            if (currentWeight + G.getNodeWeight(i) <= totalWeight / 2) {
                partition[i] = true;
                currentWeight += G.getNodeWeight(i);
            }
        }
    }
    std::vector<bool> prevPartition = partition;
    // std::cout << "Initial partition KL: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //        std::cout << partition[i] << " ";
    //    }
    // std::cout << std::endl;



    int maxIterations = numNodes / 2; // Maximum iterations as per the pseudocode
    int g_max = 1;
    int i = 0;


    do {
        std::vector<int> D(numNodes, 0); // Initialize D values
        std::vector<int> gv, av, bv;     // Lists to store gains and nodes

        computeInitialGains(G, partition, D);
        for (int n = 0; n < halfNumNodes; ++n) {
            // Compute D values for all nodes in A and B
            

            // Find a and b that maximize g = D[a] + D[b] - 2 * MatAdj[a][b]
            int bestA = -1, bestB = -1;
            int bestG = std::numeric_limits<int>::min();
            for (int a = 0; a < numNodes; ++a) {
                if (partition[a]) {
                    for (int b = 0; b < numNodes; ++b) {
                        if (!partition[b]) {
                            int g = D[a] + D[b] - 2 * G.getMatAdj()[a][b][1];
                            if (g > bestG) {
                                bestG = g;
                                bestA = a;
                                bestB = b;
                            }
                        }
                    }
                }
            }

            // Remove a and b from further consideration in this pass
            partition[bestA] = false;
            partition[bestB] = true;

            // Update gv, av, and bv
            gv.push_back(bestG);
            av.push_back(bestA);
            bv.push_back(bestB);

            // Update D values for the remaining nodes in A and B
            for (int i = 0; i < numNodes; ++i) {
                if (partition[i]) {
                    D[i] -= 2 * G.getMatAdj()[bestA][i][1] - 2 * G.getMatAdj()[bestB][i][1];
                }
                else {
                    D[i] += 2 * G.getMatAdj()[bestA][i][1] - 2 * G.getMatAdj()[bestB][i][1];
                }
            }

        }

        // Find k which maximizes g_max, the sum of gv[1], ..., gv[k]
        int maxIdx = -1;
        int maxSum = 0;
        for (int k = 0; k < halfNumNodes; ++k) {
            int currentSum = 0;
            for (int j = 0; j <= k; ++j) {
                currentSum += gv[j];
            }
            if (currentSum > maxSum) {
                maxSum = currentSum;
                maxIdx = k;
            }
        }

        // If g_max > 0, exchange nodes between partitions av[1], av[2], ..., av[k] and bv[1], bv[2], ..., bv[k]
        if (maxIdx >= 0) {
            for (int j = 0; j <= maxIdx; ++j) {
                int temp = av[j];
                av[j] = bv[j];
                bv[j] = temp;
                partition[av[j]] = true;
                partition[bv[j]] = false;
            }
        }

        // Update g_max
        g_max = maxSum;

        // Check for improvement in the partition
        hasImproved = (prevPartition != partition);
        prevPartition = partition;
        i++;

    } while (g_max > 0 && hasImproved && i < maxIterations);


    // Remaining code for calculating balance and cut size...
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

    std::cout << "Partition Balance Factor KL: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
    std::cout << "Cut size KL: " << calculateCutSize(G, partition) << std::endl;

    return partition;
}

//Recursive imeplementation of multilevel KL
std::vector<bool> multilevel_KL(Graph& G) {

    // Check if the graph is small enough to be partitioned using a sequential algorithm
    // if (G.num_of_nodes() <= SEQUENTIAL_THRESHOLD) {
    //     // G.printGraph();
    //     return kernighanLin(G);
    // }

    // // Coarsen the graph
    // Graph G1 = coarsening(G);

    // // Recursively partition the coarser graph
    // std::vector<bool> partition_coarse = multilevel_KL(G1);

    // Uncoarsen the partitioning
    // std::vector<bool> partition_uncoarse(G.num_of_nodes(), false);
    // std::pair<int, int> ids;

    // for (int i = 0; i < G1.num_of_nodes(); i++) {
    //     ids = G1.getCoarseIdsById(i);
    //     partition_uncoarse[ids.first] = partition_coarse[i];
    //     partition_uncoarse[ids.second] = partition_coarse[i];
    // }

    // // Refine the partitioning using the Kernighan-Lin algorithm
    // return kernighanLin(G);

    std::vector<std::unordered_map<int, std::pair<int, int>>> coarsenG;

    std::vector<bool> newPartition;
    //coarsenG.push_back(G);
    std::cout << "Coarsening" << std::endl;
    coarsenG.push_back(coarsenGraph(G));
    std::cout << "Fine Coarsening con " << G.num_of_nodes() << " nodes and " << G.num_of_edges() << " edges" << std::endl;
    int i = 0;
    while (coarsenG[i].size() > 50) {
        std::cout << "Coarsening" << std::endl;
        coarsenG.push_back(coarsenGraph(G));
        std::cout << "Fine Coarsening con " << G.num_of_nodes() << " nodes and " << G.num_of_edges() << " edges" << std::endl;
        i++;
    }

    //std::cout<<"Num nodes: "<<coarsenG.at(coarsenG.size() - 1).num_of_nodes()<<std::endl;
    std::vector<bool> partition = kernighanLin1(G);

    for (int i = coarsenG.size() - 1; i >= 0; i--) {

        newPartition = uncoarsening2(coarsenG[i], partition);


        // if(newPartition.size() < G.num_of_nodes()/2){
        //     partition = kernighanLin1(coarsenG.at(i-1), newPartition);
        // }
        // else{
        //     partition = newPartition;
        // }
        partition = newPartition;

    }


    return partition;

}

std::vector<bool> old_multilevel_KL(Graph& G) {
    std::vector<Graph> coarsenG;

    std::vector<bool> newPartition;
    coarsenG.push_back(G);
    int i = 0;
    while (coarsenG[i].num_of_nodes() > SEQUENTIAL_THRESHOLD) {
        coarsenG.push_back(coarsening(coarsenG.at(i)));
        // coarsenG[i].printGraph();
        i++;
    }
    std::vector<bool> partition = kernighanLin(coarsenG.at(coarsenG.size() - 1));

    for (int i = coarsenG.size() - 2; i >= 0; i--) {
        newPartition = uncoarsening(coarsenG.at(i + 1), partition, coarsenG.at(i).num_of_nodes());
        partition = kernighanLin(coarsenG.at(i), newPartition);
    }

    return partition;
}


// Original Kernighan-Lin algorithm
// std::vector<bool> kernighanLinOG(Graph& graph) {
//     int sizeNodes = graph.num_of_nodes();
//     int halfSize = sizeNodes / 2;

//     auto partitionA = generateRandomBisection(graph); // Initial partition A
//     std::cout << "Initial partition: " << std::endl;
//     for (int i = 0; i < sizeNodes; ++i) {
//         std::cout << partitionA[i] << " ";
//     }
//     std::cout << std::endl;

//     std::vector<int> netGains(sizeNodes, 0); // Net gains for each node
//     // Perform iterations until no further improvement can be made
//     bool improved = true;
//     while (improved) {
//         improved = false;
//         computeNetGains(graph, partitionA, netGains);

//         for (int i = 0; i < halfSize; ++i) {
//             int maxGainIdx = -1;
//             int maxGain = 0;

//             // Find the node in set A with the maximum net gain
//             for (int j = 0; j < sizeNodes; ++j) {
//                 if (!partitionA[j] && netGains[j] > maxGain) {
//                     maxGainIdx = j;
//                     maxGain = netGains[j];
//                 }
//             }

//             int nodeToMove = maxGainIdx;
//             partitionA[nodeToMove] = true;

//             // Find the corresponding node in set B with the maximum net gain
//             for (int j = 0; j < sizeNodes; ++j) {
//                 if (partitionA[j] && netGains[j] > maxGain) {
//                     maxGainIdx = j;
//                     maxGain = netGains[j];
//                 }
//             }

//             // Check if there is an improvement in cut size after the swap
//             int cutBefore = calculateCutSize(graph, partitionA);
//             partitionA[maxGainIdx] = false;
//             int cutAfter = calculateCutSize(graph, partitionA);
//             if (cutAfter < cutBefore) {
//                 improved = true;
//             } else {
//                 partitionA[nodeToMove] = false;
//             }
//         }
//     }
//     return partitionA;
// }

// Kernighan-Lin algorithm
// std::vector<bool> kernighanLin(Graph& graph) {
//     // auto partitionA = generateRandomBisection(graph); // Initial partition A
//     std::vector<bool> partitionA(graph.num_of_nodes(), false); // Initial partition A
//     for (int i = 0; i < graph.num_of_nodes() / 2; ++i) {
//         partitionA[i] = true;
//     }
//     std::vector<bool> lock(graph.num_of_nodes(), false); // Locks for each node (true if locked, false if unlocked
//     std::vector<int> netGains(graph.num_of_nodes(), 0);  // Net gains for each node
//     std::vector<std::vector<int>> g(graph.num_of_nodes(), std::vector<int>(graph.num_of_nodes(), 0));
//     G_Max gMax;

//     std::vector<G_Max> vecGMax;
//     int maxGain = 0;

//     bool hasImproved = true;
//     std::vector<bool> prevPartition = partitionA;

//     do {
//         gMax.gMax = INT_MIN;
//         computeNetGains(graph, partitionA, netGains);


//         for (int k = 0; k < graph.num_of_nodes() / 2; k++) {

//             for (int i = 0; i < graph.num_of_nodes(); i++) {
//                 for (int j = i + 1; j < graph.num_of_nodes(); j++) {
//                     //prova nodo unlocked
//                     if (!lock[i] && !lock[j]) {
//                         g[i][j] = netGains[i] + netGains[j] - 2 * graph.getMatAdj()[i][j][1];
//                         if (g[i][j] > gMax.gMax) {
//                             gMax.gMax = g[i][j];
//                             gMax.i = i;
//                             gMax.j = j;
//                         }
//                     }
//                 }
//             }

//             vecGMax.push_back(gMax);

//             lock[gMax.i] = true;
//             lock[gMax.j] = true;

//             for (int i = 0; i < graph.num_of_nodes(); i++) {
//                 if (!lock[i]) {
//                     if (!partitionA[i]) {
//                         netGains[i] = netGains[i] + 2 * graph.getMatAdj()[i][gMax.i][1] - 2 * graph.getMatAdj()[i][gMax.j][1];
//                     }
//                     else {
//                         netGains[i] = netGains[i] - 2 * graph.getMatAdj()[i][gMax.i][1] + 2 * graph.getMatAdj()[i][gMax.j][1];
//                     }
//                 }
//             }
//         }

//         // Find k, such that Gk = Sum(Pi) gains[i] is maximized
//         maxGain = INT_MIN;
//         int maxGainIdx = 0;
//         int sum = 0;
//         for (int k = 0; k < vecGMax.size(); ++k) {
//             // Calculate the sum of gains[0] to gains[k]

//             for (int i = 0; i <= k; ++i) {
//                 sum += vecGMax[i].gMax;
//             }
//             if (sum > maxGain) {
//                 maxGain = sum;
//                 maxGainIdx = k;
//             }
//         }

//         if (maxGain > 0) {
//             for (int i = 0; i <= maxGainIdx; i++) {
//                 partitionA[vecGMax[i].i] = !partitionA[vecGMax[i].i];
//                 partitionA[vecGMax[i].j] = !partitionA[vecGMax[i].j];
//             }
//         }

//         // Clear the locks
//         for (int i = 0; i < graph.num_of_nodes(); ++i) {
//             lock[i] = false;
//         }
//         vecGMax.clear();

//         // Check for improvement in the partition
//         hasImproved = (prevPartition != partitionA);
//         prevPartition = partitionA;

//     } while (maxGain <= 0 && hasImproved);

//     return partitionA;
// }

// TODO: fix edges (nodes are fine)
// Graph uncoarsening(Graph G1) {
//     Graph G;

//     // Get the map of nodes from G1
//     std::map<int, Node> nodes = G1.getNodes();
//     std::set<std::pair<int, int>> addedEdges;

//     // Uncoarsen the matched nodes first
//     for (const auto& nodePair : nodes) {
//         int nodeId = nodePair.first;
//         const Node& node = nodePair.second;
//         if (node.coarse != nullptr) {
//             G.setNode(node.coarse->n1, node.coarse->weight1);
//             G.setNode(node.coarse->n2, node.coarse->weight2);
//             // G.setEdge(node.coarse->n1, node.coarse->n2, node.coarse->adj[0][1][1]);
//         }
//         else {
//             G.setNode(nodeId, node.weight);
//         }
//     }

//     // Uncoarsen the unmatched nodes
//     for (const auto& nodePair : nodes) {
//         int nodeId = nodePair.first;
//         const Node& node = nodePair.second;
//         bool isMatched = false;
//         for (const auto& edgePair : G1.getEdges()) {
//             if (edgePair.n1 == nodeId || edgePair.n2 == nodeId) {
//                 isMatched = true;
//                 break;
//             }
//         }
//         if (!isMatched && node.coarse != nullptr) {
//             G.setNode(nodeId, node.weight);
//         }
//     }

//     G.setSizeNodes(G.getNodes().size());
//     // Uncoarsen the edges
//     for (const auto& node : nodes) {
//         Coarse* coarse = node.second.coarse;

//         for (int i = 0; i < G.num_of_nodes(); i++) {

//             /*
//             if (coarse->adj[0][i][0] == 1) {
//                 int n1 = coarse->n1;
//                 int n2 = coarse->n2;
//                 int weight = coarse->adj[0][i][1];
//                 if (n1 != n2) { // se il nodo 1 e 2 sono diversi
//                     if (addedEdges.find({ n1, n2 }) == addedEdges.end() && addedEdges.find({ n2, n1 }) == addedEdges.end()) {
//                         G.setEdge(n1, n2, weight);
//                         addedEdges.insert({ n1, n2 });
//                     }
//                 }
//             }*/


//             if (coarse->adj[0][i][0] == 1) {
//                 int n1 = coarse->n1;
//                 int n2 = i;
//                 int weight = coarse->adj[0][i][1];
//                 if (n1 != n2) { // se il nodo 1 e 2 sono diversi
//                     if (addedEdges.find({ n1, n2 }) == addedEdges.end() && addedEdges.find({ n2, n1 }) == addedEdges.end()) {
//                         G.setEdge(n1, n2, weight);
//                         addedEdges.insert({ n1, n2 });
//                     }
//                 }
//             }

//             if (coarse->adj.size() > 1) {
//                 if (coarse->adj[1][i][0] == 1) {
//                     int n2 = coarse->n2;
//                     int n1 = i;
//                     int weight = coarse->adj[1][i][1];
//                     if (n1 != n2) { // se il nodo 1 e 2 sono diversi
//                         if (addedEdges.find({ n1, n2 }) == addedEdges.end() && addedEdges.find({ n2, n1 }) == addedEdges.end()) {
//                             G.setEdge(n1, n2, weight);
//                             addedEdges.insert({ n1, n2 });
//                         }
//                     }
//                 }
//             }


//         }
//         // std::cout << std::endl;
//     }

//     G.setSizeEdges(G.getEdges().size());
//     // Compute the adjacency matrix and matrix degree for G
//     G.computeAdjacencyMatrix();
//     // G.computeMatrixDegree();

//     return G;
// }


// Graph uncoarsening(Graph G1) {
//     Graph G;

//     // Get the map of nodes from G1
//     std::map<int, Node> nodes = G1.getNodes();

//     // A set to keep track of the added edges
//     std::set<std::pair<int, int>> addedEdges;

//     // Uncoarsen the matched nodes first
//     for (const auto& nodePair : nodes) {
//         int nodeId = nodePair.first;
//         const Node& node = nodePair.second;
//         if (node.coarse != nullptr) {
//             G.setNode(node.coarse->n1, node.coarse->weight1);
//             G.setNode(node.coarse->n2, node.coarse->weight2);

//             // Check if the edge between the two coarse nodes is already added
//             if (!addedEdges.count({node.coarse->n1, node.coarse->n2})) {
//                 G.setEdge(node.coarse->n1, node.coarse->n2, node.coarse->adj[0][1][1]);

//                 // Add the edge to the set of added edges
//                 addedEdges.insert({node.coarse->n1, node.coarse->n2});
//             }
//         } else {
//             G.setNode(nodeId, node.weight);
//         }
//     }

//     // Uncoarsen the unmatched nodes
//     for (const auto& nodePair : nodes) {
//         int nodeId = nodePair.first;
//         const Node& node = nodePair.second;
//         bool isMatched = false;
//         for (const auto& edgePair : G1.getEdges()) {
//             if (edgePair.n1 == nodeId || edgePair.n2 == nodeId) {
//                 isMatched = true;
//                 break;
//             }
//         }
//         if (!isMatched && node.coarse != nullptr) {
//             G.setNode(nodeId, node.weight);
//         }
//     }

//     // Uncoarsen the edges
//     std::vector<Edge> edges = G1.getEdges();
//     for (const auto& edge : edges) {
//         int n1 = edge.n1;
//         int n2 = edge.n2;
//         int weight = edge.weight;
//         Coarse* coarse1 = nodes[n1].coarse;
//         Coarse* coarse2 = nodes[n2].coarse;

//         std::cout << "n1: " << n1 << " n2: " << n2 << std::endl;
//         std::cout << "coarse1: " << coarse1->n1 << "-" << coarse1->n2 << std::endl;
//         std::cout << "coarse2: " << coarse2->n1 << "-" << coarse2->n2 << std::endl;
//         for()
//         // If both nodes are part of a coarse node, add an edge between the two coarse nodes
//         if (coarse1 != nullptr && coarse2 != nullptr) {
//             int coarseN1 = coarse1->n1;
//             int coarseN2 = coarse2->n1;

//             // Check if the edge between the two coarse nodes is already added
//             if (!addedEdges.count({coarseN1, coarseN2})) {
//                 G.setEdge(coarseN1, coarseN2, coarse1->adj[0][1][1]);

//                 // Add the edge to the set of added edges
//                 addedEdges.insert({coarseN1, coarseN2});
//             }
//         }
//         // If only one node is part of a coarse node, add an edge between the node and the coarse node
//         else if (coarse1 != nullptr && coarse2 == nullptr) {
//             int coarseN1 = coarse1->n1;

//             // Check if the edge between the node and the coarse node is already added
//             if (!addedEdges.count({n2, coarseN1})) {
//                 G.setEdge(n2, coarseN1, coarse1->adj[0][1][1]);

//                 // Add the edge to the set of added edges
//                 addedEdges.insert({n2, coarseN1});
//             }
//         } else if (coarse1 == nullptr && coarse2 != nullptr) {
//             int coarseN2 = coarse2->n1;

//             // Check if the edge between the node and the coarse node is already added
//             if (!addedEdges.count({n1, coarseN2})) {
//                 G.setEdge(n1, coarseN2, coarse2->adj[0][1][1]);

//                 // Add the edge to the set of added edges
//                 addedEdges.insert({n1, coarseN2});
//             }
//         }
//     }

//     G.setSizeNodes(G.getNodes().size());
//     G.setSizeEdges(G.getEdges().size());

//     // Compute the adjacency matrix for G
//     G.computeAdjacencyMatrix();

//     return G;
// }

// Implementation of multilevel KL
// std::vector<bool> multilevel_KL(Graph& G) {

//     std::vector<bool> partition;
//     std::vector<bool> partition_uncoarse(G.num_of_nodes(), false);

//     Graph G1 = coarsening(G);
//     auto nodes = G1.getNodes();

//     partition = kernighanLin(G1);

//     std::pair<int, int> ids;

//     for (int i = 0; i < G1.num_of_nodes(); i++) {
//         ids = G1.getCoarseIdsById(i);
//         partition_uncoarse[ids.first] = partition[i];
//         partition_uncoarse[ids.second] = partition[i];
//     }

//     return partition_uncoarse;
// }