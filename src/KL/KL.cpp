#include <iostream>
#include <vector>
#include <random>
#include "../Graph/Graph.h"
#include <climits>
#include <utility>
#include <algorithm>
#include <set>

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
        maxGain = INT_MIN;
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

std::vector<std::pair<int, int>> heavyEdgeMatching(Graph G) {
    std::vector<std::pair<int, int>> M;
    std::vector<bool> visited(G.num_of_nodes(), false);

    std::vector<Edge> edges = G.getEdges();

    // Sort the edges in decreasing order of weight
    std::sort(edges.begin(), edges.end(), [](const Edge& e1, const Edge& e2) {
        return e1.weight > e2.weight;
    });

    for(auto& edge : edges) {
        if(!visited[edge.n1] && !visited[edge.n2]) {
            M.push_back(std::make_pair(edge.n1, edge.n2));
            visited[edge.n1] = true;
            visited[edge.n2] = true;
        }
    }

    return M;
}


Graph coarsening(Graph G) {
    Graph G1;
    std::map<int, Node> nodes = G.getNodes();

    std::vector<std::pair<int, int>> M = heavyEdgeMatching(G);

    std::cout << "Matching: " << std::endl;
    for (const auto& edge : M) {
        std::cout << "Node " << edge.first << " - Node " << edge.second << std::endl;
    }

    auto matAdj = G.getMatAdj();

    for(auto& edge : M) {
        Coarse* coarse = new Coarse;
        coarse->n1 = edge.first;
        coarse->n2 = edge.second;
        coarse->weight1 = nodes.find(edge.first)->second.weight;
        coarse->weight2 = nodes.find(edge.second)->second.weight;
        G1.setNode(G1.returnLastID(), coarse->weight1 + coarse->weight2, coarse);
    }

    
    //manage unmatched nodes optimized for c++ 17
    /*for (const auto& [id, node]: nodes) {
        // Check if the node is not part of any matching pair
        if (std::find_if(M.begin(), M.end(), [id](const std::pair<int, int>& edge) {
                return edge.first == id || edge.second == id;
            }) == M.end()) {
            // The node is unmatched, add it to G1 with nullptr for Coarse pointer
            G1.setNode(G1.returnLastID(), node.weight);
        }
    }*/
    
    // Process the unmatched nodes
    for (const auto& nodePair : nodes) {
    
        int nodeId = nodePair.first;
        Node& node = nodes[nodeId]; 
        // Check if the node is not part of any matching pair
        bool isMatched = false;
        for (const auto& edgePair : M) {
            if (edgePair.first == nodeId || edgePair.second == nodeId) {
                isMatched = true;
                break;
            }
        }
        // If the node is not matched, add it to G1 with nullptr for Coarse pointer
        if (!isMatched) {
            Coarse* coarse = new Coarse;
            coarse->n1 = nodeId;
            coarse->n2 = nodeId;
            coarse->weight1 = nodes[nodeId].weight;
            coarse->weight2 = nodes[nodeId].weight;
            G1.setNode(G1.returnLastID(), node.weight, coarse);
        }
    }

    std::map<int, Node> newNodes = G1.getNodes();
    std::set<std::pair<int, int>> addedEdges;

    //settare edge del nuovo grafo
    for(auto& edge : M) {

        for(int i=0; i<G.num_of_nodes(); i++){

            if(matAdj[edge.first][i][0] && matAdj[edge.second][i][0]) {
                //dobbiamo unire edge del nodo trovato con il nodo dato dai 2 uniti
                //cerco quindi l'id del nodo nuovo in G1 tramite gli id dei nodi che l'hanno formato
                //forse potremmo fare una map per rendere la ricerca costante e non sequenziale

                int id1 = G1.findNodeIdByCoarseIds(edge.first, edge.second);
                int id2 = G1.findNodeIdByCoarseSingleId(i);
                std::cout<<"id1: "<<id1<<" id2: "<<id2<<std::endl;
                if(id1!=-1 && id2!=-1){

                    if (addedEdges.find({id1, id2}) == addedEdges.end() && addedEdges.find({id2, id1}) == addedEdges.end()) {
                        G1.setEdge(id1, id2, matAdj[edge.first][i][1] + matAdj[edge.second][i][1]);
                        std::cout << "id1: " << id1 << " id2: " << id2 << " weight: " << matAdj[edge.first][i][1] + matAdj[edge.second][i][1] << std::endl;
                        addedEdges.insert({id1, id2}); // Add the edge to the set
                    }

                    
                }
                
            }
            else if(matAdj[edge.first][i][0] == 1 && matAdj[edge.second][i][0] == 0 && i!=edge.first && i!=edge.second){

                int id1 = G1.findNodeIdByCoarseSingleId(edge.first);
                int id2 = G1.findNodeIdByCoarseSingleId(i);
                if(id1!=-1 && id2!=-1){

                    if (addedEdges.find({id1, id2}) == addedEdges.end() && addedEdges.find({id2, id1}) == addedEdges.end()) {
                        G1.setEdge(id1, id2, matAdj[edge.first][i][1]);
                        addedEdges.insert({id1, id2}); // Add the edge to the set
                    }
                    
                }

            }
            else if(matAdj[edge.first][i][0] == 0 && matAdj[edge.second][i][0] == 1 && i!=edge.first && i!=edge.second){

                int id1 = G1.findNodeIdByCoarseSingleId(edge.second);
                int id2 = G1.findNodeIdByCoarseSingleId(i);
                if(id1!=-1 && id2!=-1){

                    if (addedEdges.find({id1, id2}) == addedEdges.end() && addedEdges.find({id2, id1}) == addedEdges.end()) {
                        G1.setEdge(id1, id2, matAdj[edge.second][i][1]);
                        addedEdges.insert({id1, id2}); // Add the edge to the set
                    }
                    
                }
            }
        }
    }
    G1.setSizeNodes(G1.getNodes().size());
    G1.setSizeEdges(G1.getEdges().size());
    // Compute the adjacency matrix and matrix degree for G1
    G1.computeAdjacencyMatrix();
    // G1.computeMatrixDegree();

    return G1;
}

//TODO: fix edges (nodes are fine)
Graph uncoarsening(Graph G1) {
    Graph G;

    // Get the map of nodes from G1
    std::map<int, Node> nodes = G1.getNodes();

    // Uncoarsen the matched nodes first
    for (const auto& nodePair : nodes) {
        int nodeId = nodePair.first;
        const Node& node = nodePair.second;
        if (node.coarse != nullptr) {
            G.setNode(node.coarse->n1, node.coarse->weight1);
            G.setNode(node.coarse->n2, node.coarse->weight2);
            G.setEdge(node.coarse->n1, node.coarse->n2, node.weight);
        } else {
            G.setNode(nodeId, node.weight);
        }
    }

    // Uncoarsen the unmatched nodes
    for (const auto& nodePair : nodes) {
        int nodeId = nodePair.first;
        const Node& node = nodePair.second;
        bool isMatched = false;
        for (const auto& edgePair : G1.getEdges()) {
            if (edgePair.n1 == nodeId || edgePair.n2 == nodeId) {
                isMatched = true;
                break;
            }
        }
        if (!isMatched && node.coarse != nullptr) {
            G.setNode(nodeId, node.weight);
        }
    }

    // Uncoarsen the edges
    for (const auto& edgePair : G1.getEdges()) {
        int n1 = edgePair.n1;
        int n2 = edgePair.n2;
        int weight = edgePair.weight;

        // Find the uncoarsened nodes corresponding to n1 and n2
        int uncoarsened_n1 = G1.findNodeIdByCoarseSingleId(n1);
        int uncoarsened_n2 = G1.findNodeIdByCoarseSingleId(n2);

        // Add the edge to the uncoarsened graph G
        if (uncoarsened_n1 != -1 && uncoarsened_n2 != -1) {
            G.setEdge(uncoarsened_n1, uncoarsened_n2, weight);
        }
    }

    G.setSizeNodes(G.getNodes().size());
    G.setSizeEdges(G.getEdges().size());
    // Compute the adjacency matrix and matrix degree for G
    G.computeAdjacencyMatrix();
    // G.computeMatrixDegree();

    return G;
}