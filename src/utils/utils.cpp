#include "utils.h"

double computeMedian(const Eigen::VectorXd& vector) {
    Eigen::VectorXd sortedVector = vector;
    std::sort(sortedVector.data(), sortedVector.data() + sortedVector.size());
    int size = sortedVector.size();

    if (size % 2 == 0) {
        // If the vector size is even, take the average of the middle two values
        int mid = size / 2;
        return 0.5 * (sortedVector(mid - 1) + sortedVector(mid));
    }
    else {
        // If the vector size is odd, return the middle value
        return sortedVector(size / 2);
    }
}

double parallel_computeMedian(const Eigen::VectorXd& vector) {
    Eigen::VectorXd sortedVector = vector;
    std::sort(sortedVector.data(), sortedVector.data() + sortedVector.size());
    int size = sortedVector.size();
    double medianValue;

    if (size % 2 == 0) {
        // If the vector size is even, take the average of the middle two values
        int mid = size / 2;
        std::mutex mutex;
        std::vector<double> middleValues(2);

        // Compute the middle two values in parallel
        std::vector<std::thread> threads;
        for (int i = 0; i < 2; i++) {
            threads.emplace_back([&sortedVector, &middleValues, mid, i, &mutex]() {
                std::lock_guard<std::mutex> lock(mutex);
                middleValues[i] = sortedVector(mid - 1 + i); });
        }

        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }

        medianValue = 0.5 * (middleValues[0] + middleValues[1]);
    }
    else {
        // If the vector size is odd, return the middle value
        medianValue = sortedVector(size / 2);
    }

    return medianValue;
}

// balance factor of 1.0 is perfectly balanced, 0.0 or 2.0 is completely unbalanced
double calculateBalanceFactor(Graph& graph, const std::vector<bool>& partitionA) {
    double weightA = 0.0;
    double weightB = 0.0;
    for (int i = 0; i < graph.num_of_nodes(); i++) {
        if (partitionA[i]) {
            weightA += graph.getNodeWeight(i);
        }
        else {
            weightB += graph.getNodeWeight(i);
        }
    }
    return std::min(weightA, weightB) / std::max(weightA, weightB);
}

int calculatePartitionsWeight(Graph& graph, const std::vector<bool>& partitionA) {
    int weight = 0;
    for (int i = 0; i < graph.num_of_nodes(); i++) {
        if (partitionA[i]) {
            weight += graph.getNodeWeight(i);
        }
    }
    return weight;
}

double calculateBalanceFactorPartitions(Graph& G, const std::vector<std::vector<bool>>& partitions) {
    int partitionSize = partitions.size();
    std::vector<int> partitionsWeight(partitionSize);
    int maxPartitionWeight = 0;
    int minPartitionWeight = std::numeric_limits<int>::max();

    for (int i = 0; i < partitionSize; ++i) {
        partitionsWeight[i] = calculatePartitionsWeight(G, partitions[i]);
        maxPartitionWeight = std::max(maxPartitionWeight, partitionsWeight[i]);
        minPartitionWeight = std::min(minPartitionWeight, partitionsWeight[i]);
    }

    return (double)minPartitionWeight / (double)maxPartitionWeight;
}

int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA) {
    int cutSize = 0;

    for (const Edge& edge : graph.getEdges()) {
        if (partitionA[edge.n1] != partitionA[edge.n2]) {
            cutSize += edge.weight;
            // cutSize++;
        }
    }

    return cutSize;
}

int calculateCutSizeBetweenPartitions(Graph& G, const std::vector<bool>& partitionA, const std::vector<bool>& partitionB) {
    int cutSize = 0;

    const std::vector<std::vector<std::vector<int>>>& matAdj = G.getMatAdj();

    for (int i = 0; i < G.num_of_nodes(); ++i) {
        for (int j = 0; j < G.num_of_nodes(); ++j) {
            if (partitionA[i] && partitionB[j] && matAdj[i][j][0] != -1) {
                cutSize += matAdj[i][j][1]; // Use the weight stored in the adjacency matrix
            }
        }
    }

    return cutSize;
}

int calculateCutSizePartitions(Graph& G, const std::vector<std::vector<bool>>& partitions) {
    int totalCutSize = 0;
    for (int i = 0; i < partitions.size(); ++i) {
        for (int j = i + 1; j < partitions.size(); ++j) {
            totalCutSize += calculateCutSizeBetweenPartitions(G, partitions[i], partitions[j]);
        }
    }
    return totalCutSize;
}

std::vector<int> sortIndices(const Eigen::VectorXd& vec) {
    std::vector<int> indices(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        indices[i] = i;
    }

    // Sort indices based on the values in the vector
    std::sort(indices.begin(), indices.end(), [&vec](int a, int b) { return vec[a] < vec[b]; });

    return indices;
}

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

std::vector<std::pair<int, int>> heavyEdgeMatching(Graph G) {
    std::vector<std::pair<int, int>> M;
    std::vector<bool> visited(G.num_of_nodes(), false);

    std::vector<Edge> edges = G.getEdges();

    // Sort the edges in decreasing order of weight
    std::sort(edges.begin(), edges.end(), [](const Edge& e1, const Edge& e2) { return e1.weight > e2.weight; });

    for (auto& edge : edges) {
        if (!visited[edge.n1] && !visited[edge.n2]) {
            M.push_back(std::make_pair(edge.n1, edge.n2));
            visited[edge.n1] = true;
            visited[edge.n2] = true;
        }
    }

    return M;
}

Graph coarsening(Graph& G) {
    Graph G1;
    std::map<int, Node> nodes = G.getNodes();

    std::vector<std::pair<int, int>> M = heavyEdgeMatching(G);

    auto matAdj = G.getMatAdj();

    for (auto& edge : M) {
        Coarse* coarse = new Coarse;
        coarse->n1 = edge.first;
        coarse->n2 = edge.second;
        coarse->weight1 = nodes.find(edge.first)->second.weight;
        coarse->weight2 = nodes.find(edge.second)->second.weight;
        coarse->adj.push_back(matAdj[edge.first]);
        coarse->adj.push_back(matAdj[edge.second]);
        G1.setNode(G1.returnLastID(), coarse->weight1 + coarse->weight2, coarse);
        delete coarse;
    }

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
            coarse->adj.push_back(matAdj[nodeId]);
            G1.setNode(G1.returnLastID(), node.weight, coarse);
            delete coarse;
        }
    }

    std::map<int, Node> newNodes = G1.getNodes();
    std::set<std::pair<int, int>> addedEdges;

    // settare edge del nuovo grafo
    for (auto& edge : M) {
        for (int i = 0; i < G.num_of_nodes(); i++) {

            if (i == edge.first || i == edge.second) {
                continue;
            }
            if (matAdj[edge.first][i][0] && matAdj[edge.second][i][0]) {
                // dobbiamo unire edge del nodo trovato con il nodo dato dai 2 uniti
                // cerco quindi l'id del nodo nuovo in G1 tramite gli id dei nodi che l'hanno formato
                // forse potremmo fare una map per rendere la ricerca costante e non sequenziale
                int id1 = G1.findNodeIdByCoarseIds(edge.first, edge.second);
                int id2 = G1.findNodeIdByCoarseSingleId(i);
                if (id1 != -1 && id2 != -1) {
                    if (addedEdges.find({ id1, id2 }) == addedEdges.end() && addedEdges.find({ id2, id1 }) == addedEdges.end()) {
                        G1.setEdge(id1, id2, matAdj[edge.first][i][1] + matAdj[edge.second][i][1]);
                        addedEdges.insert({ id1, id2 }); // Add the edge to the set
                    }
                }
            }
            else if (matAdj[edge.first][i][0] == 1 && matAdj[edge.second][i][0] == 0 && i != edge.first && i != edge.second) {
                int id1 = G1.findNodeIdByCoarseSingleId(edge.first);
                int id2 = G1.findNodeIdByCoarseSingleId(i);
                if (id1 != -1 && id2 != -1) {
                    if (addedEdges.find({ id1, id2 }) == addedEdges.end() && addedEdges.find({ id2, id1 }) == addedEdges.end()) {
                        G1.setEdge(id1, id2, matAdj[edge.first][i][1]);
                        addedEdges.insert({ id1, id2 }); // Add the edge to the set
                    }
                }
            }
            else if (matAdj[edge.first][i][0] == 0 && matAdj[edge.second][i][0] == 1 && i != edge.first && i != edge.second) {
                int id1 = G1.findNodeIdByCoarseSingleId(edge.second);
                int id2 = G1.findNodeIdByCoarseSingleId(i);
                if (id1 != -1 && id2 != -1) {
                    if (addedEdges.find({ id1, id2 }) == addedEdges.end() && addedEdges.find({ id2, id1 }) == addedEdges.end()) {
                        G1.setEdge(id1, id2, matAdj[edge.second][i][1]);
                        addedEdges.insert({ id1, id2 }); // Add the edge to the set
                    }
                }
            }
        }
    }

    G1.setSizeNodes(G1.getNodes().size());
    G1.setSizeEdges(G1.getEdges().size());
    G1.computeAdjacencyMatrix();

    return G1;
}

void savePartitionDataToFile(const PartitionData& partitionData) {
    std::ofstream outputFile(partitionData.fileName);
    if (outputFile.is_open()) {
        outputFile << "Execution time graph reading: " << partitionData.executionTimes[0] << " seconds" << std::endl;
        outputFile << "Total nodes weight: " << partitionData.totalNodesWeight << std::endl;
        outputFile << "Total edges weight: " << partitionData.totalEdgesWeight << std::endl
            << std::endl;
        outputFile << "Execution time partitioning: " << partitionData.executionTimes[1] << " seconds" << std::endl
            << std::endl;
        for (size_t i = 0; i < partitionData.partitions.size(); i++) {
            if (partitionData.partitions[i].size() > 0) {
                outputFile << "Partition " << i + 1 << ": ";
                for (auto j : partitionData.partitions[i]) {
                    outputFile << j << " ";
                }
                outputFile << std::endl;
                outputFile << "Balance Factor: " << partitionData.balanceFactors[i] << " | Cut Size: " << partitionData.cutSizes[i] << std::endl
                    << std::endl;
            }
        }
        if (partitionData.averageBalanceFactor != 1)
            outputFile << "Average Balance Factor: " << partitionData.averageBalanceFactor << std::endl;
        else
            outputFile << "Average Balance Factor: " << partitionData.balanceFactors[0] << std::endl;
        outputFile << "Average Cut Size: " << partitionData.averageCutSize << std::endl;
        outputFile << "Cut Size Between Partitions: " << partitionData.cutSizePartitions << std::endl;
        outputFile << "CPU time used: " << partitionData.usage.ru_utime.tv_sec << " seconds " << partitionData.usage.ru_utime.tv_usec << " microseconds" << std::endl;
        outputFile << "Memory usage: " << partitionData.usage.ru_maxrss << " kilobytes | " << partitionData.usage.ru_maxrss / 1024.0 << " MBs | " << partitionData.usage.ru_maxrss / (1024.0 * 1024.0) << " GBs" << std::endl;
        outputFile << "CPU percentage: " << partitionData.cpu_percentage << "%" << std::endl;
        outputFile.close();
    }
}

double calculateAverageBalanceFactor(const std::vector<double>& balanceFactors) {
    double sum = 0.0;
    for (double balanceFactor : balanceFactors) {
        sum += balanceFactor;
    }
    return sum / balanceFactors.size();
}

double calculateAverageCutSize(const std::vector<int>& cutSizes) {
    double sum = 0.0;
    for (int cutSize : cutSizes) {
        sum += cutSize;
    }
    return sum / cutSizes.size();
}