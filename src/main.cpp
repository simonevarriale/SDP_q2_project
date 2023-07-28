#include <iostream>
#include <fstream>
#include <random>
#include "./Graph/Graph.h"

using namespace std;
// extern void RSB(Graph *G, int p);
extern std::vector<bool> kernighanLin(Graph& graph);
extern Graph coarsening(Graph graph);
extern Graph uncoarsening(Graph G1);
extern std::vector<bool> multilevel_KL(Graph& graph);

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

    inputFile.close();
}

int main() {

    Graph G;
    std::string file1 = "./simple_graph.txt";
    std::string file2 = "./test_graph.txt";
    std::string file3 = "./connected_graph.txt";
    read_input(file1, &G);

    // Test the read input
    std::cout << "Number of nodes: " << G.num_of_nodes() << std::endl;
    std::cout << "Number of edges: " << G.num_of_edges() << std::endl;

    G.printNodes();

    G.printEdges();

    G.computeAdjacencyMatrix();

    G.printAdjacencyMatrix();
    // G.computeMatrixDegree();

    std::vector<bool> partitionA = kernighanLin(G);

    std::cout << "Final partition: " << std::endl;

    for (int i = 0; i < G.num_of_nodes(); ++i) {
        std::cout << partitionA[i] << " ";
    }
    std::cout << std::endl;

    Graph G1 = coarsening(G);

    G1.printGraph();

    std::vector<bool> partitionB = kernighanLin(G1);
    std::cout << "Final partition: " << std::endl;
    for (int i = 0; i < G1.num_of_nodes(); ++i) {
        std::cout << partitionB[i] << " ";
    }
    std::cout << std::endl;

    Graph G2 = uncoarsening(G1);
    std::cout << "Original Graph" << std::endl;
    //G.printEdges();
    G.printGraph();
    std::cout << "Coarsened Graph" << std::endl;
    G1.printGraph();
    std::cout << "Uncoarsened Graph" << std::endl;
    G2.printGraph();
    //G2.printEdges();
    // G2.printAdjacencyMatrix();
    // G.printAdjacencyMatrix();

    std::vector<bool> part = multilevel_KL(G);
    std::cout << "Final partition uncoarsen: " << std::endl;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        std::cout << part[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}