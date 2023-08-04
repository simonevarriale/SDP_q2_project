#include <iostream>
#include <fstream>
#include <random>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <climits>
#include "./Graph/Graph.h"

using namespace std;
extern void RSB(Graph& G, int p);
//extern std::vector<bool> kernighanLin(Graph& graph);
extern std::vector<bool> kernighanLin(Graph& graph, std::vector<bool> partitionA = {});
extern Graph coarsening(Graph graph);
extern std::vector<bool> uncoarsening(Graph G1, std::vector<bool> partition, int graphSize);
extern std::vector<bool> multilevel_KL(Graph& graph);
extern std::vector<bool> fiducciaMattheyses(Graph& graph, int maxIterations);
extern std::vector<bool> fm(Graph& graph);

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

int main() {

    Graph G;
    std::string file1 = "./simple_graph.txt";
    std::string file2 = "./test_graph.txt";
    std::string file3 = "./connected_graph.txt";
    std::string file4 = "./big_graph.txt";
    read_input(file4, &G);

    // Test the read input
    // std::cout << "Number of nodes: " << G.num_of_nodes() << std::endl;
    // std::cout << "Number of edges: " << G.num_of_edges() << std::endl;

    // G.printNodes();

    // G.printEdges();

    // G.printAdjacencyMatrix();
    // G.computeMatrixDegree();

    // std::vector<bool> partitionA = kernighanLin(G);
    // std::cout << "Final partition: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << partitionA[i] << " ";
    // }
    // std::cout << std::endl;

    // Graph G1 = coarsening(G);

    // G1.printGraph();

    // std::vector<bool> partitionB = kernighanLin(G1);
    // std::cout << "Final partition KL: " << std::endl;
    // for (int i = 0; i < G1.num_of_nodes(); ++i) {
    //     std::cout << partitionB[i] << " ";
    // }
    // std::cout << std::endl;

    // Graph G2 = uncoarsening(G1);
    // std::cout << "Original Graph" << std::endl;
    //G.printEdges();
    // G.printGraph();
    // std::cout << "Coarsened Graph" << std::endl;
    // G1.printGraph();
    // std::cout << "Uncoarsened Graph" << std::endl;
    // G2.printGraph();
    //G2.printEdges();
    // G2.printAdjacencyMatrix();
    // G.printAdjacencyMatrix();

    // std::vector<bool> fm = fiducciaMattheyses(G, 2);
    // std::cout << "Fm partition: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << fm[i] << " ";
    // }
    // std::cout << std::endl;

    // std::vector<bool> part = multilevel_KL(G);
    // std::cout << "Final partition uncoarsen: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << part[i] << " ";
    // }
    // std::cout << std::endl;

    //G.printGraph();
    // Graph G1 = coarsening(G);
    // auto part = kernighanLin(G1);
    // G1.printGraph();
    // std::cout << "Final partition coarsen: " << std::endl;
    // for (int i = 0; i < G1.num_of_nodes(); ++i) {
    //     std::cout << part[i] << " ";
    // }
    // std::cout << std::endl;
    // auto uncors = uncoarsening(G1, part, G.num_of_nodes());
    // std::cout << "Final partition uncoarsen: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << uncors[i] << " ";
    // }
    
    std::cout << std::endl;
    auto klPart = kernighanLin(G);
    std::cout << "Final partition KL: " << std::endl;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        std::cout << klPart[i] << " ";
    }
    std::cout << std::endl;

    auto multilevel = multilevel_KL(G);
    std::cout << "Final multilevel partition KL: " << std::endl;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
         std::cout << multilevel[i] << " ";
    }
    std::cout << std::endl;

    RSB(G,2);
    


    // auto fmPart = fm(G);
    // std::cout << "Final partition FM: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << fmPart[i] << " ";
    // }
    // std::cout << std::endl;

    // auto secondFm = fiducciaMattheyses(G, 2);
    // std::cout << "Final partition FM OG: " << std::endl;
    // for (int i = 0; i < G.num_of_nodes(); ++i) {
    //     std::cout << secondFm[i] << " ";
    // }
    // std::cout << std::endl;


    return 0;
}