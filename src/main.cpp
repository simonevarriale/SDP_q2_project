#include <iostream>
#include <fstream>
#include <random>
#include "./Graph/Graph.h"

using namespace std;
extern void RSB(Graph *G, int p);
extern void kernighanLin(Graph& graph);
extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);
extern void computeNetGains(Graph& graph, const std::vector<bool>& partitionA, std::vector<int>& netGains);

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

void read_input(const std::string& filename , Graph* G){

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
        G->setNode(n,weight);
    }

    int source;
    int destination;
    // Read the edges
    for (int i = 0; i < numEdges; ++i) {
        inputFile >> source >> destination >> weight;
        G->setEdge(source, destination,weight);
        G->incrementDegree(source);
        G->incrementDegree(destination);
    }

    inputFile.close();

}

int main() {
    
    Graph G;

    read_input("./simple_graph.txt", &G);

    // Test the read input
    std::cout << "Number of nodes: " << G.num_of_nodes() << std::endl;
    std::cout << "Number of edges: " << G.num_of_edges() << std::endl;

    G.printNodes();
    // std::cout << "Nodes: " << std::endl;
    // for (auto&  n: G.getNodes()) {
    //     std::cout<<n.first<<" "<<n.second.weight<<endl;
    // }
    // std::cout << std::endl;

    G.printEdges();
    // std::cout << "Edges:" << std::endl;
    // for (auto&  edge: G.getEdges()) {
    //     std::cout << edge.n1 << " " << edge.n2 << " " << edge.weight << std::endl;
    // }

    G.computeAdjacencyMatrix();

    auto mat = G.getMatAdj();

    G.printAdjacencyMatrix();
    // for(int i=0; i<G.num_of_nodes(); i++){
    //     for(int j=0; j<G.num_of_nodes(); j++){
    //         cout<<mat[i][j][0]<<"-"<<mat[i][j][1]<<" ";
    //     }
    //     cout<<endl;
    // }

    G.printNodes();
    // std::cout << "Nodes: " << std::endl;
    // for (auto&  n: G.getNodes()) {
    //     std::cout<<n.first<<" "<<n.second.degree<<endl;
    // }

    G.computeMatrixDegree();

    auto mat1 = G.getMatDegree();
    cout << "Degree Matrix:" << endl;
    for(int i=0; i<G.num_of_nodes(); i++){
        for(int j=0; j<G.num_of_nodes(); j++){
            cout<<mat1[i][j];
        }
        cout<<endl;
    }

    auto partition = generateRandomBisection(G);
    std::cout << "Initial partition: " << std::endl;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        std::cout << partition[i] << " ";
    }
    std::cout << std::endl;

    // RSB(&G, 2);
    std::vector<int> netGains(G.num_of_nodes(), 0); // Net gains for each node
    computeNetGains(G, partition, netGains);
    std::cout << "Net gains: " << std::endl;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        std::cout << netGains[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}