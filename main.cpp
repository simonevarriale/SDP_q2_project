#include <iostream>
#include <fstream>
#include "./Graph/Graph.h"

using namespace std;

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
    }

    inputFile.close();

}

int main() {
    
    Graph G;

    read_input("./test_graph.txt", &G);

    // Test the read input
    std::cout << "Number of nodes: " << G.num_of_nodes() << std::endl;
    std::cout << "Number of edges: " << G.num_of_edges() << std::endl;

    std::cout << "Nodes: ";
    for (auto&  n: G.getNodes()) {
        std::cout<<n.n<<" "<<n.weight<<endl;
    }
    std::cout << std::endl;

    std::cout << "Edges:" << std::endl;
    for (auto&  edge: G.getEdges()) {
        std::cout << edge.n1 << " " << edge.n2 << " " << edge.weight << std::endl;
    }

    G.computeAdjacencyMatrix();

    auto mat = G.getMatAdj();

    for(int i=0; i<G.num_of_nodes(); i++){
        for(int j=0; j<G.num_of_nodes(); j++){
            cout<<mat[i][j][0]<<"-"<<mat[i][j][1]<<" ";
        }
        cout<<endl;
    }

    return 0;
}


