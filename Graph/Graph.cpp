#include "Graph.h"

void Graph::setNode(int n, int weight){   
    Node value; 
    value.weight=weight;  
    value.degree=0;
    Nodes.insert({n,value});
}

void Graph::setEdge(int n1, int n2, int weight){   
    Edge value; 
    value.n1=n1; 
    value.n2=n2;
    value.weight=weight;  
    Edges.push_back(value);
}


void Graph::computeAdjacencyMatrix() {
    // Clear the existing adjacency matrix
    MatAdj.clear();

    // Resize the adjacency matrix to match the number of nodes
    MatAdj.resize(sizeN, std::vector<std::vector<int>>(sizeN, std::vector<int>(2, 0)));

    // Iterate over the edges and update the adjacency matrix
    for (const Edge& edge : Edges) {
        int n1 = edge.n1;
        int n2 = edge.n2;
        int weight = edge.weight;

        // Set the link in the adjacency matrix
        MatAdj[n1][n2][0] = 1;
        MatAdj[n2][n1][0] = 1;

        // Set the weight in the adjacency matrix
        MatAdj[n1][n2][1] = weight;
        MatAdj[n2][n1][1] = weight;
    }
}

 void Graph::computeMatrixDegree(){
    MatDegree.resize(sizeN, std::vector<int>(sizeN,(0)));
    for(int i=0; i<sizeN; i++){
        MatDegree[i][i] = Nodes.find(i)->second.degree;
    }
}

void Graph::incrementDegree(int idNode){
    auto it = Nodes.find(idNode);
    it->second.degree++;
}