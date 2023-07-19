#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <set>

typedef struct{

    int n;
    int weight;

} Node;

typedef struct{

    int n1;
    int n2;
    int weight;

} Edge;

class Graph {

    private:
        std::vector<Node> Nodes;
        std::vector<Edge> Edges;
        std::vector<std::vector<std::vector<int>>> Mat_Adj; //inside the adjacency matrix we have a field for the link e one for the weight
        int sizeN;
        int sizeE;

    public:
        int num_of_nodes(){return sizeN;}
        int num_of_edges(){return sizeE;}
        void setSizeNodes( int value ) { sizeN=value;}
        void setSizeEdges( int value ) { sizeE=value;}
        void setNode (int n, int weight){   
            Node value; 
            value.n=n; 
            value.weight=weight;  
            Nodes.push_back(value);
        }

        void setEdge (int n1, int n2, int weight){   
            Edge value; 
            value.n1=n1; 
            value.n2=n2;
            value.weight=weight;  
            Edges.push_back(value);
        }

        std::vector<Node> getNodes(){ return Nodes;}
        std::vector<Edge> getEdges(){ return Edges;}
        std::vector<std::vector<std::vector<int>>> getMatAdj(){ return Mat_Adj;}

        void computeAdjacencyMatrix() {
    // Clear the existing adjacency matrix
    Mat_Adj.clear();

    // Resize the adjacency matrix to match the number of nodes
    Mat_Adj.resize(sizeN, std::vector<std::vector<int>>(sizeN, std::vector<int>(2, 0)));

    // Iterate over the edges and update the adjacency matrix
    for (const Edge& edge : Edges) {
        int n1 = edge.n1;
        int n2 = edge.n2;
        int weight = edge.weight;

        // Set the link in the adjacency matrix
        Mat_Adj[n1][n2][0] = 1;
        Mat_Adj[n2][n1][0] = 1;

        // Set the weight in the adjacency matrix
        Mat_Adj[n1][n2][1] = weight;
        Mat_Adj[n2][n1][1] = weight;
    }
}

};


#endif