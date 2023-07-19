#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>

typedef struct{

    int weight;
    int degree;

} Node;

typedef struct{

    int n1;
    int n2;
    int weight;

} Edge;

class Graph {

    private:
        //the id of the map is the id of the node while the node struct contains info about the node
        std::map<int, Node> Nodes;
        std::vector<Edge> Edges;
        std::vector<std::vector<std::vector<int>>> MatAdj; //inside the adjacency matrix we have a field for the link e one for the weight
        std::vector<std::vector<int>> MatDegree;
        int sizeN;
        int sizeE;

    public:

        int num_of_nodes(){return sizeN;}
        int num_of_edges(){return sizeE;}
        void setSizeNodes( int value ) { sizeN=value;}
        void setSizeEdges( int value ) { sizeE=value;}
        
        std::map<int, Node> getNodes(){ return Nodes;}
        std::vector<Edge> getEdges(){ return Edges;}
        std::vector<std::vector<std::vector<int>>> getMatAdj(){ return MatAdj;}
        std::vector<std::vector<int>> getMatDegree(){ return MatDegree;}

        void setNode(int n, int weight);
        void setEdge(int n1, int n2, int weight);

        void computeAdjacencyMatrix();
        void computeMatrixDegree();
        void incrementDegree(int idNode);

       

};


#endif