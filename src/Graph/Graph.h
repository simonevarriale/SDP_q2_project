#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <stack>

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
        Graph(int numNodes = 0, int numEdges = 0) : sizeN(numNodes), sizeE(numEdges) {}
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

        //Debug functions
        void printNodes();
        void printEdges();
        void printAdjacencyMatrix();
        void printDegreeMatrix();


    bool isFullyConnected() {
        int sizeNodes = sizeN;
        std::vector<bool> visited(sizeNodes, false);

        // Perform DFS starting from each node
        for (int startNode = 0; startNode < sizeNodes; ++startNode) {
            if (!visited[startNode]) {
                std::stack<int> S;
                S.push(startNode);
                visited[startNode] = true;
                while (!S.empty()) {
                    int node = S.top();
                    S.pop();
                    for (const auto& neighbor : MatAdj[node]) {
                        int neighborNode = neighbor[0]; // Get the neighbor node ID
                        if (!visited[neighborNode]) {
                            visited[neighborNode] = true;
                            S.push(neighborNode);
                        }
                    }
                }
            }
        }

        // Check if all nodes are visited
        for (bool v : visited) {
            if (!v) {
                return false;
            }
        }
        return true;
    }

};


#endif