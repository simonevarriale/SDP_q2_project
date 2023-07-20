#include "../Graph/Graph.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

double computeMedian(const Eigen::VectorXd& vector) {
    Eigen::VectorXd sortedVector = vector;
    std::sort(sortedVector.data(), sortedVector.data() + sortedVector.size());
    int size = sortedVector.size();

    if (size % 2 == 0) {
        // If the vector size is even, take the average of the middle two values
        int mid = size / 2;
        return 0.5 * (sortedVector(mid - 1) + sortedVector(mid));
    } else {
        // If the vector size is odd, return the middle value
        return sortedVector(size / 2);
    }
}

void RSB(Graph *G, int p){

    int sizeNodes = G->num_of_nodes();
    Eigen::MatrixXd L(sizeNodes,sizeNodes);

    auto matDeg = G->getMatDegree();
    auto matAdj = G->getMatAdj();

    //compute Laplacian matrix
    std::cout << "Laplacian matrix:" << std::endl;
    for(int i=0; i<sizeNodes; i++){
        for(int j=0; j<sizeNodes; j++){
            L(i,j) = matDeg[i][j] - matAdj[i][j][0]*matAdj[i][j][1];
            std::cout << L(i,j) << " ";
        }
        std::cout << std::endl;
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    // Eigen::VectorXd eigenvalues = solver.eigenvalues().real();


    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real(); //ascending order
        std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
        // Eigen::MatrixXd sortedEigenvectors = eigenvectors.rowwise().reverse(); //descending order
        // std::cout << "Sorted eigenvectors:\n" << sortedEigenvectors << std::endl;
        std::cout << eigenvectors.col(1) << std::endl;
        // Find the Fiedler vector (second smallest eigenvector)
        Eigen::VectorXd fiedlerVector = eigenvectors.col(1);

        // Find the median value of the Fiedler vector
        double medianValue = computeMedian(fiedlerVector);
        std::vector<int> partition(L.rows());
        for (int i = 0; i < L.rows(); ++i) {
            if (fiedlerVector(i) <= medianValue) {
                partition[i] = 0; // Assign node i to partition 0
            } else {
                partition[i] = 1; // Assign node i to partition 1
            }
        }

        // Print the partitioning result
        std::cout << "Partitioning result:\n";
        for (int i = 0; i < L.rows(); ++i) {
            std::cout << "Node " << i << " in Partition " << partition[i] << "\n";
        }
    } else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    // Eigen::VectorXi sortedIndices = eigenvalues.argsort();
    
}