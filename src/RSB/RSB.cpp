#include "../Graph/Graph.h"
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/LU>
#include <cmath>

// Function to calculate the cut size between two sets A and B
extern int calculateCutSize(Graph& graph, const std::vector<bool>& partitionA);

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

std::vector<bool> RSB(Graph& G, int p) {

    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    // compute Laplacian matrix
    //std::cout << "Laplacian matrix:" << std::endl;
    for (int i = 0; i < sizeNodes; i++) {
        for (int j = 0; j < sizeNodes; j++) {
            L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
            //std::cout << L(i, j) << " ";
        }
        //std::cout << std::endl;
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    // Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real(); // ascending order
        //std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
        // Eigen::MatrixXd sortedEigenvectors = eigenvectors.rowwise().reverse(); //descending order
        // std::cout << "Sorted eigenvectors:\n" << sortedEigenvectors << std::endl;
        //std::cout << eigenvectors.col(1) << std::endl;
        // Find the Fiedler vector (second smallest eigenvector)
        Eigen::VectorXd fiedlerVector = eigenvectors.col(1);

        // Find the median value of the Fiedler vector
        double medianValue = computeMedian(fiedlerVector);
        std::vector<bool> partition(L.rows());
        for (int i = 0; i < L.rows(); ++i) {
            if (fiedlerVector(i) <= medianValue) {
                partition[i] = 0; // Assign node i to partition 0
            }
            else {
                partition[i] = 1; // Assign node i to partition 1
            }
        }

        double weightA = 0.0;
        double weightB = 0.0;
    for (int i = 0; i < G.num_of_nodes(); i++) {
        if (partition[i]) {
            weightA += G.getNodeWeight(i);
        }
        else {
            weightB += G.getNodeWeight(i);
        }
    }

     std::cout<<"Partition Balance Factor: " <<std::min(weightA, weightB) / std::max(weightA, weightB)<<std::endl;
     std::cout<<"Cut size RSB: " << calculateCutSize(G,partition) <<std::endl;


        return partition;

    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
    // Eigen::VectorXi sortedIndices = eigenvalues.argsort();
}


extern Graph coarsening(Graph G);

Eigen::VectorXd interpolate(Eigen::VectorXd fv1, Eigen::MatrixXd L, int sizeNodes){
    Eigen::VectorXd fv(sizeNodes);

    int sum=0;
    int num=0;

    for(int i=0; i<fv1.size(); i++){
        fv[i] = fv1[i];
    }

    for(int i=fv1.size(); i<sizeNodes; i++){

        for(int j=0; j<sizeNodes; j++){

            if(L(i,j) != 0 ){
                sum+=fv[j];
                num++;
            }

        }

        fv[i] = sum/num;
        sum=0;
        num=0;

    }

    return fv;
}

Eigen::VectorXd rqi(Eigen::VectorXd fv, Eigen::MatrixXd L, int sizeNodes){

    float theta = fv.transpose() * L * fv;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sizeNodes, sizeNodes);
    double p;

    do{
        Eigen::VectorXd x = (L-theta*I).lu().solve(fv);
        fv = x/x.norm();
        theta = fv.transpose() * L * fv;
        p = ((L*fv).transpose() * (L*fv) - theta*theta);
        p = sqrt(p);

    }while(p<0.1);

    return fv;

}



Eigen::VectorXd fiedler(Graph G){

    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);
    Eigen::VectorXd fv;
    Eigen::VectorXd fv1;

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    // compute Laplacian matrix
    //std::cout << "Laplacian matrix:" << std::endl;
    for (int i = 0; i < sizeNodes; i++) {
        for (int j = 0; j < sizeNodes; j++) {
            L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
            //std::cout << L(i, j) << " ";
        }
        //std::cout << std::endl;
    }
   

    if(sizeNodes > 10){ //grandezza grafo maggiore di un certo numero di nodi
        Graph G1 = coarsening(G);
        fv1 =  fiedler(G1);
        fv = interpolate(fv1, L, sizeNodes);
        fv = rqi(fv, L, sizeNodes);
        
    }
    else{
        
    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    // Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real(); // ascending order
        //std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
        // Eigen::MatrixXd sortedEigenvectors = eigenvectors.rowwise().reverse(); //descending order
        // std::cout << "Sorted eigenvectors:\n" << sortedEigenvectors << std::endl;
        //std::cout << eigenvectors.col(1) << std::endl;
        // Find the Fiedler vector (second smallest eigenvector)
            fv = eigenvectors.col(1);
        }

    }

    return fv;
}


std::vector<bool> MLRSB(Graph& G, int p){

    Eigen::VectorXd fiedlerV;

    fiedlerV = fiedler(G);

    double medianValue = computeMedian(fiedlerV);
        std::vector<bool> partition(G.num_of_nodes());
        for (int i = 0; i < G.num_of_nodes(); ++i) {
            if (fiedlerV(i) <= medianValue) {
                partition[i] = 0; // Assign node i to partition 0
            }
            else {
                partition[i] = 1; // Assign node i to partition 1
            }
        }

        double weightA = 0.0;
        double weightB = 0.0;
    for (int i = 0; i < G.num_of_nodes(); i++) {
        if (partition[i]) {
            weightA += G.getNodeWeight(i);
        }
        else {
            weightB += G.getNodeWeight(i);
        }
    }

     std::cout<<"Partition Balance Factor: " <<std::min(weightA, weightB) / std::max(weightA, weightB)<<std::endl;
     std::cout<<"Cut size RSB: " << calculateCutSize(G,partition) <<std::endl;


        return partition;


}