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
extern Graph coarsening(Graph& G);
bool hasImproved = true;

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

        std::cout << "Partition Balance Factor RSB: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
        std::cout << "Cut size RSB: " << calculateCutSize(G, partition) << std::endl;


        return partition;

    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
    // Eigen::VectorXi sortedIndices = eigenvalues.argsort();
}

Eigen::VectorXd interpolate(Eigen::VectorXd fv1, Eigen::MatrixXd L, int sizeNodes) {
    // std::cout << "Interpolate" << std::endl;
    Eigen::VectorXd fv(sizeNodes);

    int sum = 0;
    int num = 0;

    for (int i = 0; i < fv1.size(); i++) {
        fv[i] = fv1[i];
    }

    for (int i = fv1.size(); i < sizeNodes; i++) {

        for (int j = 0; j < sizeNodes; j++) {

            if (L(i, j) != 0) {
                sum += fv[j];
                num++;
            }

        }

        fv[i] = sum / num;
        sum = 0;
        num = 0;

    }

    return fv;
}

Eigen::VectorXd rqi(Eigen::VectorXd fv, Eigen::MatrixXd L, int sizeNodes) {
    // std::cout << "RQI" << std::endl;
    float theta = fv.transpose() * L * fv;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sizeNodes, sizeNodes);
    double p;
    do {
        Eigen::VectorXd x = (L - theta * I).lu().solve(fv);
        fv = x / x.norm();
        theta = fv.transpose() * L * fv;
        p = ((L * fv).transpose() * (L * fv) - theta * theta);
        p = sqrt(p);
        // std::cout << "p: " << p << std::endl;
    } while (p < 0.0000001);

    return fv;

}

Eigen::VectorXd fiedler(Graph& G) {

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

    if (sizeNodes > 50 && hasImproved) { //grandezza grafo maggiore di un certo numero di nodi
        // std::cout << "Coarsening" << std::endl;
        Graph G1 = coarsening(G);
        if (G1.num_of_nodes() == G.num_of_nodes() - 1)
            hasImproved = false;
        // std::cout << "Fine Coarsening con " << G1.num_of_nodes() << " nodes." << std::endl;
        fv1 = fiedler(G1);
        fv = interpolate(fv1, L, sizeNodes);
        fv = rqi(fv, L, sizeNodes);
        // std::cout << "Fine RQI" << std::endl;
    }
    else {
        // std::cout << "Laplacian matrix:" << std::endl;
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

std::vector<bool> MLRSB(Graph& G, int p) {

    Eigen::VectorXd fiedlerV;

    fiedlerV = fiedler(G);

    double medianValue = computeMedian(fiedlerV);
    std::vector<bool> partition(G.num_of_nodes());
    // std::cout << "Start partitioning" << std::endl;
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        if (fiedlerV(i) <= medianValue) {
            partition[i] = 0; // Assign node i to partition 0
        }
        else {
            partition[i] = 1; // Assign node i to partition 1
        }
    }
    // std::cout << "End partitioning" << std::endl;

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

    std::cout << "Partition Balance Factor MLRSB: " << std::min(weightA, weightB) / std::max(weightA, weightB) << std::endl;
    std::cout << "Cut size MLRSB: " << calculateCutSize(G, partition) << std::endl;


    return partition;


}

std::vector<int> sortIndices(const Eigen::VectorXd& vec) {
    std::vector<int> indices(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        indices[i] = i;
    }

    // Sort indices based on the values in the vector
    std::sort(indices.begin(), indices.end(), [&vec](int a, int b) {
        return vec[a] < vec[b];
        });

    return indices;
}

int calculateCutSizeBetweenPartitions(Graph& G, const std::vector<bool>& partitionA, const std::vector<bool>& partitionB) {
    int cutSize = 0;

    const std::vector<std::vector<std::vector<int>>>& matAdj = G.getMatAdj();

    for (int i = 0; i < G.num_of_nodes(); ++i) {
        for (int j = 0; j < G.num_of_nodes(); ++j) {
            if (partitionA[i] && partitionB[j] && matAdj[i][j][0] != -1) {
                cutSize += matAdj[i][j][1]; // Use the weight stored in the adjacency matrix
            }
        }
    }

    return cutSize;
}


// Calculate the total cut size across multiple partitions
int calculateCutSizePartitions(Graph& G, const std::vector<std::vector<bool>>& partitions) {
    int totalCutSize = 0;
    for (int i = 0; i < partitions.size(); ++i) {
        for (int j = i + 1; j < partitions.size(); ++j) {
            totalCutSize += calculateCutSizeBetweenPartitions(G, partitions[i], partitions[j]);
        }
    }
    return totalCutSize;
}

std::vector<std::vector<bool>> pWayPartition(Graph& G, int p) {
    int sizeNodes = G.num_of_nodes();
    Eigen::MatrixXd L(sizeNodes, sizeNodes);

    G.computeMatrixDegree();
    auto matDeg = G.getMatDegree();
    auto matAdj = G.getMatAdj();

    // Compute Laplacian matrix
    for (int i = 0; i < sizeNodes; i++) {
        for (int j = 0; j < sizeNodes; j++) {
            L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

    if (eigenSolver.info() == Eigen::Success) {
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();
        Eigen::VectorXd fiedlerVector = eigenvectors.col(1);

        // Sort fiedlerVector to get partitioning indices
        std::cout << "Sorting indices" << std::endl;
        std::vector<int> sortedIndices = sortIndices(fiedlerVector);
        std::cout << "Finished sorting indices" << std::endl;

        // Calculate the size of each sub-partition
        int subPartitionSize = sizeNodes / p;
        std::vector<std::vector<bool>> partitions(p, std::vector<bool>(sizeNodes, false));

        for (int i = 0; i < p; ++i) {
            int startIndex = i * subPartitionSize;
            int endIndex = (i == p - 1) ? sizeNodes : (i + 1) * subPartitionSize;

            for (int j = startIndex; j < endIndex; ++j) {
                partitions[i][sortedIndices[j]] = true;
            }
        }

        std::vector<double> partitionWeights(p, 0.0);
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < G.num_of_nodes(); ++j) {
                if (partitions[i][j]) {
                    partitionWeights[i] += G.getNodeWeight(j);
                }
            }
        }

        for (int i = 0; i < p; ++i) {
            std::cout << "Partition " << i << " weight: " << partitionWeights[i] << std::endl;
        }

        std::cout << "Cut size pRSB: " << calculateCutSizePartitions(G, partitions) << std::endl;

        return partitions;
    }
    else {
        std::cout << "Failed to compute eigenvectors." << std::endl;
    }

    return {};
}

std::vector<std::vector<bool>> pMLRSB(Graph& G, int p) {

    Eigen::VectorXd fiedlerV;

    fiedlerV = fiedler(G);

    std::vector<double> partitionValues(G.num_of_nodes());

    for (int i = 0; i < G.num_of_nodes(); ++i) {
        partitionValues[i] = fiedlerV(i);
    }

    std::vector<double> partitionThresholds(p - 1);
    std::vector<std::vector<bool>> partitions(p, std::vector<bool>(G.num_of_nodes(), false));

    // Sort the partition values
    std::sort(partitionValues.begin(), partitionValues.end());

    // Calculate partition thresholds
    for (int i = 1; i < p; ++i) {
        partitionThresholds[i - 1] = partitionValues[i * G.num_of_nodes() / p];
    }

    // Assign nodes to partitions based on thresholds
    for (int i = 0; i < G.num_of_nodes(); ++i) {
        bool assigned = false;
        for (std::size_t j = 0; j < p - 1; ++j) {
            if (fiedlerV(i) <= partitionThresholds[j]) {
                partitions[j][i] = true;
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            partitions[p - 1][i] = true;
        }
    }

    // Calculate partition statistics and print results
    for (std::size_t part = 0; part < p; ++part) {
        double partitionWeight = 0.0;
        int numNodes = 0;
        for (int i = 0; i < G.num_of_nodes(); ++i) {
            if (partitions[part][i]) {
                partitionWeight += G.getNodeWeight(i);
                numNodes++;
            }
        }
        std::cout << "Partition " << part << " Weight: " << partitionWeight << "; Number of nodes: " << numNodes << std::endl;
    }

    std::cout << "Cut size pMLRSB: " << calculateCutSizePartitions(G, partitions) << std::endl;

    return partitions;
}
