#include "MLRSB.h"

bool hasImproved = true;

Eigen::VectorXd interpolate(Eigen::VectorXd fv1, Eigen::MatrixXd L, int sizeNodes) {
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
    float theta = fv.transpose() * L * fv;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sizeNodes, sizeNodes);
    double p;
    do {
        Eigen::VectorXd x = (L - theta * I).lu().solve(fv);
        fv = x / x.norm();
        theta = fv.transpose() * L * fv;
        p = ((L * fv).transpose() * (L * fv) - theta * theta);
        p = sqrt(p);
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
    for (int i = 0; i < sizeNodes; i++) {
        for (int j = 0; j < sizeNodes; j++) {
            L(i, j) = matDeg[i][j] - matAdj[i][j][0] * matAdj[i][j][1];
        }
    }

    if (sizeNodes > 50 && hasImproved) { //grandezza grafo maggiore di un certo numero di nodi
        Graph G1 = coarsening(G);
        if (G1.num_of_nodes() == G.num_of_nodes() - 1)
            hasImproved = false;
        fv1 = fiedler(G1);
        fv = interpolate(fv1, L, sizeNodes);
        fv = rqi(fv, L, sizeNodes);
    }
    else {
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(L);

        if (eigenSolver.info() == Eigen::Success) {
            Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real(); // ascending order

            // Find the Fiedler vector (second smallest eigenvector)
            fv = eigenvectors.col(1);
        }

    }

    return fv;
}

std::vector<bool> MLRSB(Graph& G) {

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

    return partition;

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

    return partitions;
}
