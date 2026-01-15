#include "helpers.h"

double dist(Eigen::Vector3d xi, Eigen::Vector3d xj){
    return (xi - xj).norm();
}


double E(Eigen::Vector3d xi, Eigen::Vector3d xj){
    double d = dist(xi,xj);
    
    return 1 / ((d * d) + 1);
}
//dxiE(xi,xj) = - dxjE(xi,xj)
Eigen::Vector3d dE(Eigen::Vector3d xi, Eigen::Vector3d xj){
    Eigen::Vector3d r = xi - xj;
    double s = r.squaredNorm();
    return -2.0 * r / ((s+1)*(s+1));
}
//this is Hii Hij = -Hii; Hji = - Hii; Hjj = Hii;
Eigen::Matrix3d H(Eigen::Vector3d xi, Eigen::Vector3d xj){
    Eigen::Vector3d r = xi - xj;
    double s = r.squaredNorm();
    double denom1 = s + 1.0;

    // Hessian scalars
    double a = -4.0 / (denom1 * denom1);
    double b =  8.0 / (denom1 * denom1 * denom1);

    // 3x3 Hessian block
    Eigen::Matrix3d H = a * Eigen::Matrix3d::Identity() + b * (r * r.transpose());
    return H;
}



//minSeparation defines how far a pt should be away for it to be calculedted in the energy. (ie minSeparation = 0 -> all neighbors are looked at)
double distanceEnergy(Eigen::VectorXd DoFs, int minSeparation){
    double energy = 0;
    int n_pts = DoFs.size() * 0.25;
    std::vector<Eigen::Vector3d> pts = DoFsToPos(DoFs,n_pts);
    for (int i = 0; i <  n_pts; ++i){
        for ( int j = i+1; j <  n_pts; j++){
            int topoDist = std::min(j - i,  n_pts - (j - i));
            if (topoDist < minSeparation) continue;

            energy += E(pts[i],pts[j]);
        }
    }
    return energy;
}

Eigen::VectorXd distanceGradient(Eigen::VectorXd DoFs, int minSeparation){
    Eigen::VectorXd g = Eigen::VectorXd::Zero(DoFs.size());
    int n_pts = DoFs.size() * 0.25;
    std::vector<Eigen::Vector3d> pts = DoFsToPos(DoFs,n_pts);
    for (int i = 0; i <  n_pts; ++i){
        for ( int j = i+1; j <  n_pts; j++){
            int topoDist = std::min(j - i,  n_pts - (j - i));
            if (topoDist < minSeparation) continue;

            Eigen::Vector3d pt_grad = dE(pts[i],pts[j]);

            g.segment<3>(3*i) += pt_grad;
            g.segment<3>(3*j) -= pt_grad;

        }
    }
    return g;
}

Eigen::SparseMatrix<double> distanceHessian(Eigen::VectorXd DoFs, int minSeparation){
    typedef Eigen::Triplet<double> T;
    std::vector<Eigen::Triplet<double>> triplets;
    
    int n_pts = DoFs.size() * 0.25;
    std::vector<Eigen::Vector3d> pts = DoFsToPos(DoFs,n_pts);
    for (int i = 0; i <  n_pts; ++i){
        for ( int j = i+1; j <  n_pts; j++){
            int topoDist = std::min(j - i,  n_pts - (j - i));
            if (topoDist < minSeparation) continue;

            Eigen::Matrix3d H_block = H(pts[i],pts[j]);
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) {
                    triplets.push_back(T(3*i+a, 3*i+b, H_block(a,b)));    // Hii
                    triplets.push_back(T(3*j+a, 3*j+b, H_block(a,b)));    // Hjj
                    triplets.push_back(T(3*i+a, 3*j+b, -H_block(a,b)));   // Hij
                    triplets.push_back(T(3*j+a, 3*i+b, -H_block(a,b)));   // Hji
    }

        }
    }
    Eigen::SparseMatrix<double> H_sparse(DoFs.size(),DoFs.size());
    H_sparse.setFromTriplets(triplets.begin(), triplets.end());

    return H_sparse;
}