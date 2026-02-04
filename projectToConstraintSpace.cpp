#include "projectToConstraintSpace.h"

//f_i(dofs) = ||v_i+1 - v_i||² = 1
//d_f_i(dofs) = (0,...,0, -2(v_i+1-v_i), 2(v_i+1 - v_i),0,...,0)
Eigen::SparseMatrix<double> springJacobian(Eigen::VectorXd& dofs){
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t i = 0; i < dofs.size()-3; i+=3){
        Eigen::Vector3d v1 =  dofs.segment<3>(i);
        Eigen::Vector3d v2 =  dofs.segment<3>(i+3);
        Eigen::Vector3d e = v2-v1;
        for( size_t k = 0; k< 3; ++k){
            triplets.emplace_back(i, i + k, -2*e(k));
            triplets.emplace_back(i, i + k + 3, 2*e(k));
        }
    }
    //edge -1 0
    size_t id_last = dofs.size()-4;
    Eigen::Vector3d v1 =  dofs.segment<3>(id_last);   //last
    Eigen::Vector3d v2 =  dofs.segment<3>(0);         //first
    Eigen::Vector3d e = v2-v1;
    for( size_t k = 0; k< 3; ++k){
        triplets.emplace_back(dofs.size()-1, id_last + k, -2*e(k));
        triplets.emplace_back(dofs.size()-1, 0 + k + 3, 2*e(k));
    }

    Eigen::SparseMatrix<double> J(dofs.size(),dofs.size());
    J.setFromTriplets(triplets.begin(), triplets.end());
    

    return J;    
}

//delta = d0​−XT (XXT)−1 X*d0  
//d0 disired direction
//X is constraints
//returns direction that preseves constraints
Eigen::VectorXd projectToTangentSpace(
    const Eigen::SparseMatrix<double>& X,  // (m x n)
    const Eigen::VectorXd& d0               // (n)
){
    // v = X * d0
    Eigen::VectorXd v = X * d0;

    // A = X * Xᵀ
    Eigen::SparseMatrix<double> A = X * X.transpose();

    // Solve A * lambda = v
    Eigen::VectorXd lambda;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    lambda = solver.solve(v);

    // delta = d0 - Xᵀ * lambda   (n)
    Eigen::VectorXd delta = d0 - X.transpose() * lambda;

    return delta;
}


