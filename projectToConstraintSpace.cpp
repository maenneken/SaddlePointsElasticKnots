#include "projectToConstraintSpace.h"


Eigen::SparseMatrix<double> stackJacobians(
    const Eigen::SparseMatrix<double>& J1,
    const Eigen::SparseMatrix<double>& J2)
{
    // Dimensions
    const int m1 = J1.rows();
    const int m2 = J2.rows();
    const int n  = J1.cols();   // must match J2.cols()

    assert(J2.cols() == n && "Jacobians must have same number of columns (DOFs)");

    Eigen::SparseMatrix<double> J(m1 + m2, n);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(J1.nonZeros() + J2.nonZeros());

    // Copy J1
    for (int k = 0; k < J1.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(J1, k); it; ++it)
            triplets.emplace_back(it.row(), it.col(), it.value());

    // Copy J2 (row-shifted)
    for (int k = 0; k < J2.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(J2, k); it; ++it)
            triplets.emplace_back(it.row() + m1, it.col(), it.value());

    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

//f_i(dofs) = ||v_i+1 - v_i||² = 1
//d_f_i(dofs) = (0,...,0, -2(v_i+1-v_i), 2(v_i+1 - v_i),0,...,0)
Eigen::SparseMatrix<double> stretchJacobian(Eigen::VectorXd& dofs){
    std::vector<Eigen::Triplet<double>> triplets;
    size_t n_edges = dofs.size() / 3;
    for (size_t i = 0; i < n_edges; ++i){
        size_t i1 = (i+1) % n_edges; //wrap to start
        Eigen::Vector3d v1 =  dofs.segment<3>(3*i);
        Eigen::Vector3d v2 =  dofs.segment<3>(3*i1); 
        Eigen::Vector3d e = v2-v1;
        for( size_t k = 0; k< 3; ++k){
            triplets.emplace_back(i, 3*i + k, -2*e(k));
            triplets.emplace_back(i, 3*i1 + k, 2*e(k));
        }
    }
    Eigen::SparseMatrix<double> J(n_edges,dofs.size());
    J.setFromTriplets(triplets.begin(), triplets.end());
    

    return J;    
}
//fi​=∥(vi+1​−vi​)−(vi​−vi−1​)∥² = fi​=∥vi+1​−2vi​+vi−1​∥²
//c_i = vi+1​−2vi​+vi−1
//∂f/∂v_{i-1} = +2 c_i
//∂f/∂v_i = -4 c_i
//∂f/∂v_{i+1} = +2 c_i
//wire wants to be straight
//todo e_i dot e_i+1 ableiten sollte konstante bend an jedem punkt
//abgeleitet aus discrete curvature aus discrete Elastic Rod paper
Eigen::SparseMatrix<double> bendJacobian(Eigen::VectorXd& dofs){
    std::vector<Eigen::Triplet<double>> triplets;
    size_t n_edges = dofs.size() / 3;
    for (size_t i = 0; i < n_edges; ++i){
        size_t im1 = (i + n_edges - 1) % n_edges; //wrap to start
        size_t i1 = (i+1) % n_edges; //wrap to start
        
        Eigen::Vector3d v0 =  dofs.segment<3>(3*im1);
        Eigen::Vector3d v1 =  dofs.segment<3>(3*i);
        Eigen::Vector3d v2 =  dofs.segment<3>(3*i1); 
        Eigen::Vector3d c= v2 - 2*v1 + v0;
        for( size_t k = 0; k< 3; ++k){
            triplets.emplace_back(i, 3*im1 + k, 2*c(k));
            triplets.emplace_back(i, 3*i + k, -4*c(k));
            triplets.emplace_back(i, 3*i1 + k, 2*c(k));
        }
    }
    Eigen::SparseMatrix<double> J(n_edges,dofs.size());
    J.setFromTriplets(triplets.begin(), triplets.end());
    
    return J;    
}
//fi​(dofs)=(vi+2​−2vi+1​+vi​)−(vi+1​−2vi​+vi−1​)
//fi​=vi+2​−3vi+1​+3vi​−vi−1
//curvature between neighbors should stay the same
//not realy twist it is more a bend

Eigen::SparseMatrix<double> twistJacobian(Eigen::VectorXd& dofs){
    std::vector<Eigen::Triplet<double>> triplets;
    size_t n_edges = dofs.size() / 3;
    for (size_t i = 0; i < n_edges; ++i){
        size_t im1 = (i + n_edges - 1) % n_edges; //wrap to start
        size_t i1 = (i+1) % n_edges; //wrap to start
        size_t i2 = (i+2) % n_edges; //wrap to start
        
        for( size_t k = 0; k< 3; ++k){
            triplets.emplace_back(i, 3*im1 + k,     -1.0);
            triplets.emplace_back(i, 3*i   + k,      3.0);
            triplets.emplace_back(i, 3*i1  + k,     -3.0);
            triplets.emplace_back(i, 3*i2  + k,      1.0);
        }
    }
    Eigen::SparseMatrix<double> J(n_edges,dofs.size());
    J.setFromTriplets(triplets.begin(), triplets.end());
    
    return J;    
}

//twist bend stretch together only lets the knot be rotated... maybe weighted least squares?

//contactforce has shape n_vertecies x 3
Eigen::SparseMatrix<double> contactJacobian(Eigen::MatrixXd& contactForces){
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t i = 0; i < contactForces.rows(); ++i){
        Eigen::Vector3d v = contactForces.row(i).transpose();
        for( size_t k = 0; k< 3; ++k){
            triplets.emplace_back(0, 3*i + k, v(k));
        }
    }
    Eigen::SparseMatrix<double> J(1,contactForces.rows()*3);
    J.setFromTriplets(triplets.begin(), triplets.end());

    return J;    
}
Eigen::SparseMatrix<double> gradientJacobian(Eigen::VectorXd& gradient){
    std::vector<Eigen::Triplet<double>> triplets;
    size_t n_vertecies = gradient.size() / 3;

    for (size_t i = 0; i < n_vertecies; ++i){
        Eigen::Vector3d v = gradient.segment<3>(3*i);
        for( size_t k = 0; k< 3; ++k){
            triplets.emplace_back(0, 3*i + k, v(k));
        }
    }
    Eigen::SparseMatrix<double> J(1,n_vertecies*3);
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


