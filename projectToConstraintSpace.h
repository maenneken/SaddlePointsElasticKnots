#pragma once 
#include "helpers.h"
Eigen::SparseMatrix<double> stackJacobians(const Eigen::SparseMatrix<double>& J1, const Eigen::SparseMatrix<double>& J2);
Eigen::SparseMatrix<double> stretchJacobian(Eigen::VectorXd& dofs);
Eigen::SparseMatrix<double> bendJacobian(Eigen::VectorXd& dofs);
Eigen::SparseMatrix<double> contactJacobian(Eigen::MatrixXd& contactForces);
Eigen::SparseMatrix<double> twistJacobian(Eigen::VectorXd& dofs);
Eigen::SparseMatrix<double> gradientJacobian(Eigen::VectorXd& gradient);
Eigen::VectorXd projectToTangentSpace(const Eigen::SparseMatrix<double>& X, const Eigen::VectorXd& d0);