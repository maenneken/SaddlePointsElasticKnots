#pragma once 
#include "helpers.h"

Eigen::SparseMatrix<double> springJacobian(Eigen::VectorXd& dofs);
Eigen::VectorXd projectToTangentSpace(const Eigen::SparseMatrix<double>& X, const Eigen::VectorXd& d0);