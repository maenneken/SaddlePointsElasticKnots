#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

// ElasticRods
#include "3rdparty/ElasticKnots/3rdparty/ElasticRods/ElasticRod.hh"
#include "3rdparty/ElasticKnots/3rdparty/ElasticRods/PeriodicRod.hh"
#include "3rdparty/ElasticKnots/3rdparty/ElasticRods/RodMaterial.hh"

// ElasticKnots
#include "3rdparty/ElasticKnots/PeriodicRodList.hh"
#include "3rdparty/ElasticKnots/ContactProblem.hh"

struct HessianAndGradient {
    Eigen::SparseMatrix<double, 0, int> H;
    Eigen::VectorXd g;
};

PeriodicRod define_periodic_rod(std::vector<Eigen::Vector3d> pts, RodMaterial material);
std::vector<Eigen::Vector3d> read_nodes_from_file(std::string &filename);
std::vector<Eigen::Vector3d> DoFsToPos(Eigen::VectorXd dofs, uint n_pts);
Eigen::VectorXd DoFsToTwist(Eigen::VectorXd dofs, uint n_pts);
Eigen::MatrixXd toEigenDense(SuiteSparseMatrix& ssm);
Eigen::SparseMatrix<double>  toEigenSparse(SuiteSparseMatrix& ssm);
Eigen::SparseMatrix<double> computeHessian(ContactProblem& cp);
std::vector<Eigen::Vector3d> reduce_knot_resolution(std::vector<Eigen::Vector3d> pts, size_t factor);
HessianAndGradient removeTwist(Eigen::SparseMatrix<double, 0, int> H_sparse, Eigen::VectorXd g);
HessianAndGradient removeTheta(Eigen::SparseMatrix<double, 0, int> H_sparse, Eigen::VectorXd g);