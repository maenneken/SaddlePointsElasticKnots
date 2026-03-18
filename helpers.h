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

#include "KnotVisualizer.h"


#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RESET   "\033[0m"

struct HessianAndGradient {
    Eigen::MatrixXd H;
    Eigen::VectorXd g;

    HessianAndGradient() = default;

    HessianAndGradient(Eigen::MatrixXd& H_, Eigen::VectorXd& g_){
    H = H_;
    g = g_;
}
};

PeriodicRod define_periodic_rod(std::vector<Eigen::Vector3d> pts, RodMaterial material);
std::vector<Eigen::Vector3d> read_nodes_from_file(const std::string &filename);
std::vector<Eigen::Vector3d> DoFsToPos(Eigen::VectorXd dofs, uint n_pts);
Eigen::VectorXd DoFsToTwist(Eigen::VectorXd dofs, uint n_pts);
Eigen::MatrixXd toEigenDense(SuiteSparseMatrix& ssm);
Eigen::SparseMatrix<double>  toEigenSparse(SuiteSparseMatrix& ssm);
Eigen::SparseMatrix<double> computeHessian(ContactProblem& cp);
std::vector<Eigen::Vector3d> reduce_knot_resolution(std::vector<Eigen::Vector3d> pts, size_t factor);
HessianAndGradient makeSmaller(Eigen::MatrixXd &H , Eigen::VectorXd& g, size_t k);
HessianAndGradient insertInBiggerHg(Eigen::MatrixXd & H_big, Eigen::VectorXd& g_big, Eigen::MatrixXd & H_small, Eigen::VectorXd& g_small);
void savePathTxt(const std::string& filename, const std::vector<Eigen::VectorXd>& path);
std::vector<Eigen::VectorXd>loadPathTxt(const std::string& filename);
void showPath(std::vector<Eigen::VectorXd>& path, ContactProblem& cp, KnotVisualizer& Viewer);
void relax_start_goal(ContactProblem& cp, Eigen::VectorXd& start_dofs,  Eigen::VectorXd& goal_dofs);