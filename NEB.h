#pragma once
#include "helpers.h"

Eigen::MatrixXd listToMatrix(std::vector<Eigen::VectorXd> vecs);
Eigen::VectorXd tangent(Eigen::VectorXd R_pre, Eigen::VectorXd R_next);
Eigen::VectorXd springForce(double spring_constant,Eigen::VectorXd R_pre, Eigen::VectorXd R, Eigen::VectorXd R_next, Eigen::VectorXd t);
Eigen::VectorXd perpForce(Eigen::VectorXd d_R, Eigen::VectorXd t);
Eigen::VectorXd nebForce(double spring_constant, Eigen::VectorXd d_R, Eigen::VectorXd R_pre, Eigen::VectorXd R, Eigen::VectorXd R_next);
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path, double spring_constant, double step_size);
void fitt_globalNebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path, double spring_constant, double step_size);
void globalNebGradientStep(Eigen::MatrixXd& gradient, Eigen::MatrixXd& path,double spring_constant,double step_size);