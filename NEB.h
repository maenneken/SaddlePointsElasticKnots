#pragma once
#include "helpers.h"

Eigen::VectorXd tangent(Eigen::VectorXd R_pre, Eigen::VectorXd R_next);
Eigen::VectorXd springForce(double spring_constant,Eigen::VectorXd R_pre, Eigen::VectorXd R, Eigen::VectorXd R_next, Eigen::VectorXd t);
Eigen::VectorXd perpForce(Eigen::VectorXd d_R, Eigen::VectorXd t);
Eigen::VectorXd nebForce(double spring_constant, Eigen::VectorXd d_R, Eigen::VectorXd R_pre, Eigen::VectorXd R, Eigen::VectorXd R_next);
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path, double spring_constant, double step_size);