#pragma once
#include "helpers.h"


Eigen::MatrixXd listToMatrix(const std::vector<Eigen::VectorXd>& vecs);
double lineStepNewton(ContactProblem& cp, const Eigen::VectorXd& R,const Eigen::VectorXd& d,const Eigen::VectorXd& F, double eps);
Eigen::VectorXd tangent(const Eigen::VectorXd& R_pre, const  Eigen::VectorXd & R_next);
Eigen::VectorXd springForce(double spring_constant,const Eigen::VectorXd& R_pre, const Eigen::VectorXd& R,const Eigen::VectorXd& R_next,const Eigen::VectorXd& t);
Eigen::VectorXd perpForce(const Eigen::VectorXd& d_R, const Eigen::VectorXd& t);
Eigen::VectorXd nebForce(double spring_constant, const Eigen::VectorXd & d_R,const Eigen::VectorXd& R_pre,const Eigen::VectorXd& R,const Eigen::VectorXd& R_next);
Eigen::VectorXd climbingForce(const Eigen::VectorXd& d_R,const  Eigen::VectorXd& R_pre, const Eigen::VectorXd& R_next);
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F_neb, double spring_constant, double step_size);
void nebConjugateGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F,std::vector<Eigen::VectorXd>& d, double spring_constant, double step_size);
void ciNebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F_neb, double spring_constant, double step_size, size_t max_id);
void globalNebGradientStep(Eigen::MatrixXd& gradient, Eigen::MatrixXd& path, double spring_constant, double step_size);