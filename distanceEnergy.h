#include "helpers.h"


double distanceEnergy(Eigen::VectorXd DoFs, int minSeparation = 3);
Eigen::VectorXd distanceGradient(Eigen::VectorXd DoFs, int minSeparation = 3);
Eigen::SparseMatrix<double> distanceHessian(Eigen::VectorXd DoFs, int minSeparation = 3);