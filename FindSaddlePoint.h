#include "helpers.h"
#include <iostream>
#include <vector>


#include "KnotVisualizer.h"

Eigen::VectorXd reflectGradient(Eigen::VectorXd g, Eigen::SparseMatrix<double, 0, int> H_sparse, int saddleType, double etol, bool flipAllNegEig);
Eigen::VectorXd stepTowardsSaddle(ContactProblem& cp, double step, int saddleType, bool useTwist, double etol, bool clampedEnds, bool flipAllNegEig);
Eigen::VectorXd stepTowardsSaddleNewton(ContactProblem& cp, double step, int saddleType, bool useTwist, double etol, bool clampedEnds, bool flipAllNegEig);