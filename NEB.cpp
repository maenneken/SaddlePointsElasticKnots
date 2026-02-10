#include "NEB.h"

Eigen::VectorXd tangent(Eigen::VectorXd R_pre, Eigen::VectorXd R_next){
    Eigen::VectorXd t = R_next - R_pre;
    t.normalize();
    return t;
}
//F_spring = k * (||R_{i+1} - R_i|| - ||R_i - R_{i-1}||) * tangent_i
Eigen::VectorXd springForce(double spring_constant,Eigen::VectorXd R_pre, Eigen::VectorXd R, Eigen::VectorXd R_next, Eigen::VectorXd t){
    Eigen::VectorXd F = spring_constant*((R_next - R).norm() - (R - R_pre).norm()) * t;
    return F;
}
//F_perp =-d_R_{i} + d_R_{i}°tangent_i*tanget_i
Eigen::VectorXd perpForce(Eigen::VectorXd d_R, Eigen::VectorXd t){
    Eigen::VectorXd F = - d_R + (d_R.dot(t) * t);
    return F;
}
Eigen::VectorXd nebForce(double spring_constant, Eigen::VectorXd d_R, Eigen::VectorXd R_pre, Eigen::VectorXd R, Eigen::VectorXd R_next){
    Eigen::VectorXd t = tangent(R_pre, R_next);
    Eigen::VectorXd F = springForce(spring_constant,R_pre,R,R_next,t) + perpForce(d_R, t);
    return F;
}

//gradient Step
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path, double spring_constant, double step_size){
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next   = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient();
        Eigen::VectorXd F_neb = nebForce(spring_constant,d_R,R_pre,R,R_next);
        path[i] += step_size * F_neb;
    }
}

//TODO Global L-BFGS look at paper
//R_j+1 = R_j + F_j*H_j^−1
//R_j is the whole path, F_j are all forces H_j are all Hessians of the path
