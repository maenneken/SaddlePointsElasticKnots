#include "NEB.h"

Eigen::MatrixXd listToMatrix(std::vector<Eigen::VectorXd> vecs){
    size_t N = vecs.size();
    size_t d = vecs[0].size();

    Eigen::MatrixXd M(N, d);
    for (size_t i = 0; i < N; ++i) {
        M.row(i) = vecs[i].transpose();
    }
    return M;
}
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
//this force only acts on the image with the highest Energy
Eigen::VectorXd climbingForce(Eigen::VectorXd d_R, Eigen::VectorXd R_pre, Eigen::VectorXd R_next){
    Eigen::VectorXd t = tangent(R_pre, R_next);
    Eigen::VectorXd F = -d_R + 2 * (d_R.dot(t)*t);
    return F;
}

//gradient Step
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path, double spring_constant, double step_size){
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next   = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient(true);
        Eigen::VectorXd F_neb = nebForce(spring_constant,d_R,R_pre,R,R_next);
        path[i] += step_size * F_neb;
    }
}

//gradient Step Nudged elastic Band with climbing image
void ciNebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path, double spring_constant, double step_size, size_t max_id){
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next   = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient(true);


        //we are at image with highest Energy -> go towards Saddle
        if(i == max_id){
            Eigen::VectorXd F_ci = climbingForce(d_R,R_pre,R_next);
            path[i] += step_size * F_ci;
        }
        //else do normal NEB
        else{
            Eigen::VectorXd F_neb = nebForce(spring_constant,d_R,R_pre,R,R_next);
            path[i] += step_size * F_neb;
        }
        
    }
}
void globalNebGradientStep(
    Eigen::MatrixXd& gradient,  // N x D
    Eigen::MatrixXd& path,      // N x D
    double spring_constant,
    double step_size)
{
    int N = path.rows();

    // ---- Tangents ----
    Eigen::MatrixXd tangents = path.bottomRows(N-2) - path.topRows(N-2);
    for (int i = 0; i < tangents.rows(); ++i)
        tangents.row(i).normalize();

    // ---- Spring force ----
    Eigen::MatrixXd d_forward = path.middleRows(2, N-2) 
                              - path.middleRows(1, N-2);

    Eigen::MatrixXd d_backward = path.middleRows(1, N-2) 
                               - path.middleRows(0, N-2);

    Eigen::VectorXd len_f = d_forward.rowwise().norm();
    Eigen::VectorXd len_b = d_backward.rowwise().norm();

    Eigen::VectorXd delta = spring_constant * (len_f - len_b);
    Eigen::MatrixXd F_spring = delta.asDiagonal() * tangents;

    // ---- Perpendicular force ----
    Eigen::MatrixXd grad_mid = gradient.middleRows(1, N-2);
    Eigen::VectorXd proj = (grad_mid.cwiseProduct(tangents)).rowwise().sum();
    Eigen::MatrixXd F_perp = -grad_mid + proj.asDiagonal() * tangents;

    // ---- Total NEB force ----
    Eigen::MatrixXd F_neb = F_spring + F_perp;

    // ---- Update path ----
    path.middleRows(1, N-2) += step_size * F_neb;
}

//TODO Global L-BFGS look at paper
//R_j+1 = R_j + F_j*H_j^−1
//R_j is the whole path, F_j are all forces H_j are all Hessians of the path
