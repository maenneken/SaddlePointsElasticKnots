#include "NEB.h"

Eigen::MatrixXd listToMatrix(const std::vector<Eigen::VectorXd>& vecs){
    size_t N = vecs.size();
    size_t d = vecs[0].size();

    Eigen::MatrixXd M(N, d);
    for (size_t i = 0; i < N; ++i) {
        M.row(i) = vecs[i].transpose();
    }
    return M;
}
double lineStepNewton(ContactProblem& cp, const Eigen::VectorXd& R,const Eigen::VectorXd& d,const Eigen::VectorXd& F, double eps)
{
    // first derivative along direction
    double f1 = F.dot(d);

    // finite difference force evaluation
    Eigen::VectorXd R_eps = R + eps * d;
    cp.setVars(R_eps);
    Eigen::VectorXd F_eps = -cp.gradient(true);

    double f2 = (F_eps.dot(d) - f1) / eps;

    if(std::abs(f2) < 1e-14)
        return 0.0;

    // single Newton step
    return -f1 / f2;
}

Eigen::VectorXd tangent(const Eigen::VectorXd& R_pre, const  Eigen::VectorXd & R_next){
    Eigen::VectorXd t = R_next - R_pre;
    t.normalize();
    return t;
}
//F_spring = k * (||R_{i+1} - R_i|| - ||R_i - R_{i-1}||) * tangent_i
Eigen::VectorXd springForce(double spring_constant,const Eigen::VectorXd& R_pre, const Eigen::VectorXd& R,const Eigen::VectorXd& R_next,const Eigen::VectorXd& t){
    Eigen::VectorXd F = spring_constant*((R_next - R).norm() - (R - R_pre).norm()) * t;
    return F;
}
//F_perp =-d_R_{i} + d_R_{i}°tangent_i*tanget_i
Eigen::VectorXd perpForce(const Eigen::VectorXd& d_R, const Eigen::VectorXd& t){
    Eigen::VectorXd F = - d_R + (d_R.dot(t) * t);
    return F;
}
Eigen::VectorXd nebForce(double spring_constant, const Eigen::VectorXd & d_R,const Eigen::VectorXd& R_pre,const Eigen::VectorXd& R,const Eigen::VectorXd& R_next){
    Eigen::VectorXd t = tangent(R_pre, R_next);
    Eigen::VectorXd F = springForce(spring_constant,R_pre,R,R_next,t) + perpForce(d_R, t);
    return F;
}
//this force only acts on the image with the highest Energy
Eigen::VectorXd climbingForce(const Eigen::VectorXd& d_R, const Eigen::VectorXd& R_pre, const Eigen::VectorXd& R_next){
    Eigen::VectorXd t = tangent(R_pre, R_next);
    Eigen::VectorXd F = -d_R + 2 * (d_R.dot(t)*t);
    return F;
}

//gradient Step
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F_neb, double spring_constant, double step_size){
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next  = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient(true);
        F_neb[i] = nebForce(spring_constant,d_R,R_pre,R,R_next);
        path[i] += step_size * F_neb[i];
    }
}

//gradient Step Nudged elastic Band with climbing image
void ciNebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F_neb, double spring_constant, double step_size, size_t max_id){
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next   = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient(true);


        //we are at image with highest Energy -> go towards Saddle
        if(i == max_id){
            F_neb[i] = climbingForce(d_R,R_pre,R_next);

        }
        //else do normal NEB
        else{
            F_neb[i] = nebForce(spring_constant,d_R,R_pre,R,R_next);
        }
        path[i] += step_size * F_neb[i];
        
    }
}
//Conjugate gradient step, to improve performance
//directly out of Optimization methods for finding minimum energy paths
void nebConjugateGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F,std::vector<Eigen::VectorXd>& d, double spring_constant, double step_size){
    const double fd_eps = 1e-6;
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next  = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient(true);
        // CURRENT force
        Eigen::VectorXd F_old = nebForce(spring_constant,d_R,R_pre,R,R_next);

        // initialize direction
        if(d[i].norm() < 1e-12)
            d[i] = F_old;
        
        //calculate step size
        double lambda = lineStepNewton(cp,R,d[i],F[i],fd_eps);
        lambda = std::clamp(lambda,-1.0,1.0);

        //calculate step
        path[i] += lambda * d[i];
        
        //apply to contact Problem
        cp.setVars(path[i]);
        d_R = cp.gradient(true);

        //new Force
        Eigen::VectorXd F_new = nebForce(spring_constant,d_R, path[i-1],path[i],path[i+1]);

        // Polak–Ribiere
        double denom = std::max(F_old.squaredNorm(),1e-14);

        double gamma =
            F_new.dot(F_new - F_old) / denom;

        gamma = std::max(0.0, gamma);

        //Restart CG to stop stalling
        if(F_new.dot(d[i]) <= 0.0){
            //Restart
            std::cout << "restarting";
            d[i] = F_new;
        }
        else {
            // new direction
            d[i] = F_new + gamma * d[i];
        }
        

        // update stored force
        F[i] = F_new;
    }
}

void globalNebGradientStep(Eigen::MatrixXd& gradient, Eigen::MatrixXd& path, double spring_constant, double step_size){
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




void globalNebLBFGS(Eigen::MatrixXd& d_R, Eigen::MatrixXd& R, double spring_constant, bool climbing_image, size_t idx_max){

    int N = R.rows();

    // ---- Tangents ----
    Eigen::MatrixXd tangents = R.bottomRows(N-2) - R.topRows(N-2);
    for (int i = 0; i < tangents.rows(); ++i)
        tangents.row(i).normalize();

    // ---- Spring force ----
    Eigen::MatrixXd d_forward = R.middleRows(2, N-2) 
                              - R.middleRows(1, N-2);

    Eigen::MatrixXd d_backward = R.middleRows(1, N-2) 
                               - R.middleRows(0, N-2);

    Eigen::VectorXd len_f = d_forward.rowwise().norm();
    Eigen::VectorXd len_b = d_backward.rowwise().norm();

    Eigen::VectorXd delta = spring_constant * (len_f - len_b);
    Eigen::MatrixXd F_spring = delta.asDiagonal() * tangents;

    // ---- Perpendicular force ----
    Eigen::MatrixXd grad_mid = d_R.middleRows(1, N-2);
    Eigen::VectorXd proj = (grad_mid.cwiseProduct(tangents)).rowwise().sum();
    Eigen::MatrixXd F_perp = -grad_mid + proj.asDiagonal() * tangents;

    // ---- Total NEB force ----
    Eigen::MatrixXd F_neb = F_spring + F_perp;

    if(climbing_image){
        F_neb.row(idx_max-1) = climbingForce(d_R.row(idx_max),R.row(idx_max-1),R.row(idx_max+1));
    }


}
