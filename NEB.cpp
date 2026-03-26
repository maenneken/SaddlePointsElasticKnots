#include "NEB.h"

RowMatrix listToMatrix(const std::vector<Eigen::VectorXd>& vecs){
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
RowMatrix gradientMatrix(ContactProblem& cp, const RowMatrix& R){
    RowMatrix G(R.rows(), R.cols());
    //#pragma omp parallel for
    for (size_t i = 0; i < R.rows();++i){
            cp.setVars(R.row(i));
             G.row(i) = cp.gradient(true).transpose(); 
    }
    return G;
}
Eigen::VectorXd tangent(const Eigen::VectorXd& R_pre, const  Eigen::VectorXd & R_next){
    Eigen::VectorXd t = R_next - R_pre;
    t.normalize();
    return t;
}
//F_spring = k * (||R_{i+1} - R_i|| - ||R_i - R_{i-1}0.2||) * tangent_i
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
RowMatrix globalNEBForce(RowMatrix& d_R, RowMatrix& R, double spring_constant, bool climbing_image, size_t idx_max){
    int N = R.rows();

    // ---- Tangents ----
    RowMatrix tangents = R.bottomRows(N-2) - R.topRows(N-2);
    for (int i = 0; i < tangents.rows(); ++i)
        tangents.row(i).normalize();

    // ---- Spring force ----
    RowMatrix d_forward = R.middleRows(2, N-2) 
                              - R.middleRows(1, N-2);

    RowMatrix d_backward = R.middleRows(1, N-2) 
                               - R.middleRows(0, N-2);

    Eigen::VectorXd len_f = d_forward.rowwise().norm();
    Eigen::VectorXd len_b = d_backward.rowwise().norm();

    Eigen::VectorXd delta = spring_constant * (len_f - len_b);
    RowMatrix F_spring = delta.asDiagonal() * tangents;

    // ---- Perpendicular force ----
    RowMatrix grad_mid = d_R.middleRows(1, N-2);
    Eigen::VectorXd proj = (grad_mid.cwiseProduct(tangents)).rowwise().sum();
    RowMatrix F_perp = -grad_mid + proj.asDiagonal() * tangents;

    // ---- Total NEB force ----
    RowMatrix F_neb = F_spring + F_perp;

    if(climbing_image){
        F_neb.row(idx_max-1) = climbingForce(d_R.row(idx_max),R.row(idx_max-1),R.row(idx_max+1));
    }

    return F_neb;
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
        Eigen::VectorXd F_old = F[i];

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

        double gamma = F_new.dot(F_new - F_old) / denom;

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

//todo add line search to make it converge faster.
void globalNebLBFGSStep(ContactProblem& cp, RowMatrix& R, LBGFhistory& history, double spring_constant, bool climbing_image, size_t idx_max, double max_step){
    RowMatrix d_R = gradientMatrix(cp,R);
    RowMatrix F_mat = globalNEBForce(d_R,R,spring_constant,climbing_image,idx_max);
    
    //flatten 
    Eigen::VectorXd F = Eigen::Map<Eigen::VectorXd>(F_mat.data(), F_mat.size());
    int k = history.s_list.size();
    std::vector<double> alpha(k);
    Eigen::VectorXd q = F;


    //backwards
    for (int i = k-1; i >= 0; --i){
         alpha[i] = history.rho_list[i] * history.s_list[i].dot(q);
        q -= alpha[i] * history.y_list[i]; 
    }

    // Initial Hessian approximation
    Eigen::VectorXd z = q;
    if(k > 0){
        double last_sy = history.y_list[k-1].dot(history.s_list[k-1]);
        double last_yy = history.y_list[k-1].dot(history.y_list[k-1]);
        double gamma = last_sy / last_yy;
        z *= gamma;
    }

    // Second loop (forwards)
    for(int i = 0; i < k; ++i){
        double beta = history.rho_list[i] * history.y_list[i].dot(z);
        z += history.s_list[i] * (alpha[i] - beta);
    }

    // L-BFGS_Hess step (Eq. 8)
    //FH⁻1
    Eigen::VectorXd step = z;

    // Optional: step size control (highly recommended for NEB)
    if (step.norm() > max_step)
        step *= max_step / step.norm();

    // Update coordinates

    //R includes ends F does not
    int rows_internal = R.rows() - 2;
    int cols = R.cols();

    // Map internal part directly (NO COPY)
    Eigen::Map<Eigen::VectorXd> R_vec(R.middleRows(1, rows_internal).data(), rows_internal * cols);

    Eigen::VectorXd R_old = R_vec;

    //todo instead of claming step, calk best stpe size
    //R_vec += step;

    // Line search to find optimal step size
    double step_size = 1.0;
    const double c = 1e-4;  // Armijo condition constant
    const int max_ls_iter = 20;
    
    Eigen::VectorXd R_trial = R_vec + step_size * step;
    cp.setVars(Eigen::Map<const Eigen::MatrixXd>(R_trial.data(), rows_internal, cols));
    RowMatrix d_R_trial = gradientMatrix(cp, R);
    Eigen::VectorXd F_trial = Eigen::Map<Eigen::VectorXd>(
        globalNEBForce(d_R_trial, R, spring_constant, climbing_image, idx_max).data(),
        rows_internal * cols);
    
    double initial_reduction = -c * step_size * F.dot(step);
    
    for (int ls_iter = 0; ls_iter < max_ls_iter && F_trial.dot(step) > initial_reduction; ++ls_iter) {
        step_size *= 0.5;
        R_trial = R_vec + step_size * step;
        cp.setVars(Eigen::Map<const Eigen::MatrixXd>(R_trial.data(), rows_internal, cols));
        d_R_trial = gradientMatrix(cp, R);
        F_trial = Eigen::Map<Eigen::VectorXd>(
            globalNEBForce(d_R_trial, R, spring_constant, climbing_image, idx_max).data(),
            rows_internal * cols);
    }
    
    R_vec += step_size * step;
    cp.setVars(Eigen::Map<const Eigen::MatrixXd>(R_vec.data(), rows_internal, cols));

    //upodate history
    d_R = gradientMatrix(cp,R);
    Eigen::VectorXd F_old = F;
    F_mat = globalNEBForce(d_R,R,spring_constant,climbing_image,idx_max);
    F = Eigen::Map<Eigen::VectorXd>(F_mat.data(), F_mat.size());

    Eigen::VectorXd s = R_vec - R_old;
    Eigen::VectorXd y = F - F_old;

    double ys = y.dot(s);

    // Curvature condition
    if (ys <= 1e-10)
        return;

    double rho = 1.0 / ys;

    // limited memory: remove oldest entry
    if (history.s_list.size() == history.n_memory) {
        history.s_list.erase(history.s_list.begin());
        history.y_list.erase(history.y_list.begin());
        history.rho_list.erase(history.rho_list.begin());
    }

    history.s_list.push_back(s);
    history.y_list.push_back(y);
    history.rho_list.push_back(rho);
}
