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
RowMatrix gradientMatrix(ContactProblem& cp, const RowMatrix& R){
    RowMatrix G(R.rows(), R.cols());
    for (size_t i = 0; i < R.rows(); ++i) {
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
//HENKELMANN: This is the original, simple tangent definition. It can lead to kinks in the path, especially near the saddle point. An improved version is implemented in the function "improved_tangent" below, which takes into account the energy differences between neighboring images to produce a smoother tangent.
Eigen::VectorXd improved_tangent(const Eigen::VectorXd& R_pre, const Eigen::VectorXd& R_curr,
                         const Eigen::VectorXd& R_next, double E_pre, double E_curr, double E_next)
{
    Eigen::VectorXd t_plus  = R_next - R_curr;
    Eigen::VectorXd t_minus = R_curr - R_pre;

    Eigen::VectorXd t;
    if (E_next > E_curr && E_curr > E_pre) {
        t = t_plus;
    } else if (E_next < E_curr && E_curr < E_pre) {
        t = t_minus;
    } else {
        // Sattelpunkt → gewichtete Interpolation
        double dE_max = std::max(std::abs(E_next - E_curr), std::abs(E_pre - E_curr));
        double dE_min = std::min(std::abs(E_next - E_curr), std::abs(E_pre - E_curr));
        if (E_next > E_pre)
            t = t_plus * dE_max + t_minus * dE_min;
        else
            t = t_plus * dE_min + t_minus * dE_max;
    }
    t.normalize();
    return t;
}
double nebEnergy(ContactProblem& cp, const std::vector<Eigen::VectorXd>& path, double spring_constant){
    for(size_t i = 1; i < path.size()-1;++i){
        cp.setVars(path[i]);
        double E = cp.energy();
        double E_spring = 0.5 * spring_constant * ((path[i+1] - path[i]).norm() - (path[i] - path[i-1]).norm());
        E += E_spring;
        return E;
    }
    return 0;
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
//with improve tangent
Eigen::VectorXd nebForce_improved(double spring_constant, const Eigen::VectorXd& d_R,
                          const Eigen::VectorXd& R_pre, const Eigen::VectorXd& R,
                          const Eigen::VectorXd& R_next,
                          double E_pre, double E_curr, double E_next)
{
    Eigen::VectorXd t = improved_tangent(R_pre, R, R_next, E_pre, E_curr, E_next);
    return springForce(spring_constant, R_pre, R, R_next, t) + perpForce(d_R, t);
}
//this force only acts on the image with the highest Energy
Eigen::VectorXd climbingForce(const Eigen::VectorXd& d_R, const Eigen::VectorXd& R_pre, const Eigen::VectorXd& R_next){
    Eigen::VectorXd t = tangent(R_pre, R_next);
    Eigen::VectorXd F = -d_R + 2 * (d_R.dot(t)*t);
    return F;
}
//with improved tangent
Eigen::VectorXd climbingForce_improved(const Eigen::VectorXd& d_R, const Eigen::VectorXd& R_pre,
                               const Eigen::VectorXd& R_next,
                               double E_pre, double E_curr, double E_next)
{
    Eigen::VectorXd t = improved_tangent(R_pre, 
                            0.5*(R_pre+R_next), // R_curr approximation
                            R_next, E_pre, E_curr, E_next);
    return -d_R + 2.0 * (d_R.dot(t) * t);
}
RowMatrix globalNEBForce(RowMatrix& d_R, RowMatrix& R, 
                          double spring_constant, bool climbing_image, 
                          size_t idx_max)
{
    int N    = R.rows();
    int cols = R.cols();
    RowMatrix F = RowMatrix::Zero(N, cols);

    // ---- Endpoints: gradient only ----
    F.row(0)     = -d_R.row(0);
    F.row(N - 1) = -d_R.row(N - 1);

    // ---- Tangents (internal only) ----
    RowMatrix tangents = R.bottomRows(N - 2) - R.topRows(N - 2);
    for (int i = 0; i < tangents.rows(); ++i)
        tangents.row(i).normalize();

    // ---- Spring force ----
    RowMatrix d_forward  = R.middleRows(2, N - 2) - R.middleRows(1, N - 2);
    RowMatrix d_backward = R.middleRows(1, N - 2) - R.middleRows(0, N - 2);
    Eigen::VectorXd len_f = d_forward.rowwise().norm();
    Eigen::VectorXd len_b = d_backward.rowwise().norm();
    Eigen::VectorXd delta = spring_constant * (len_f - len_b);
    RowMatrix F_spring    = delta.asDiagonal() * tangents;

    // ---- Perpendicular force ----
    RowMatrix grad_mid    = d_R.middleRows(1, N - 2);
    Eigen::VectorXd proj  = (grad_mid.cwiseProduct(tangents)).rowwise().sum();
    RowMatrix F_perp      = -grad_mid + proj.asDiagonal() * tangents;

    // ---- Assign internal forces ----
    F.middleRows(1, N - 2) = F_spring + F_perp;

    // ---- Climbing image ----
    if (climbing_image) {
        F.row(idx_max) = climbingForce(
            d_R.row(idx_max), 
            R.row(idx_max - 1), 
            R.row(idx_max + 1));
    }

    return F;
}

//gradient Step
void nebGradientStep(ContactProblem& cp, std::vector<Eigen::VectorXd>& path,std::vector<Eigen::VectorXd>& F_neb, double spring_constant, double step_size, size_t max_id, bool climbingImage){
    for( size_t i = 1; i< path.size()-1 ; ++i){
        Eigen::VectorXd R_pre   = path[i-1];
        Eigen::VectorXd R       = path[i];
        Eigen::VectorXd R_next  = path[i+1];
        cp.setVars(R);
        Eigen::VectorXd d_R = cp.gradient(true);
        double E_curr = cp.energy();
        cp.setVars(R_pre);
        double E_pre = cp.energy();
        cp.setVars(R_next);
        double E_next = cp.energy();
        //we are at image with highest Energy -> go towards Saddle
        if(climbingImage && i == max_id){
            //F_neb[i] = climbingForce(d_R,R_pre,R_next);
            F_neb[i] = climbingForce_improved(d_R,R_pre,R_next,E_pre,E_curr,E_next);

        }
        //else do normal NEB
        else{
            //F_neb[i] = nebForce(spring_constant,d_R,R_pre,R,R_next);
            F_neb[i] = nebForce_improved(spring_constant,d_R,R_pre,R,R_next,E_pre,E_curr,E_next);
        }
        path[i] += step_size * F_neb[i];
    }
}

void globalNebLineSearchStep(ContactProblem& cp, RowMatrix& R,
                              double spring_constant, bool climbing_image,
                              size_t idx_max, double& step_size)
{
    int n    = R.rows();
    int cols = R.cols();

    RowMatrix d_R   = gradientMatrix(cp, R);
    RowMatrix F_mat = globalNEBForce(d_R, R, spring_constant,
                                     climbing_image, idx_max);

    Eigen::VectorXd F = Eigen::Map<Eigen::VectorXd>(
        F_mat.data(), F_mat.size());
    const double F_norm_sq = F.squaredNorm();

    const double c      = 1e-2;
    const double decay  = 0.5;
    const int    max_ls = 10;

    double best_size    = step_size;
    double best_norm_sq = std::numeric_limits<double>::infinity();

    // Save full path including endpoints
    Eigen::VectorXd R_old = Eigen::Map<Eigen::VectorXd>(
        R.data(), n * cols);

    double alpha = step_size;
    for (int ls = 0; ls < max_ls; ++ls) {
        RowMatrix R_trial = R;
        Eigen::Map<Eigen::VectorXd>(
            R_trial.data(), n * cols) = R_old + alpha * F;

        RowMatrix d_R_trial   = gradientMatrix(cp, R_trial);
        RowMatrix F_trial_mat = globalNEBForce(
            d_R_trial, R_trial, spring_constant, climbing_image, idx_max);

        double F_trial_norm_sq = F_trial_mat.squaredNorm();

        if (F_trial_norm_sq < best_norm_sq) {
            best_norm_sq = F_trial_norm_sq;
            best_size    = alpha;
        }
        if (F_trial_norm_sq <= F_norm_sq * (1.0 - 2.0 * c * alpha))
            break;

        alpha *= decay;
    }

    // Apply best step to full path
    Eigen::Map<Eigen::VectorXd>(R.data(), n * cols)
        = R_old + best_size * F;

    // Adaptive step size
    if (best_size >= step_size * 0.8)
        step_size = std::min(step_size * 1.2, 0.5);
    else if (best_size <= step_size * 0.2)
        step_size = std::max(step_size * 0.7, 1e-4);
}

void globalNebLBFGSHesStep(ContactProblem& cp, RowMatrix& R, LBGFhistory& history,
                            double spring_constant, bool climbing_image,
                            size_t idx_max, double max_step)
{
    int rows_internal = R.rows() - 2;
    int cols          = R.cols();

    // ── Kraft am aktuellen Punkt ──────────────────────────────────────────────
    RowMatrix        d_R   = gradientMatrix(cp, R);
    RowMatrix        F_mat = globalNEBForce(d_R, R, spring_constant, climbing_image, idx_max);
    Eigen::VectorXd  F     = Eigen::Map<Eigen::VectorXd>(F_mat.data(), F_mat.size());

    // ── L-BFGS two-loop recursion ─────────────────────────────────────────────
    int k = static_cast<int>(history.s_list.size());
    std::vector<double> alpha(k);
    Eigen::VectorXd q = F;

    for (int i = k - 1; i >= 0; --i) {
        alpha[i] = history.rho_list[i] * history.s_list[i].dot(q);
        q -= alpha[i] * history.y_list[i];
    }

    Eigen::VectorXd z = q;
    if (k > 0) {
        double sy = history.y_list[k-1].dot(history.s_list[k-1]);
        double yy = history.y_list[k-1].dot(history.y_list[k-1]);
        double gamma = (yy > 1e-14) ? std::clamp(sy / yy, 0.01, 10.0) : 1.0;
        z *= gamma;
    }

    for (int i = 0; i < k; ++i) {
        double beta = history.rho_list[i] * history.y_list[i].dot(z);
        z += history.s_list[i] * (alpha[i] - beta);
    }

    // ── Schritt begrenzen ─────────────────────────────────────────────────────
    Eigen::VectorXd step = z;
    double norm = step.norm();
    if (norm > max_step)
        step *= max_step / norm;

    // ── R aktualisieren ────────────────────────────────────────────────────────
    Eigen::VectorXd F_old = F;
    RowMatrix R_internal_old = R.middleRows(1, rows_internal);

    R.middleRows(1, rows_internal) +=
        Eigen::Map<const RowMatrix>(step.data(), rows_internal, cols);

    // ── Neue Kraft ─────────────────────────────────────────────────────────────
    d_R   = gradientMatrix(cp, R);
    F_mat = globalNEBForce(d_R, R, spring_constant, climbing_image, idx_max);
    F     = Eigen::Map<Eigen::VectorXd>(F_mat.data(), F_mat.size());

    // ── History-Update mit Powell Damping ────────────────────────────────────
    Eigen::VectorXd s = step;
    Eigen::VectorXd y = F - F_old;
    double ys = y.dot(s);
    double ss = s.dot(s);

    const double damping_threshold = 0.2;
    double threshold = damping_threshold * ss;

    if (ys < threshold) {
        double denom = ss - ys;
        if (std::abs(denom) > 1e-14) {
            double theta = (ss - threshold) / denom;
            theta = std::clamp(theta, 0.0, 1.0);
            y = theta * y + (1.0 - theta) * s;
            ys = y.dot(s);
        }
    }

    std::cout << "History size: " << history.s_list.size()
              << ", ys: " << ys << std::endl;

    if (ys <= 1e-10)
        return;

    if (history.s_list.size() == history.n_memory) {
        history.s_list.erase(history.s_list.begin());
        history.y_list.erase(history.y_list.begin());
        history.rho_list.erase(history.rho_list.begin());
    }

    history.s_list.push_back(s);
    history.y_list.push_back(y);
    history.rho_list.push_back(1.0 / ys);
}