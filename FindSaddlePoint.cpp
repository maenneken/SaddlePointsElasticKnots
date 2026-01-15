#include <iostream>
#include <vector>

#include "helpers.h"
#include "KnotVisualizer.h"
#include "FindSaddlePoint.h"
#include "distanceEnergy.h"


double eigReflFactor = -1.0;
int distEnergy_minSeperation = 3;
double distEnergy_weight = 0.1;

Eigen::VectorXd reflectGradient(Eigen::VectorXd g, Eigen::SparseMatrix<double, 0, int> H_sparse, int saddleType, double etol, bool flipAllNegEig){

    //nothing gets reflected
    if (saddleType <=0) return g;
    
    Eigen::MatrixXd H_dense = Eigen::MatrixXd(H_sparse);


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_dense);
    auto u = solver.eigenvalues();
    auto U = solver.eigenvectors();

    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigen decomposition failed\n";
    }

    Eigen::VectorXd g_p = U.transpose() * g;

    //Eigenvalues are sorted so 0 is always the smallest 
    int e = 0;
    if(flipAllNegEig){
        while(u(e) < 0){ //it is a negative eigenvalue
            g_p(e) *= eigReflFactor;
            ++e;
            //std::cout << u(e) <<" (neg),  ";
        }
        std::cout << "|neg eigenvalues| = " << e << std::endl;
    }
    //only goes in this loop if flipAllNegEig == false or we did not flip enough eig for the desired saddletype
    for (int i = 0; i < g_p.size() && e < saddleType; ++i){   
        if( std::abs(u(i)) > etol){//ignore the ones close to zero
            g_p(i) *= eigReflFactor;
            ++e;
            std::cout << u(i) <<", ";
        }
            
    }

    Eigen::VectorXd g_new = U * g_p;

    return g_new;
}

Eigen::VectorXd stepTowardsSaddle(ContactProblem& cp, double step, int saddleType, bool useTwist, double etol, bool clampedEnds, bool flipAllNegEig){

    auto R = cp.getVars();
    auto g = cp.gradient() + distanceGradient(R);
    auto H_sparse = computeHessian(cp) + distanceHessian(R);

    Eigen::VectorXd g_new;

    if(!useTwist){
        HessianAndGradient Hg = removeTwist(H_sparse,g);
        Eigen::VectorXd g_new_small = reflectGradient(Hg.g, Hg.H,saddleType, etol, flipAllNegEig);
        //std::cout << g.size() << " " << Hg.g.size() << std::endl;
        g_new = g;
        for (size_t i = 0; i < g_new_small.size();++i){
            g_new(i) = g_new_small(i);
        }

    } else if(!clampedEnds) { //when keeping theta constant the knot is free to rotate (look supplemental.pdf)
        //remove it for correct reflection calk
        HessianAndGradient Hg = removeTheta(H_sparse,g);
        Eigen::VectorXd g_new_small = reflectGradient(Hg.g, Hg.H,saddleType,etol,flipAllNegEig);
        g_new = g;
        for (size_t i = 0; i < g_new_small.size();++i){
            g_new(i) = g_new_small(i);
        }
        g_new(g_new.size()-1) = 0; //keep theta constant

    } else {
        g_new = reflectGradient(g, H_sparse,saddleType,etol, flipAllNegEig);
    }

    R = R - step * g_new;
    
    cp.setVars(R);
    
    return g_new;

}

Eigen::VectorXd stepTowardsSaddleNewton(ContactProblem& cp, double step, int saddleType, bool useTwist, double etol, bool clampedEnds, bool flipAllNegEig){

    auto R = cp.getVars();
    auto g = cp.gradient();
    auto H_sparse = computeHessian(cp);

    auto H_dist = distanceHessian(R,distEnergy_minSeperation);
    auto g_dist = distanceGradient(R,distEnergy_minSeperation);

    g+= distEnergy_weight * g_dist;
    H_sparse += distEnergy_weight * H_dist;

    /*
    SuiteSparseMatrix H = cp.hessian();
    auto H_sparse = toEigenSparse(H);
    */
    Eigen::VectorXd g_new;

    if(!useTwist){   
        HessianAndGradient Hg = removeTwist(H_sparse,g);
        Eigen::VectorXd g_new_small = reflectGradient(Hg.g, Hg.H,saddleType,etol,flipAllNegEig);
        g_new = g;
        for (size_t i = 0; i < g_new_small.size();++i){
            g_new(i) = g_new_small(i);
        }

    } else if(!clampedEnds) { //when keeping theta constant the knot is free to rotate (look supplemental.pdf)
        //remove it for correct reflection calk
        HessianAndGradient Hg = removeTheta(H_sparse,g);
        Eigen::VectorXd g_new_small = reflectGradient(Hg.g, Hg.H,saddleType,etol,flipAllNegEig);
        g_new = g;
        for (size_t i = 0; i < g_new_small.size();++i){
            g_new(i) = g_new_small(i);
        }
        g_new(g_new.size()-1) = 0; //keep theta constant

    }
    else {
        g_new = reflectGradient(g, H_sparse,saddleType, etol,flipAllNegEig);
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> newtonSolver;
    newtonSolver.compute(H_sparse);

    if (newtonSolver.info() != Eigen::Success) {
        throw std::runtime_error("Factorization failed");
    }

    Eigen::VectorXd newt_step = newtonSolver.solve(g_new);

    if (newtonSolver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
    }
    R = R - step * newt_step;

    cp.setVars(R);
    
    return g_new;

}
int main(int argc, char** argv) {
    std::string file = "../data/NoCollision/reduced0033.obj";
    double rod_radius = 0.2;
    int reductionFactor = 4;
    bool hasCollisions = true;
    int contactStiffness = 10000;

    // Parse command line arguments
    if (argc >= 2) {
        file = argv[1];
    }
    if (argc >= 3) {
        reductionFactor = std::stod(argv[2]);
    }
    if (reductionFactor < 1){
        reductionFactor = 1;
    }


    std::vector<double> params = {rod_radius, rod_radius};

    //Create material
    RodMaterial material(
        "ellipse",
        2000,     // Young's modulus
        0.3,      // Poisson's ratio
        params
    );

    // Read centerline nodes
    std::vector<Eigen::Vector3d> centerline = reduce_knot_resolution(read_nodes_from_file(file), reductionFactor);

    int n_pts = centerline.size();  

    PeriodicRod pr = define_periodic_rod(centerline,material);

    std::cout << "created pr"<< std::endl;
    
    // 4. Wrap into PeriodicRodList
    PeriodicRodList rod_list = PeriodicRodList(pr);

    std::cout << "created rodlist"<< std::endl;
    // 5. Setup problem options
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = hasCollisions;
    problemOptions.contactStiffness = contactStiffness;
    problemOptions.dHat = 2 * rod_radius;

    std::cout << "set pr options"<< std::endl;
    
    // 6. Create contact problem
    ContactProblem cp(rod_list, problemOptions);

    //Set Visulizer 
    KnotVisualizer Viewer = KnotVisualizer();
    Viewer.setKnot(centerline,0.01 * rod_radius);
    std::vector<Eigen::Vector3d> theta_point = {Eigen::Vector3d(0,0,0), Eigen::Vector3d(0,1000,0)};
    Viewer.setTheta(theta_point, 0.01 * rod_radius);
    //Save first Knot Vars
    auto startKnot = cp.getVars();

    static int iterations = 10000;
    static double stepsize = 0.001;
    static double eigen_tol = 0; //Tolerance to when an eigenvalue is treated like zero
    static int saddletype = 1;
    static bool running = false;
    static bool useNewton = true;
    static bool useTwist = true;
    static bool clampedEnds = true;
    static bool flipAllNegEig = false;
    static size_t i = 0;
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add new Energy based on different distance fkt to add information. ex: 1/(d²+1) or something with e^x
        //todo add option for distance fkt where you can decide how many neighbors are included in the calkulation (direct neighbors might not be important to know where they are)
        //todo add gradient and Hessian of the energy fkt
        //Idea is to find negative Eigenvalues
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("stepsize", &stepsize,(0.0001),(0.001),"%.7f");
        ImGui::InputInt("Saddle Type", &saddletype);
        ImGui::InputDouble("Eigen Zero Tolerance", &eigen_tol, (0.01),(0.1),"%.3f"); 
        ImGui::InputDouble("Eigenvalue Reflection Factor", &eigReflFactor); 
        ImGui::InputDouble("Distance Energy weight", &distEnergy_weight, (0.01),(0.1),"%.3f"); 
        ImGui::InputInt("Distance Energy minseperation", &distEnergy_minSeperation);
        ImGui::Checkbox("Use Newton", &useNewton);
        ImGui::Checkbox("Use Twist", &useTwist);
        ImGui::Checkbox("Clamped Ends", &clampedEnds);
        ImGui::Checkbox("Reflect all negative Eigenvalues", &flipAllNegEig);
        if (ImGui::Button("Find Saddle")) {
            running = true;
            i = 0;                           
        }
        if (ImGui::Button("Reset")){
            cp.setVars(startKnot);
            auto pts = DoFsToPos(cp.getVars(),n_pts);
            Viewer.updateKnot(pts);
        }
        ImGui::End();   

    });
    
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        if (running){
            Eigen::VectorXd g_new;
            Eigen::VectorXd g_old = cp.gradient();
            if(useNewton){
                g_new = stepTowardsSaddleNewton(cp,stepsize,saddletype,useTwist,eigen_tol,clampedEnds, flipAllNegEig);
            } else {
                g_new = stepTowardsSaddle(cp,stepsize,saddletype,useTwist,eigen_tol,clampedEnds, flipAllNegEig);
            }
            
            if(i % 100 == 0){
                auto pts = DoFsToPos(cp.getVars(),n_pts);
                auto twist = DoFsToTwist(cp.getVars(),n_pts);
                Viewer.updateKnot(pts);
                Viewer.colorTwist(twist);
                Viewer.showNodeGradient(DoFsToPos(g_old ,n_pts));
                Viewer.showNodeGradientModified(DoFsToPos(g_new ,n_pts));
                
                auto last = cp.getVars().size()-1;
                Viewer.updateTheta(cp.getVars()[last], g_old[last], g_new[last]);
                //Viewer.showContactForce(DoFsToPos(cp.contactForces(),n_pts));
                cp.contactEnergy();
                Viewer.frameTick();
                auto Dofs = cp.getVars();
                std::cout       << std::endl << std::endl
                                << "    It: " << i
                                << ", current Energy: " << cp.energy() 
                                << ", Gradient: " << g_old.norm() 
                                << ", reflected Gradient: " << g_new.norm()
                                << ", difference: " << (g_old- g_new).norm()
                                << ", Position: " << cp.getVars().norm()
                                << ", Twist Min/Max: " << twist.minCoeff() << " / " <<twist.maxCoeff()
                                << ", dist Energy: " << distEnergy_weight * distanceEnergy(Dofs,distEnergy_minSeperation)
                                << ", dist gradient: " << distEnergy_weight * distanceGradient(Dofs,distEnergy_minSeperation).norm()
                                << std::endl << std::endl;

            }
            /** 
            if(i % 1000 == 0){
                std::cout << g_old - g_new << std::endl << std::endl;
            }
            */
            i++;
        }
        if(i > iterations){
            i = 0;
            running = false;
        }
        Viewer.frameTick();
    }

    return 0;
}
