#include "helpers.h"
#include <iostream>
#include <vector>


#include "KnotVisualizer.h"




Eigen::VectorXd smashKnotStep(ContactProblem& cp, double step, double smashPlane, int axis){

    auto R = cp.getVars();
    auto g = cp.gradient();

    Eigen::VectorXd g_new;

    g_new = g;

    double k=1;
    
    for (int i = axis; i < R.size() * 0.75; i+=3){
        g_new(i) += std::exp(R(i) - smashPlane);
        g_new(i) -= std::exp(-R(i) - smashPlane);

    }

    R = R - step * g_new;
    cp.setVars(R);
    
    return g_new;

}

Eigen::VectorXd smashKnotStepNewton(ContactProblem& cp, double step, double smashPlane, int axis){

    auto R = cp.getVars();
    auto g = cp.gradient();
    auto H_sparse = computeHessian(cp);

    Eigen::VectorXd g_new;

    g_new = g;


    for (int i = axis; i < R.size() * 0.75; i+=3){
        g_new(i) += std::exp(R(i) - smashPlane);
        g_new(i) -= std::exp(-R(i) - smashPlane);

        H_sparse.coeffRef(i, i) +=
            std::exp(R(i) - smashPlane)
            + std::exp(-R(i) - smashPlane);
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
    std::string file = "../data/NoCollision/reduced0001.obj";
    double rod_radius = 0.2;
    int reductionFactor = 1;
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

    static int iterations = 100000;
    static double stepsize = 1;
    static double smashSpeed = 0.001;
    static bool running = false;
    static bool useNewton = true;
    static int smashAxis = 1;
    double smashPlane = 30;

    static size_t i = 0;
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("stepsize", &stepsize,(0.0001),(0.001),"%.7f");
        ImGui::InputDouble("Smash speed", &smashSpeed,(0.0001),(0.001),"%.7f");
        ImGui::InputInt("smashing axis", &smashAxis);
        ImGui::InputDouble("Smashing Plane", &smashPlane);
        ImGui::Checkbox("Use Newton", &useNewton);

        if (ImGui::Button("Smash")) {
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
                g_new = smashKnotStepNewton(cp,stepsize, smashPlane, smashAxis);
            } else {
                g_new = smashKnotStep(cp,stepsize, smashPlane, smashAxis);
            }

            smashPlane -= smashSpeed;
            if(smashPlane <= 5){
                smashPlane = 5;
                smashSpeed = 0;
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
                Viewer.frameTick();
                std::cout       << "    It: " << i
                                << ", current Energy: " << cp.energy() 
                                << ", Gradient: " << g_old.norm() 
                                << ", reflected Gradient: " << g_new.norm()
                                << ", difference: " << (g_old- g_new).norm()
                                << ", Position: " << cp.getVars().norm()
                                << ", Twist Min/Max: " << twist.minCoeff() << " / " <<twist.maxCoeff()
                                << ", ContactEnergy" << cp.contactEnergy()
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
