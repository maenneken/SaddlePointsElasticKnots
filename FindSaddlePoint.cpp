#include "helpers.h"
#include <iostream>
#include <vector>


#include "KnotVisualizer.h"


Eigen::VectorXd stepTowardsSaddle(ContactProblem& cp, double step, int saddleType){
    auto g = cp.gradient();

    auto H_sparse = computeHessian(cp);
    //when there is conntact hessian throws error, because contact creates new entrys in contact hessian. When trying to add them together using  addWithSubSparsity it is not allowed
    auto R = cp.getVars();


    Eigen::MatrixXd H_dense = Eigen::MatrixXd(H_sparse);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_dense);
    auto u = solver.eigenvalues();
    auto U = solver.eigenvectors();

    Eigen::VectorXd g_p = U.transpose() * g;

    //Eigenvalues are sorted so 0 is always the smallest
    for (int i = 0; i < saddleType; ++i){
        g_p(i) *= -1.0;
    }
    

    Eigen::VectorXd g_new = U * g_p;

    R = R - step * g_new;

    cp.setVars(R);
    
    return u;

}

Eigen::VectorXd stepTowardsSaddleNewton(ContactProblem& cp, double step, int saddleType){
    auto g = cp.gradient();

    auto H_sparse = computeHessian(cp);
    //when there is conntact hessian throws error, because contact creates new entrys in contact hessian. When trying to add them together using  addWithSubSparsity it is not allowed
    auto R = cp.getVars();


    Eigen::MatrixXd H_dense = Eigen::MatrixXd(H_sparse);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_dense);
    auto u = solver.eigenvalues();
    auto U = solver.eigenvectors();

    Eigen::VectorXd g_p = U.transpose() * g;

    //Eigenvalues are sorted so 0 is always the smallest
    for (int i = 0; i < saddleType; ++i){
        g_p(i) *= -1.0;
    }
    

    Eigen::VectorXd g_new = U * g_p;

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
    
    return u;

}
int main(int argc, char** argv) {
    std::string file = "../data/NoCollision/reduced0033.obj";
    double rod_radius = 0.2;
    int reductionFactor = 4;
    bool hasCollisions = true;
    int contactStiffness = 1000;

    // Parse command line arguments
    if (argc >= 2) {
        file = argv[1];
    }
    if (argc >= 3) {
        reductionFactor = std::stod(argv[2]);
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

    //Save first Knot Vars
    auto startKnot = cp.getVars();

    static int iterations = 1000;
    static double stepsize = 0.001;
    static int saddletype = 1;
    static bool running = false;
    static bool useNewton = false;
    static size_t i = 0;
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options + reduceKnotResolution
        //todo show twist by changing the color of the edge
        //todo show gradient and modified gradient
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("stepsize", &stepsize,(0.0001),(0.001),"%.7f");
        ImGui::InputInt("Saddle Type", &saddletype);
        ImGui::Checkbox("Use Newton", &useNewton);
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
            Eigen::VectorXd eigenvalues;
            if(useNewton){
                eigenvalues = stepTowardsSaddleNewton(cp,stepsize,saddletype);
            } else {
                eigenvalues = stepTowardsSaddle(cp,stepsize,saddletype);
            }
            
            if(i % 100 == 0){
                auto pts = DoFsToPos(cp.getVars(),n_pts);
                auto twist = DoFsToTwist(cp.getVars(),n_pts);
                Viewer.updateKnot(pts);
                Viewer.colorTwist(twist);
                Viewer.showNodeGradient(DoFsToPos(cp.gradient(),n_pts));
                //Viewer.showContactForce(DoFsToPos(cp.contactForces(),n_pts));
                Viewer.frameTick();
                std::cout       <<"It: " << i
                                << ", current Energy: " << cp.energy() 
                                << ", Gradient: " << cp.gradient().norm() 
                                << ", Position: " << cp.getVars().norm()
                                << ", Twist Min/Max: " << twist.minCoeff() << " / " <<twist.maxCoeff()
                                << ", Smallest Eigenvalue: " << eigenvalues(0)
                                << std::endl << std::endl;

            }
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
