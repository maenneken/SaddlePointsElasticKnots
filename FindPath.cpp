#include "rrt.h"


int main(int argc, char** argv) {
    std::string start_file = "../data/L400-r0.2-UpTo9Crossings/4_1/0033.obj";
    std::string goal_file = "../data/L400-r0.2-UpTo9Crossings/4_1/0001.obj";
    double rod_radius = 0.2;
    int reductionFactor = 1;
    bool hasCollisions = true;
    int contactStiffness = 10000;


    // Parse command line arguments
    if (argc >= 3) {
        start_file = argv[1];
        goal_file = argv[2];
    }
    if (argc >= 4) {
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
    std::vector<Eigen::Vector3d> start_centerline = reduce_knot_resolution(read_nodes_from_file(start_file), reductionFactor);
    std::vector<Eigen::Vector3d> goal_centerline = reduce_knot_resolution(read_nodes_from_file(goal_file), reductionFactor);

    int n_pts = start_centerline.size();  

    PeriodicRod start_pr = define_periodic_rod(start_centerline,material);
    PeriodicRod goal_pr = define_periodic_rod(goal_centerline,material);

    Eigen::VectorXd start_dofs = start_pr.getDoFs();
    Eigen::VectorXd goal_dofs = goal_pr.getDoFs();

    if(start_dofs.size() != goal_dofs.size()){
        throw std::runtime_error("start and goal are not the same size");
    }

    std::cout << "created pr"<< std::endl;
    
    // 4. Wrap into PeriodicRodList
    PeriodicRodList rod_list = PeriodicRodList(start_pr);
    

    std::cout << "created rodlist"<< std::endl;
    // 5. Setup problem options
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = hasCollisions;
    problemOptions.contactStiffness = contactStiffness;
    problemOptions.dHat = 2 * rod_radius;

    std::cout << "set pr options"<< std::endl;
    
    // 6. Create contact problem
    ContactProblem cp(rod_list, problemOptions);

    //find path by rrt
    std::cout << "Start Energy: " << cp.energy()<< std::endl;
    cp.setVars(goal_dofs);
    std::cout << "Goal Energy: " << cp.energy()<< std::endl;
   


    KnotVisualizer Viewer = KnotVisualizer();
    Viewer.setKnot(start_centerline,0.01 * rod_radius);
    std::vector<Eigen::Vector3d> theta_point = {Eigen::Vector3d(0,0,0), Eigen::Vector3d(0,1000,0)};
    Viewer.setTheta(theta_point, 0.01 * rod_radius);
    //Save first Knot Vars
    auto startKnot = cp.getVars();

    static int iterations = 1000;
    static double maxEnergy = 1000;
    static double stepsize = 0.01;
    static double steplength = 0.1;
    static double goalBias = 0.1;


   
    bool running = false;


    //TODO drehe indexe von goal knot so, dass der abstand zwischen start und goal minmal ist. 0-> muss nicht auf 0 zeigen

    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("Max Energy", &maxEnergy,(0.001),(1),"%.2f");
        ImGui::InputDouble("stepsize", &stepsize,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("stepLength", &steplength,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("goal Bias", &goalBias,(0.001),(0.01),"%.4f");
        


        if (ImGui::Button("Find Path")) {
            running=true;
        }
  
        ImGui::End();   

    });
    
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        if(running){
            RRT rrt(start_dofs,goal_dofs, maxEnergy, steplength, goalBias, stepsize);

            auto path = rrt.findPath(cp,iterations,Viewer);

            std::cout << "Found a path of size: " << path.size() << std::endl; 

            running=false;
        }
    }


    return 0;
}