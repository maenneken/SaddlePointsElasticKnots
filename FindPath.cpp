#include "rrt.h"
#include <algorithm>
#include <thread>   // for std::this_thread::sleep_for
#include <chrono>   // for std::chrono::milliseconds





int main(int argc, char** argv) {
    //std::string start_file = "../data/L400-r0.2-UpTo9Crossings/4_1/0001.obj";
    //std::string goal_file = "../data/L400-r0.2-UpTo9Crossings/4_1/0033.obj";
    std::string start_file = "../data/NoCollision/reduced0001.obj";
    std::string goal_file = "../data/NoCollision/reduced0033.obj";
    double rod_radius = 0.2;
    int reductionFactor = 4;
    bool hasCollisions = true;
    int contactStiffness = 10000;


    // Parse command line arguments
    if (argc >= 3) {
        start_file = argv[1];
        goal_file = argv[2];
    }
    if (argc >= 4) {
        reductionFactor = std::stod(argv[3]);
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

    double start_energy = cp.contactEnergy();
    std::cout << "Start contact Energy: " << cp.contactEnergy()<< std::endl;
    cp.setVars(goal_dofs);
    double goal_energy = cp.contactEnergy();
    std::cout << "Goal contact Energy: " << cp.contactEnergy()<< std::endl;

    double edgeLength = (start_dofs.segment<3>(0) - start_dofs.segment<3>(3)).norm();
    std::cout << "Start Knot Edge length: " << edgeLength <<std::endl;

    edgeLength = (goal_dofs.segment<3>(0) - goal_dofs.segment<3>(3)).norm();
    std::cout << "Goal Knot Edge length: " << edgeLength <<std::endl;



   


    TreeVisualizer Viewer = TreeVisualizer();
    Viewer.setKnot(start_centerline,0.01 * rod_radius);

    //show the goal state
    Viewer.setGoalKnot(goal_centerline,0.01 * rod_radius);
    //Save first Knot Vars
    auto startKnot = cp.getVars();

    static int iterations = 5000;
    static int pruningInterval = 1000;
    static double maxEnergy = std::max({start_energy,goal_energy})*2 +1;
    static double stepsize = 1;
    static double steplength = 10;
    static double goalBias = 0.2;
    static bool oneRandDirection = true;
    static double constraintStiffness = 1;

   
    bool running = false;
    bool show_path = false;
    std::vector<Eigen::VectorXd> path;
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputInt("Pruning Interval", &pruningInterval,1000,10000);
        ImGui::InputDouble("Max Energy", &maxEnergy,(0.001),(1),"%.2f");
        ImGui::InputDouble("stepsize", &stepsize,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("stepLength", &steplength,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("goal Bias", &goalBias,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("Constraint Stiffness", &constraintStiffness,(0.001),(0.01),"%.4f");
        ImGui::Checkbox("oneRandDirection", &oneRandDirection);
    
        


        if (ImGui::Button("Find Path")) {
            running=true;
        }
        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            savePathTxt("foundPath.txt",path);
        }

  
        ImGui::End();   

    });
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        if(running){
            RRT rrt(start_dofs, goal_dofs, maxEnergy, steplength, goalBias, stepsize, pruningInterval,oneRandDirection,constraintStiffness);

            path = rrt.findConstrainedPath(cp,iterations,Viewer);
            //path = rrt.findPath(cp,iterations,Viewer);

            std::cout << "Found a path of size: " << path.size() << std::endl; 
            showPath(path,cp,Viewer);
            running=false;
        }
        if(show_path){
            showPath(path,cp,Viewer);
            show_path = false;
        }

    }


    return 0;
}