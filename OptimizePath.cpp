#include "NEB.h"
#include <algorithm>
#include <thread>   // for std::this_thread::sleep_for
#include <chrono>   // for std::chrono::milliseconds



int main(int argc, char** argv) {
       std::string path_file = "foundPath.txt";
    double rod_radius = 0.2;
    bool hasCollisions = true;
    int contactStiffness = 10000;

    // Parse command line arguments
    if (argc >= 2) {
        path_file = argv[1];
    }
    

    std::vector<double> params = {rod_radius, rod_radius};

    //Create material
    RodMaterial material(
        "ellipse",
        2000,     // Young's modulus
        0.3,      // Poisson's ratio
        params
    );
    std::vector<Eigen::VectorXd> path = loadPathTxt(path_file);

    int n_pts = path[0].size()/4; 
    // Read centerline nodes
    std::vector<Eigen::Vector3d> centerline = DoFsToPos(path[0],n_pts);

     

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

    
    
    //relax ends of path and add them
    auto optimizerOptions = NewtonOptimizerOptions();
    optimizerOptions.niter = 10000;
    optimizerOptions.gradTol = 1e-6;

    std::cout << "Compute true start min" << std::endl;
    cp.setVars(path[0]);
    compute_equilibrium(cp.m_rods,problemOptions,optimizerOptions); 
    Eigen::VectorXd start =  cp.m_rods.getDoFs();

     std::cout << "Compute true goal min" << std::endl;
    cp.setVars(path[path.size()-1]);
    compute_equilibrium(cp.m_rods,problemOptions,optimizerOptions); 
    Eigen::VectorXd goal =  cp.m_rods.getDoFs();

    //recenter data. otherwise the knot is shifted in 3d space
    Eigen::Vector3d start_centroid(0,0,0);
    Eigen::Vector3d goal_centroid(0,0,0);
    for (size_t i = 0; i < n_pts; ++i) {
        start_centroid += start.segment<3>(3*i);
        goal_centroid += goal.segment<3>(3*i);
    }
    start_centroid /= n_pts;
    goal_centroid /= n_pts;

    // 2. Subtract centroid from each point
    for (size_t i = 0; i < n_pts; ++i) {
        start.segment<3>(3*i) -= start_centroid;
        goal.segment<3>(3*i) -= goal_centroid;
    }
    
    //put the true mins in the path
    path.emplace_back(goal);
    path.insert(path.begin(), start);

    static int iterations = 1000;
    static double stepsize = 0.001;
    static double spring_constant = 1;
    int path_idx = 0;


   
    bool running = false;
    bool show_path = false;
    size_t it = 0;
    std::vector<Eigen::Vector3d> Positon = DoFsToPos(path[path_idx], n_pts);
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("stepsize", &stepsize,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("Spring constant", &spring_constant,(0.001),(0.01),"%.4f");

        ImGui::SliderInt("path_idx", &path_idx, 0, path.size()-1);

        path_idx = std::clamp(path_idx, 0, (int)path.size()-1);
        if(!show_path) Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
        


        if (ImGui::Button("Optimize Path")) {
            running=true;
            it=0;
        }
        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            savePathTxt("optimized",path);
        }

  
        ImGui::End();   

    });
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        if(running&& it < iterations){

            nebGradientStep(cp,path,spring_constant,stepsize);
            ++it;
            if(it % 10 == 0){
                cp.setVars(path[path_idx]);
                std::cout       << std::endl << std::endl
                                << BLUE << "It: " << it << RESET
                                << ", current Energy: " << cp.energy() 
                                << GREEN << ", contact Energy: " << cp.contactEnergy() << RESET
                                << ", Gradient: " << cp.gradient().norm()
                                << ", Position: " << cp.getVars().norm() 
                                << std::endl;

            }
            
        }
        if(it >= iterations){
            running=false;
            it = 0;
        }
        if(show_path){
            showPath(path,cp,Viewer);
            show_path = false;
        }
    }


    return 0;
}