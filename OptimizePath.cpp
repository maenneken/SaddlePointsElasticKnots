#include "NEB.h"
#include <algorithm>
#include <thread>   // for std::this_thread::sleep_for
#include <chrono>   // for std::chrono::milliseconds

size_t find_max_energy(ContactProblem& cp, std::vector<Eigen::VectorXd> path){
    size_t max_id = 0;
    double max_energy = 0;
    for (size_t i = 0; i < path.size();++i){
        cp.setVars(path[i]);
        if(cp.energy()> max_energy){
            max_id = i;
            max_energy = cp.energy();
        }
    }
    return max_id;
}

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
    //read path from file
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

    //energy is is not consistent after storing and reaplying vars
    /*
    //relax ends of path and add them
    auto optimizerOptions = NewtonOptimizerOptions();
    optimizerOptions.niter = 10000;
    optimizerOptions.gradTol =  1e-6;

    std::cout << "Compute true start min" << std::endl;
    cp.setVars(path[0]);
    compute_equilibrium(cp.m_rods,problemOptions,optimizerOptions);
    //needed otherwise Vars are wrong 
    cp.updateCachedVars();
    std::cout << GREEN << "Start Gradient: " << cp.gradient().norm() <<" Energy: " << cp.energy() <<" contact Energy: " << cp.contactEnergy() <<RESET<<std::endl;
    Eigen::VectorXd start =  cp.getVars();


    std::cout << "Compute true goal min" << std::endl;
    cp.setVars(path[path.size()-1]);
    compute_equilibrium(cp.m_rods,problemOptions,optimizerOptions); 
    cp.updateCachedVars();
    std::cout << GREEN << "Goal Gradient: " << cp.gradient().norm() <<" Energy: " << cp.energy()<<" contact Energy: " << cp.contactEnergy() <<RESET<<std::endl;
    Eigen::VectorXd goal =  cp.getVars();

    //saving and then cp.setVars() changes gradient and energy
    cp.setVars(start);
    std::cout << GREEN << "Start Gradient after Calling it again: " << cp.gradient().norm()<<" Energy: " << cp.energy() <<" contact Energy: " << cp.contactEnergy()<<RESET<<std::endl;
    cp.setVars(goal);
    std::cout << GREEN << "Goal Gradient after Calling it again " << cp.gradient().norm() <<" Energy: " << cp.energy() <<" contact Energy: " << cp.contactEnergy()<<RESET<<std::endl;

    cp.setVars(path[7]);
    std::cout << RED << "Difference between cp.gradient() and cp.gradient(true): "<<(cp.gradient() -  cp.gradient(true)).norm() << RESET << std::endl;
    
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
    cp.setVars(start);
    std::cout << GREEN << "Start Gradient after Centering: " << cp.gradient().norm()<<" Energy: " << cp.energy() <<RESET<<std::endl;
    cp.setVars(goal);
    std::cout << GREEN << "Goal Gradient after Centering: " << cp.gradient().norm() <<" Energy: " << cp.energy()<<RESET<<std::endl;
    

    //put the true mins in the path
    path.emplace_back(goal);
    path.insert(path.begin(), start);
    */


    //minimize twist
    /*
    for (size_t i = 0; i<path.size(); ++i){
        cp.setVars(path[i]);
        minimize_twist(*cp.m_rods.getRod(0));
    }
    */

    //matricies of path and its gradient
    Eigen::MatrixXd pathMatrix = listToMatrix(path);
    Eigen::MatrixXd gradientMatrix(path.size(), path[0].size());

    static int iterations = 1000;
    static double stepsize = 0.001;
    static double spring_constant = 1;
    int path_idx = 0;
    int old_idx = 0;
    bool globalNebGradient = false;
    bool climbingImage = true;

    //find id of max energy in path
    size_t max_id = find_max_energy(cp,path);
    std::cout << RED << "Max Energy at: " << max_id << RESET << std::endl; 


    bool running = false;
    bool show_path = false;
    size_t it = 0;
    std::vector<Eigen::Vector3d> Positon = DoFsToPos(path[path_idx], n_pts);
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("stepsize", &stepsize,(0.001),(0.01),"%.7f");
        ImGui::InputDouble("Spring constant", &spring_constant,(0.001),(0.01),"%.4f");

        ImGui::SliderInt("path_idx", &path_idx, 0, path.size()-1);
        path_idx = std::clamp(path_idx, 0, (int)path.size()-1);
        if(!show_path && old_idx != path_idx) {
            cp.setVars(path[path_idx]);
            Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
            Viewer.showNodeGradient(DoFsToPos(cp.gradient(true), n_pts));
            old_idx = path_idx;
        }
        ImGui::Checkbox("Global NEB Gradient", &globalNebGradient);
        ImGui::Checkbox("Use climbing Image", &climbingImage);


        if (ImGui::Button("Optimize Path")) {
            running=true;
            it=0;
        }
        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            savePathTxt("optimized.txt",path);
        }

  
        ImGui::End();   

    });

    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        if(running&& it < iterations){

            //global seems to want smaller step size (0.001 and smaller)
            if(globalNebGradient){
                //update gradient Matrix
                for (size_t i = 0; i < path.size();++i){
                cp.setVars(path[i]);
                gradientMatrix.row(i) = cp.gradient(true).transpose(); 
                }

                globalNebGradientStep(gradientMatrix,pathMatrix,spring_constant,stepsize);

            } 
            else if (climbingImage)
            {
                ciNebGradientStep(cp,path,spring_constant,stepsize, max_id);
            }
            
            else{
                //for local step of 0.01 works great
                nebGradientStep(cp,path,spring_constant,stepsize);
            }
            
            //do gradient decend on ends
            //start Knot
            cp.setVars(path[0]);
            path[0] -= stepsize * cp.gradient(true);

            //goal Knot
            cp.setVars(path[path.size()-1]);
            path[path.size()-1] -= stepsize * cp.gradient(true);
    

            ++it;
            if(it % 100 == 0){

                //update path list
                if(globalNebGradient){
                    for (size_t i = 0; i < path.size();++i){
                        path[i]= pathMatrix.row(i).transpose();
                    }
                }
                

                //print Energy...
                cp.setVars(path[path_idx]);
                Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
                Viewer.showNodeGradient(DoFsToPos(cp.gradient(true), n_pts));
                std::cout       << std::endl << std::endl
                                << BLUE << "It: " << it << RESET
                                << ", current Energy: " << cp.energy() 
                                <<", contact Energy: " << cp.contactEnergy()
                                << GREEN << ", Gradient: " << cp.gradient().norm() <<RESET
                                << ", Position: " << cp.getVars().head(3*n_pts).norm()
                                << ", Twist: " << cp.getVars().tail(n_pts +1).norm()
                                << std::endl;

            }
            
        }
        //done running
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