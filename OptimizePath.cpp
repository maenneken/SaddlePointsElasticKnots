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

    //matricies of path and its gradient
    RowMatrix pathMatrix = listToMatrix(path);

    static int iterations = 1000;
    static double stepsize = 0.001;
    static double max_step = 1;
    static double spring_constant = 100;
    int path_idx = 0;
    int old_idx = 0;
    bool globalNeb = true;
    bool climbingImage = false;

    //find id of max energy in path
    size_t max_id = find_max_energy(cp,path);
    std::cout << RED << "Max Energy at: " << max_id << RESET << std::endl; 

    //Force Vector for visualistation
    std::vector<Eigen::VectorXd> F_neb;
    F_neb.resize(path.size(), Eigen::VectorXd::Zero(path[0].size()));

    std::vector<Eigen::VectorXd> d_search;
    d_search.resize(path.size(), Eigen::VectorXd::Zero(path[0].size()));

    LBGFhistory history;

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
        ImGui::InputDouble("max step size", &max_step,(0.001),(0.01),"%.7f");
        ImGui::InputDouble("Spring constant", &spring_constant,(0.001),(0.01),"%.4f");

        ImGui::SliderInt("path_idx", &path_idx, 0, path.size()-1);
        path_idx = std::clamp(path_idx, 0, (int)path.size()-1);
        if(!show_path && old_idx != path_idx) {
            cp.setVars(path[path_idx]);
            Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
            Viewer.showNodeGradient(DoFsToPos(-cp.gradient(true), n_pts));
            Viewer.showNodeGradientModified(DoFsToPos(F_neb[path_idx], n_pts));
            old_idx = path_idx;
        }
        ImGui::Checkbox("Global NEB Gradient", &globalNeb);
        ImGui::Checkbox("Use climbing Image", &climbingImage);

        if(ImGui::Button("Relax start and goal")){
            auto start = path[0];
            auto goal = path.back();
            relax_start_goal(cp, start, goal);

            path.insert(path.begin(), start);
            path.push_back(goal);
            Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));

            pathMatrix = listToMatrix(path);
            F_neb.resize(path.size(), Eigen::VectorXd::Zero(path[0].size()));
            d_search.resize(path.size(), Eigen::VectorXd::Zero(path[0].size()));
            
        }
        if (ImGui::Button("Optimize Path")) {
            running=true;
            it=0;
        }
        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            std::string header = "Iterations: " + std::to_string(iterations) + ", stepsize: " + std::to_string(stepsize) + ", max step size: " + std::to_string(max_step) + ", spring constant: " + std::to_string(spring_constant) + ", global NEB: " + (globalNeb ? "true" : "false") + ", climbing image: " + (climbingImage ? "true" : "false");
            savePathTxt("optimized.txt",path);
        }

  
        ImGui::End();   

    });

    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        if(running&& it < iterations){

            //global seems to want smaller step size (0.001 and smaller)
            if(globalNeb){
                //update gradient Matrix
                globalNebLBFGSStep(cp,pathMatrix,history,spring_constant,climbingImage,max_id,max_step);

            } 
            else if (climbingImage)
            {
                ciNebGradientStep(cp,path,F_neb, spring_constant,stepsize, max_id);
            }
            
            else{
                //for local step of 0.01 works great
                //nebGradientStep(cp,path,F_neb,spring_constant,stepsize);
                nebGradientStep(cp, path, F_neb, spring_constant,stepsize);
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
                if(globalNeb){
                    for (size_t i = 0; i < path.size();++i){
                        path[i]= pathMatrix.row(i).transpose();
                    }
                }
                

                //print Energy...
                cp.setVars(path[path_idx]);
                Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
                Viewer.showNodeGradient(DoFsToPos(-cp.gradient(true), n_pts));
                Eigen::VectorXd current_F_neb = Eigen::VectorXd::Zero(path[0].size());
                if(path_idx > 0 && path_idx < path.size() -1){
                    current_F_neb = nebForce(spring_constant,cp.gradient(),path[path_idx-1],path[path_idx],path[path_idx+1]);
                }
                Viewer.showNodeGradientModified(DoFsToPos(current_F_neb, n_pts));
                std::cout       << std::endl << std::endl
                                << BLUE << "It: " << it << RESET
                                << ", current Energy: " << cp.energy() 
                                <<", contact Energy: " << cp.contactEnergy()
                                << GREEN << ", Gradient: " << cp.gradient().norm() << RESET
                                << ", Position: " << cp.getVars().head(3*n_pts).norm()
                                << ", Twist: " << cp.getVars().tail(n_pts +1).norm()
                                << YELLOW <<", NEB Force: " << current_F_neb.norm() << RESET
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