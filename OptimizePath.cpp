#include "NEB.h"
#include "KnotResolution.h"
#include <algorithm>
#include <thread>   // for std::this_thread::sleep_for
#include <chrono>   // for std::chrono::milliseconds

size_t find_max_energy(ContactProblem& cp, std::vector<Eigen::VectorXd> path){
    size_t max_id = 0;
    double max_energy = 0;
    for (size_t i = 1; i < path.size()-1;++i){
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
    int resolutionDoubleFactor = 0;

    // Parse command line arguments
    if (argc >= 2) {
        path_file = argv[1];
    }
    if (argc >= 3){
        resolutionDoubleFactor = std::stod(argv[2]);
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

    //we want to increase the resolution of the path
    if(resolutionDoubleFactor >0){
        //double the resolution
        std::vector<Eigen::VectorXd> doublePath;
        doublePath.reserve(path.size());
        for(int i = 0; i< resolutionDoubleFactor; ++i){
            //for every image in the path
            for(size_t j = 0; j < path.size(); ++j){
                doublePath.emplace_back(doubleKnotResolution(path[j]));
            }
             // update path for next iteration
            path = std::move(doublePath);
        }
       
    }

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
    int current_iteration = 0;
    static double stepsize = 0.001;
    double max_step = 1;
    static double spring_constant = 1;
    int path_idx = 0;
    int old_idx = 0;
    bool globalNebHessian = false;
    bool climbingImage = false;
    bool globalNebLineSearch = true;
    bool globalLBFGSNebLineSearch = false;
    static double max_neb_force_threshold = 0.01;
    float max_neb_force = 1000;

    std::chrono::high_resolution_clock::time_point start;

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
        ImGui::InputDouble("Max NEB Force Threshold", &max_neb_force_threshold,(0.001),(0.01),"%.4f");

        ImGui::SliderInt("path_idx", &path_idx, 0, path.size()-1);
        path_idx = std::clamp(path_idx, 0, (int)path.size()-1);
        if(!show_path && old_idx != path_idx) {
            cp.setVars(path[path_idx]);
            Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
            Viewer.showNodeGradient(DoFsToPos(-cp.gradient(true), n_pts));
            Viewer.showNodeGradientModified(DoFsToPos(F_neb[path_idx], n_pts));
            old_idx = path_idx;
        }
        ImGui::Checkbox("Global Linesearch", &globalNebLineSearch);
        ImGui::Checkbox("Global LBFGS Hessian", &globalNebHessian);
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
            start = std::chrono::high_resolution_clock::now();
            
        }
        if(ImGui::Button("Increase radius")){
            increaseRadius(cp,path);
            Viewer.setKnot(DoFsToPos(path[path_idx], n_pts),0.01 * cp.m_options.dHat/2);
            pathMatrix = listToMatrix(path);
        }
        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            std::string header = "Iterations: " + std::to_string(current_iteration) + ", stepsize: " + std::to_string(stepsize) + ", max step size: " + std::to_string(max_step) + ", spring constant: " + std::to_string(spring_constant) + ", global NEB Hessian: " + (globalNebHessian ? "true" : "false") + ", global NEB Linesearch: " + (globalNebLineSearch ? "true" : "false") + ", climbing image: " + (climbingImage ? "true" : "false");
            std::string filename = path_file;
            current_iteration++;
            size_t dot = filename.find_last_of('.');
            if (dot != std::string::npos) {
                filename.insert(dot, "_" + std::to_string(current_iteration));
            }

            savePathTxt(filename, path, header);
        }

  
        ImGui::End();   

    });

    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        
        if(running && it < iterations){
            current_iteration = it;
            //global seems to want smaller step size (0.001 and smaller)
            if(globalNebLineSearch){
                globalNebLineSearchStep(cp, pathMatrix, spring_constant,climbingImage, max_id,stepsize);
            }   
            else if(globalNebHessian){
                globalNebLBFGSHesStep(cp, pathMatrix, history, spring_constant, climbingImage, max_id, max_step);
            }       
            else{
                //for local step of 0.01 works great
                nebGradientStep(cp, path, F_neb, spring_constant,stepsize, max_id, climbingImage);
                //nebGradientStepLS(cp, path, F_neb, spring_constant, stepsize, max_id, climbingImage);
            }
            
            //do gradient decend on ends
            //start Knot
            // cp.setVars(path[0]);
            // path[0] -= stepsize * cp.gradient(true);

            // //goal Knot
            // cp.setVars(path[path.size()-1]);
            // path[path.size()-1] -= stepsize * cp.gradient(true);
    

            ++it;
            

        
        
            if(it % 100 == 0){

                //update path list
                if(globalNebHessian || globalNebLineSearch || globalLBFGSNebLineSearch){
                     for (size_t i = 0; i < path.size();++i){
                        path[i]= pathMatrix.row(i).transpose();
                    }
                }
                //calk max neb force for all images
                Eigen::VectorXd current_F_neb = Eigen::VectorXd::Zero(path[0].size());
                max_neb_force = 0;
                for (size_t i = 1; i < path.size()-1;++i){
                    cp.setVars(path[i]);
                    if(climbingImage && i == max_id){
                        current_F_neb = climbingForce(cp.gradient(),path[i-1],path[i+1]);
                    }
                    else{
                        current_F_neb = nebForce(spring_constant,cp.gradient(),path[i-1],path[i],path[i+1]);
                    }
                    if(current_F_neb.norm() > max_neb_force){
                        max_neb_force = current_F_neb.norm();
                    }
                }
                

                //print Energy...
                cp.setVars(path[path_idx]);
                Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
                Viewer.showNodeGradient(DoFsToPos(-cp.gradient(true), n_pts));
                current_F_neb = Eigen::VectorXd::Zero(path[0].size());
                if(path_idx > 0 && path_idx < path.size() -1){
                    if(climbingImage && path_idx == max_id){
                        current_F_neb = climbingForce(cp.gradient(),path[path_idx-1],path[path_idx+1]);
                    }
                    else{
                        current_F_neb = nebForce(spring_constant,cp.gradient(),path[path_idx-1],path[path_idx],path[path_idx+1]);
                    }
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
                                << RED << ", Max NEB Force: " << max_neb_force << RESET
                                << std::endl;

            }
            
        }

        //done running
        if(running && (it >= iterations || max_neb_force <= max_neb_force_threshold)){
            auto end = std::chrono::high_resolution_clock::now();
            running = false;
            it = 0;
            std::cout << GREEN << "Finished optimization!" << RESET << std::endl;
            Eigen::VectorXd current_F_neb = Eigen::VectorXd::Zero(path[0].size());
            max_neb_force = 0;
            for (size_t i = 1; i < path.size()-1;++i){
                cp.setVars(path[i]);
                current_F_neb = nebForce(spring_constant,cp.gradient(),path[i-1],path[i],path[i+1]);
                if(current_F_neb.norm() > max_neb_force){
                    max_neb_force = current_F_neb.norm();
                }
                std::cout       << std::endl << std::endl
                                << BLUE << "Image: " << i << RESET
                                << ", current Energy: " << cp.energy() 
                                <<", contact Energy: " << cp.contactEnergy()
                                << GREEN << ", Gradient: " << cp.gradient().norm() << RESET
                                << ", Position: " << cp.getVars().head(3*n_pts).norm()
                                << ", Twist: " << cp.getVars().tail(n_pts +1).norm()
                                << YELLOW <<", NEB Force: " << current_F_neb.norm() << RESET
                                << std::endl;  
            }  
            std::cout << YELLOW << "Max NEB Force: " << max_neb_force << RESET << std::endl;  
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            std::cout << "Optimization took: " << duration << " ms" << std::endl;          
        }
        if(show_path){
            showPath(path,cp,Viewer);
            show_path = false;
        }
    }


    return 0;
}