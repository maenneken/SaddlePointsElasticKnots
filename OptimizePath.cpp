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
void savePath(std::string path_file, std::vector<Eigen::VectorXd> path, int current_iteration, std::string header){
    std::string filename = path_file;
    size_t dot = filename.find_last_of('.');
    if (dot != std::string::npos) {
        filename.insert(dot, "_" + std::to_string(current_iteration));
    }
    savePathTxt(filename, path, header);
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

    static int iterations = 1000000;
    int current_iteration = 0;
    static double stepsize = 0.0001;
    double max_step = 1;
    static double spring_constant = 0.00001;
    int path_idx = 0;
    int old_idx = 0;
    bool globalNebHessian = false;
    bool climbingImage = false;
    bool globalNebLineSearch = true;
    bool globalLBFGSNebLineSearch = false;
    static double max_neb_force_threshold = 1.0;
    float max_neb_force = 1000;
    double new_radius = rod_radius;

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

    bool running = true;
    bool benchmarking = true;
    bool show_path = false;
    std::vector<double> benchmark_max_neb_forces = {1.0, 0.5, 0.1, 0.05, 0.01};
    size_t benchmark_run = 0;
    size_t it = 0;
    std::vector<Eigen::Vector3d> Positon = DoFsToPos(path[path_idx], n_pts);
    

    auto start_knot = path[0];
    auto goal_knot = path.back();
    relax_start_goal(cp, start_knot, goal_knot);

    path.insert(path.begin(), start_knot);
    path.push_back(goal_knot);
    Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));

    pathMatrix = listToMatrix(path);
    F_neb.resize(path.size(), Eigen::VectorXd::Zero(path[0].size()));
    d_search.resize(path.size(), Eigen::VectorXd::Zero(path[0].size()));


    //set buttons
    start = std::chrono::high_resolution_clock::now();
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputDouble("stepsize", &stepsize,(0.001),(0.01),"%.7f");
        ImGui::InputDouble("max step size", &max_step,(0.001),(0.01),"%.7f");
        ImGui::InputDouble("Spring constant", &spring_constant,(0.001),(0.01),"%.7f");
        ImGui::InputDouble("Max NEB Force Threshold", &max_neb_force_threshold,(0.001),(0.01),"%.7f");
                

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
        if(ImGui::Button("Benchmark")){
            running=true;
            benchmarking = true;
            it=0;
            max_neb_force_threshold = benchmark_max_neb_forces[0];
            start = std::chrono::high_resolution_clock::now();
        }
        ImGui::InputDouble("new radius", &new_radius);
        if(ImGui::Button("Increase radius")){
            while(new_radius >= cp.m_options.dHat/2){
                increaseRadius(cp,path);
                Viewer.setKnot(DoFsToPos(path[path_idx], n_pts),0.01 * cp.m_options.dHat/2);
            } 
            pathMatrix = listToMatrix(path);
        }
        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            std::string header = "Iterations: " + std::to_string(current_iteration) + ", stepsize: " + std::to_string(stepsize) + ", max step size: " + std::to_string(max_step) + ", spring constant: " + std::to_string(spring_constant) + ", global NEB Hessian: " + (globalNebHessian ? "true" : "false") + ", global NEB Linesearch: " + (globalNebLineSearch ? "true" : "false") + ", climbing image: " + (climbingImage ? "true" : "false");
            savePath(path_file, path, current_iteration, header);
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

                    //do gradient decend on ends
                    //start Knot
            }
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
                double max_spring_force = 0;
                Eigen::VectorXd current_F_spring = Eigen::VectorXd::Zero(path[0].size());
                for (size_t i = 1; i < path.size()-1;++i){
                    cp.setVars(path[i]);
                    current_F_spring = springForce(spring_constant,path[i-1],path[i],path[i+1],tangent(path[i-1],path[i+1]));
                    if(current_F_spring.norm()> max_spring_force){
                        max_spring_force = current_F_spring.norm();
                    }
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
                                << RED << ", MAx spring Force: " << max_spring_force << RESET
                                << std::endl;

            }
            
        }

        //done running
        if(running && (it >= iterations || max_neb_force <= max_neb_force_threshold)){
            auto end = std::chrono::high_resolution_clock::now();
            running = false;
            std::vector<double> energies;
            std::vector<double> neb_forces;
            std::vector<double> gradients;
            std::vector<double> spring_forces;
            std::vector<double> bend_energies;
            std::vector<double> stretch_energies;
            std::vector<double> twist_energies;
            std::cout << GREEN << "Finished optimization!" << RESET << std::endl;
            Eigen::VectorXd current_F_neb = Eigen::VectorXd::Zero(path[0].size());
            max_neb_force = 0;
            cp.setVars(path[0]);
            gradients.push_back(cp.gradient(true).norm());
            energies.push_back(cp.energy());
            neb_forces.push_back(0);
            spring_forces.push_back(0);
            

            bend_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Bend));
            stretch_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Stretch));
            twist_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Twist));

            for (size_t i = 1; i < path.size()-1;++i){
                cp.setVars(path[i]);
                gradients.push_back(cp.gradient(true).norm());
                current_F_neb = nebForce(spring_constant,cp.gradient(),path[i-1],path[i],path[i+1]);
                neb_forces.push_back(current_F_neb.norm());
                energies.push_back(cp.energy());
            
                
                bend_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Bend));
                stretch_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Stretch));
                twist_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Twist));

                auto current_F_spring = springForce(spring_constant,path[i-1],path[i],path[i+1],tangent(path[i-1],path[i+1]));
                spring_forces.push_back(current_F_spring.norm());
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
            cp.setVars(path[path.size()-1]);
            gradients.push_back(cp.gradient(true).norm());
            energies.push_back(cp.energy());
            neb_forces.push_back(0);
            spring_forces.push_back(0);

            bend_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Bend));
            stretch_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Stretch));
            twist_energies.push_back(cp.m_rods.energy(PeriodicRod::EnergyType::Twist));

            double neb_energy = nebEnergy(cp,path,spring_constant);
            std::cout << GREEN << "Optimization finished after " << current_iteration << " iterations." << RESET << std::endl; 
            std::cout << YELLOW << "NEB Energy: " << neb_energy << RESET << std::endl; 
            std::cout << YELLOW << "Max NEB Force: " << max_neb_force << RESET << std::endl;  
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            std::cout << "Optimization took: " << duration << " ms" << std::endl; 
            
            if(true){
                auto now = std::chrono::system_clock::now();
                auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( now.time_since_epoch()).count();
                std::string filename = path_file;
                current_iteration++;
                size_t dot = filename.find_last_of('.');
                if (dot != std::string::npos) {
                    filename.insert(dot, "_" + std::to_string(current_iteration) + "_optimizationResult" );
                }
                std::ofstream optimization_result(filename);
                optimization_result << "Optimization finished after " << current_iteration << " iterations." << std::endl; 
                optimization_result << "NEB Energy: " << neb_energy << std::endl; 
                optimization_result << "Max NEB Force: " << max_neb_force << std::endl;  
                optimization_result << "Optimization took: " << duration << " ms" << std::endl; 
                optimization_result << "Energies: " << std::endl;
                for(size_t i = 0; i < energies.size();++i){
                    optimization_result << energies[i] << ", ";
                }
                optimization_result << std::endl;

                optimization_result << "NEB Forces: " << std::endl;
                for(size_t i = 0; i < neb_forces.size();++i){
                    optimization_result << neb_forces[i] << ", ";
                }
                optimization_result << std::endl;

                optimization_result << "Gradients: " << std::endl;
                for(size_t i = 0; i < gradients.size();++i){
                    optimization_result << gradients[i] << ", "; 
                }
                optimization_result << std::endl;

                optimization_result << "Spring Forces: " << std::endl;
                for(size_t i = 0; i < spring_forces.size();++i){
                    optimization_result << spring_forces[i] << ", ";
                }
                optimization_result << std::endl;

                optimization_result << "Bend Energies: " << std::endl;
                for(size_t i = 0; i < bend_energies.size();++i){
                    optimization_result << bend_energies[i] << ", ";
                }
                optimization_result << std::endl;

                optimization_result << "Twist Energies: " << std::endl;
                for(size_t i = 0; i < twist_energies.size();++i){
                    optimization_result << twist_energies[i] << ", ";
                }
                optimization_result << std::endl;

                optimization_result << "Stretch Energies: " << std::endl;
                for(size_t i = 0; i < stretch_energies.size();++i){
                    optimization_result << stretch_energies[i] << ", ";
                }
                optimization_result << std::endl;



                std::string header = "Iterations: " + std::to_string(current_iteration) + ", stepsize: " + std::to_string(stepsize) + ", max step size: " + std::to_string(max_step) + ", spring constant: " + std::to_string(spring_constant) + ", global NEB Hessian: " + (globalNebHessian ? "true" : "false") + ", global NEB Linesearch: " + (globalNebLineSearch ? "true" : "false") + ", climbing image: " + (climbingImage ? "true" : "false");
                savePath(path_file,path,current_iteration,header);
            }

            if(benchmarking){
                running = true;
                benchmark_run++;
                if(benchmark_run >= benchmark_max_neb_forces.size()){
                    benchmarking = false;
                    std::exit(0);
                }
                max_neb_force_threshold = benchmark_max_neb_forces[benchmark_run];
            }

        }
        if(show_path){
            showPath(path,cp,Viewer);
            show_path = false;
        }
    }

    return 0;
}