#include <iostream>
#include <vector>

#include "helpers.h"
#include "KnotVisualizer.h"
#include "KnotResolution.h"

int printNumNegEigenvalues(ContactProblem& cp){
    //look at eigenvalues
    auto H_sparse = computeHessian(cp);
    Eigen::MatrixXd H_dense = Eigen::MatrixXd(H_sparse);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_dense);
    auto u = solver.eigenvalues();
    auto U = solver.eigenvectors();

    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigen decomposition failed\n";
    }
    int e = 0;
    while(u(e) < 0){ //it is a negative eigenvalue
        ++e;
    }
    std::cout << YELLOW <<"|neg eigenvalues| = " << e << RESET << std::endl;
    return e;
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

    //Save first Knot Vars
    auto startKnot = cp.getVars();

    std::vector<Eigen::VectorXd> Knots_with_neg_eigenvalues;

    static size_t i = 0;
    //set buttons
    int path_idx = 0;
    int old_idx = 0;
    bool show_path = false;
    int sleep_mil = 0;
    double step_size = 0.5;
    int intersecion_counter = 0;
    Viewer.setUserCallback([&]() {
        ImGui::Begin("Controls");
        ImGui::PushItemWidth(200);
        ImGui::SliderInt("path_idx", &path_idx, 0, path.size()-1);
        ImGui::PopItemWidth();

        ImGui::SameLine();
        if (ImGui::Button("-")) path_idx--;

        ImGui::SameLine();
        if (ImGui::Button("+")) path_idx++;

        path_idx = std::clamp(path_idx, 0, (int)path.size()-1);
        if(!show_path && path_idx != old_idx){
            Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
            cp.setVars(path[path_idx]);
            std::cout           <<BLUE "path id: " << path_idx << RESET
                                << ", current Energy: " << cp.energy() 
                                << ", contact Energy: " << cp.contactEnergy() 
                                << GREEN << ", Gradient: " << cp.gradient().norm() << RESET
                                << ", Position: " << cp.getVars().head(3*n_pts).norm()
                                << ", Twist: " << cp.getVars().tail(n_pts+1).norm()<< std::endl;

            printNumNegEigenvalues(cp);
            old_idx = path_idx;
        } 

        ImGui::InputInt("sleep between frame", &sleep_mil);
        ImGui::InputDouble("step size", &step_size);
        if(ImGui::Button("show interpolated path")){
            show_path=true;
        }
        if(ImGui::Button("Calk Num neg Eigenvalues")){
            cp.setVars(path[path_idx]);
            printNumNegEigenvalues(cp);
        }
        if(ImGui::Button("Save all Knots with neg eigenvalue")){
            savePathTxt("negativeEigenvalues.txt",Knots_with_neg_eigenvalues);
        }
        if(ImGui::Button("Save Path")){
            savePathTxt("SavedPath.txt",path);
        }
        ImGui::End();   

    });
    
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        if(show_path){
            for( size_t i = path_idx; i< path.size()-1;++i){
                Eigen::VectorXd direction = path[i+1] - path[i];
                direction.normalize();
                auto current = path[i];
                while((current - path[i+1]).norm()> step_size){
                    current += step_size*direction;
                    cp.setVars(current);
                    std::cout       << std::endl << std::endl
                                << BLUE << "It: " << i << RESET
                                << ", current Energy: " << cp.energy() 
                                << GREEN << ", contact Energy: " << cp.contactEnergy() << RESET
                                << ", Gradient: " << cp.gradient().norm()
                                << ", Position: " << cp.getVars().head(3*n_pts).norm()
                                << ", Twist: " << cp.getVars().tail(n_pts+1).norm()
                                << std::endl;

                    if(printNumNegEigenvalues(cp)>0){
                        Knots_with_neg_eigenvalues.push_back(current);
                    }
                    


                    if(cp.contactEnergy() > 1000){
                        std::cout << RED << "possible selfintersection at " << i <<  RESET << std::endl;
                        ++intersecion_counter;
                    }
                    Viewer.updateKnot(DoFsToPos(current, n_pts));
                    Viewer.frameTick();
                    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_mil));
                }
            }
            show_path=false;
            std::cout << RED << "There are  " << intersecion_counter << " intersections" << RESET << std::endl;
        }
        Viewer.frameTick();
    }

    return 0;
}
