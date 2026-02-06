#include <iostream>
#include <vector>

#include "helpers.h"
#include "KnotVisualizer.h"


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

    //Save first Knot Vars
    auto startKnot = cp.getVars();

    static size_t i = 0;
    //set buttons
    int path_idx = 0;
    bool show_path = false;
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
        if(!show_path) Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));

        if(ImGui::Button("show interpolated path")){
            show_path=true;
        }
        ImGui::End();   

    });
    
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        if(show_path){
            for( size_t i = 0; i< path.size()-1;++i){
                Eigen::VectorXd direction = path[i+1] - path[i];
                direction.normalize();
                auto current = path[i];
                while((current - path[i+1]).norm()> 0.1){
                    current += 0.1*direction;
                    Viewer.updateKnot(DoFsToPos(current, n_pts));
                    Viewer.frameTick();
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }
            }
            show_path=false;
        }
        Viewer.frameTick();
    }

    return 0;
}
