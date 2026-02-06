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
        //todo add new Energy based on different distance fkt to add information. ex: 1/(d²+1) or something with e^x
        //todo add option for distance fkt where you can decide how many neighbors are included in the calkulation (direct neighbors might not be important to know where they are)
        //todo add gradient and Hessian of the energy fkt
        //Idea is to find negative Eigenvalues
        //todo add controls to load a Knot and set up a contactproblem with all options
        //todo instead of converting to sparse and dense and back use dense.
        ImGui::Begin("Controls");
        ImGui::SliderInt("Step", &path_idx, 0, path.size()-1);
        Viewer.updateKnot(DoFsToPos(path[path_idx], n_pts));
        ImGui::End();   

    });
    
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
    }

    return 0;
}
