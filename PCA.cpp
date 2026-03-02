#include "helpers.h"
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>


//does not work because it will not minmize it
Eigen::VectorXd generateKnot(ContactProblem& cp, size_t n_dofs){
    size_t n_pts = n_dofs / 4;
    Eigen::VectorXd rand = Eigen::VectorXd::Zero(n_dofs);
    rand.head(3*n_pts) = Eigen::VectorXd::Random(3*n_pts).normalized() * 50;
    
    //compute equilibirium of rand
    auto optimizerOptions = NewtonOptimizerOptions();
    optimizerOptions.niter = 100;
    optimizerOptions.gradTol =  1e-6;
    
    cp.setVars(rand);
    
    compute_equilibrium(cp.m_rods,cp.m_options,optimizerOptions);

    //needed otherwise Vars are wrong 
    cp.updateCachedVars();
    std::cout << GREEN << "Gradient: " << cp.gradient().norm() <<RESET<<" Energy: " << cp.energy() <<" contact Energy: " << cp.contactEnergy() <<std::endl;
    Eigen::VectorXd relaxed =  cp.getVars();

    return relaxed;
}
std::vector<Eigen::VectorXd> generateKnots(ContactProblem& cp, size_t n, size_t n_dofs){
    std::vector<Eigen::VectorXd> knots;
    for(size_t i = 0; i < n; ++i){
        Eigen::VectorXd knot = generateKnot(cp,n_dofs);
        knots.emplace_back(knot);
    }
    return knots;
    
}

std::vector<Eigen::VectorXd> loadAllKnots(std::string& folder,RodMaterial material, int reductionFactor)
{
    namespace fs = std::filesystem;

    std::vector<fs::path> files;

    // Collect files
    for (const auto& entry : fs::directory_iterator(folder))
    {
        if (entry.is_regular_file())
            files.push_back(entry.path());
    }

    // Deterministic ordering (important for PCA!)
    std::sort(files.begin(), files.end());

    std::vector<Eigen::VectorXd> knots;

    // Load files
    for (const auto& file : files)
    {
        std::cout << "Loading knot: " << file << std::endl;

        auto nodes = read_nodes_from_file(file);

        std::vector<Eigen::Vector3d> centerline = reduce_knot_resolution(nodes, reductionFactor);

        PeriodicRod pr = define_periodic_rod(centerline, material);

        knots.emplace_back(pr.getDoFs());
    }

    return knots;
}
std::vector<Eigen::VectorXd> createDataset(std::string folder,  int reductionFactor){
    std::string file = "../data/CubicLatticeKnots/L100/0_1.txt";
    
    double rod_radius = 0.2;
    std::vector<double> params = {rod_radius, rod_radius};
    // 1. Create material
    RodMaterial material(
        "ellipse",
        2000,     // Young's modulus
        0.3,      // Poisson's ratio
        params
    );

    // 2. Read centerline nodes
    std::vector<Eigen::Vector3d> centerline = reduce_knot_resolution(read_nodes_from_file(file), reductionFactor);


    PeriodicRod pr = define_periodic_rod(centerline,material);

    std::cout << "created pr"<< pr.getDoFs().size()<< std::endl;


    
    // 4. Wrap into PeriodicRodList
    PeriodicRodList rod_list = PeriodicRodList(pr);

    std::cout << "created rodlist"<< std::endl;
    // 5. Setup problem options
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = false;
    problemOptions.contactStiffness = 10000;
    problemOptions.dHat = 2 * rod_radius;

    std::cout << "set pr options"<< std::endl;
    
    // 6. Create contact problem
    ContactProblem cp(rod_list, problemOptions);
    std::cout << "created CP"<< std::endl;
    //load all cubic knots
    std::vector<Eigen::VectorXd> cubic_knots = loadAllKnots(folder,material,reductionFactor);
    savePathTxt("cubicKnots.txt",cubic_knots);


    //gradient decend on the knots
    std::vector<Eigen::VectorXd> relaxed_knots;
    for (size_t i = 0; i< cubic_knots.size();++i){
        cp.setVars(cubic_knots[i]);
        std::cout << "relaxing Knot: " << i << std::endl;
        for (size_t j = 0; j < 100000; ++j){
            Eigen::VectorXd current = cp.getVars();
            Eigen::VectorXd g = cp.gradient();
            current -= 0.001 * g;
            cp.setVars(current);
        }
        relaxed_knots.emplace_back(cp.getVars());
    }
    savePathTxt("relaxedKnots.txt",relaxed_knots);
    return relaxed_knots;
}


int main(int argc, char** argv) {
    std::string folder = "../data/CubicLatticeKnots/L100/";
    int reductionFactor = 4;
    // Parse command line arguments
    if (argc >= 2) {
        folder = argv[1];
    }
    if (argc >= 3) {
        reductionFactor = std::stod(argv[3]);
    }
    if (reductionFactor < 1){
        reductionFactor = 1;
    }
    std::vector<Eigen::VectorXd> knot_list;
    if(false){
        knot_list = createDataset(folder, reductionFactor);
    } else {
        knot_list = loadPathTxt("../data/PCA/25V_Dataset.txt");
    }


    //center data and remove all twist values -> normalise and put into a matrix
    return 0;
}