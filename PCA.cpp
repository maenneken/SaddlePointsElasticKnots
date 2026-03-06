#include "PCA.h"

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
    std::vector<Eigen::VectorXd> relaxedknots;
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

std::vector<Eigen::VectorXd> removeTwist(std::vector<Eigen::VectorXd> &data){
    size_t n_pts = data[0].size() /4;
    std::vector<Eigen::VectorXd> noTwist;
    for(Eigen::VectorXd current : data){
        noTwist.emplace_back(current.head(3*n_pts));
    }
    return noTwist;
}
Eigen::MatrixXd toMatrix(const std::vector<Eigen::VectorXd>& data)
{
    size_t N = data.size();
    size_t d = data[0].size();

    Eigen::MatrixXd M(N, d);

    for (size_t i = 0; i < N; ++i)
        M.row(i) = data[i].transpose();

    return M;
}



void PCA::fit(const std::vector<Eigen::VectorXd>& samples){
    Eigen::MatrixXd X = toMatrix(samples);

    mean = X.colwise().mean();
    Eigen::MatrixXd centered = X.rowwise() - mean;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            centered,
            Eigen::ComputeThinV
    );

    components = svd.matrixV();   // all components
}

Eigen::VectorXd PCA::project(const Eigen::VectorXd& v, int k = 3){
    Eigen::VectorXd centered = v - mean.transpose();
    return components.leftCols(k).transpose() * centered;
}


