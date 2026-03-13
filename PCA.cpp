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

std::vector<Eigen::VectorXd> removeTwist(const std::vector<Eigen::VectorXd> &data){
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
    data_dim = samples[0].size();
    Eigen::MatrixXd X = toMatrix(samples);

    mean = X.colwise().mean();
    Eigen::MatrixXd centered = X.rowwise() - mean;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            centered,
            Eigen::ComputeThinV
    );
    
    components = svd.matrixV();   // all components

    // store singular values
    auto singular_values = svd.singularValues();

    // compute explained variance
    Eigen::VectorXd var = singular_values.array().square();

    double total_var = var.sum();

    auto variance_ratio = var / total_var;

    for(int i = 0; i < variance_ratio.size(); ++i){
    std::cout << "PC" << i+1 << ": "
              << variance_ratio(i) * 100
              << "%\n";
    }
}

Eigen::VectorXd PCA::project(const Eigen::VectorXd& v, int k){
    Eigen::VectorXd centered;
    if(v.size() > data_dim){
        centered = removeTwist({v})[0] - mean.transpose();
    } else {
        centered = v - mean.transpose();
    }
   
    return components.leftCols(k).transpose() * centered;
}

Eigen::VectorXd PCA::reconstruct(const Eigen::VectorXd& y, int k){
    return mean.transpose() + components.leftCols(k) * y;
}

double PCA::dist (const Eigen::VectorXd& x, const Eigen::VectorXd& y, int k){
    Eigen::VectorXd x_short = removeTwist({x})[0];
    Eigen::VectorXd y_short = removeTwist({y})[0];

    Eigen::VectorXd x_proj = project(x_short,k);
    Eigen::VectorXd y_proj = project(y_short,k);

    double d = (x_proj - y_proj).squaredNorm();  

    return d;
}


