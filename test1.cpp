#include "helpers.h"
#include <iostream>
#include <vector>


#include "KnotVisualizer.h"
Eigen::MatrixXd readHessian(const std::string& filename, int n) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Cannot open file");
    }

    Eigen::MatrixXd H(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            in >> H(i, j);
        }
    }
    return H;
}

int main() {
    std::string file = "../data/L400-r0.2-UpTo9Crossings/4_1/0033.obj";
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
    std::vector<Eigen::Vector3d> centerline = reduce_knot_resolution(read_nodes_from_file(file), 4);

    //Set Visulizer
    KnotVisualizer Viewer = KnotVisualizer();
    Viewer.setKnot(centerline,0.01 * rod_radius);
    Viewer.frameTick();


    PeriodicRod pr = define_periodic_rod(centerline,material);

    std::cout << "created pr"<< std::endl;


    
    // 4. Wrap into PeriodicRodList
    PeriodicRodList rod_list = PeriodicRodList(pr);

    std::cout << "created rodlist"<< std::endl;
    // 5. Setup problem options
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = true;
    problemOptions.contactStiffness = 10000;
    problemOptions.dHat = 2 * rod_radius;

    std::cout << "set pr options"<< std::endl;
    
    // 6. Create contact problem
    ContactProblem cp(rod_list, problemOptions);
    std::cout << "created CP"<< std::endl;

    // 7. Output gradient and DoF sizes
    std::cout << "Gradient length: " << cp.gradient().size()
              << ", DoFs length: " << pr.getDoFs().size()
              << ", same as Vars: " << cp.getVars().size()
              << ", Energy: " << cp.energy()
              << ", contact Energy: " << cp.contactEnergy()
              << std::endl;


    SuiteSparseMatrix H = cp.hessian();
    auto H_sparse = toEigenSparse(H);

    
    Eigen::MatrixXd H_dense = toEigenDense(H);
    Eigen::MatrixXd H_dense_from_sparse = Eigen::MatrixXd(H_sparse);

    std::cout << "Diff between toEigenDense(H) and Eigen::MatrixXd(toEigenSparse(H)):" << (H_dense - H_dense_from_sparse).norm() << std::endl;


    double asym = (H_dense - H_dense.transpose()).norm();
    std::cout << "asymmetry = " << asym << std::endl;

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
    std::cout << "|neg eigenvalues| = " << e << std::endl;

    Eigen::MatrixXd H_dense_numpy = readHessian("../hessian.txt",cp.getVars().size());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> npsolver(H_dense_numpy );
    auto u_np = npsolver.eigenvalues();
    auto U_np = npsolver.eigenvectors();

    if (npsolver.info() != Eigen::Success) {
        std::cerr << "Eigen decomposition failed\n";
    }
    e = 0;
    while(u_np(e) < 0){ //it is a negative eigenvalue
        ++e;
    }
    std::cout << "|Numpy Hessian neg eigenvalues| = " << e << std::endl;

    std::cout << "Diffrence between Numpy Hessian and Eigen hessian:" << (H_dense - H_dense_numpy).norm() << std::endl;
    std::cout << "the diffrence exists because the np Hessian was not constructed correctly" << std::endl;
    return 0;
}
