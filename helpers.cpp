#include "helpers.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>



/** from 3rdparty/ElasticKnots/python/helpers.py
    def define_periodic_rod(pts, material, rest_curv_rad=np.inf, total_opening_angle=0, minimize_twist=False):
    duplicated_0 = np.linalg.norm(pts[0, :] - pts[-2, :]) < 1e-12
    duplicated_1 = np.linalg.norm(pts[1, :] - pts[-1, :]) < 1e-12
    if not duplicated_0 and not duplicated_1:
        pts = np.vstack((pts, pts[0, :], pts[1, :]))
    elif duplicated_0 != duplicated_1:
        raise ValueError("Only one of the first two nodes was duplicated.")
        
    pr = elastic_rods.PeriodicRod(pts, zeroRestCurvature=True)  # always set rest curvature to zero, then eventually modify restKappas
    pr.setMaterial(material)
    pr.totalOpeningAngle = total_opening_angle
    
    # Set rest curvature
    if rest_curv_rad != np.inf:
        rest_lengths = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        rest_kappas = compute_rest_kappas(rest_curv_rad=rest_curv_rad, rest_lengths=rest_lengths)
        pr.rod.setRestKappas(rest_kappas)
        
    # Set the bending energy type to match the definition from [Bergou et al. 2010]
    # The bending energy in [Bergou et al. 2008] is technically non-physical.
    pr.rod.bendingEnergyType = elastic_rods.BendingEnergyType.Bergou2010
    
    if minimize_twist:
        # Minimize the twisting energy 
        # (i.e. run an optimization on the \theta variables only, 
        # leaving the ends of the rod free to untwist)
        elastic_knots.minimize_twist(pr)
    
    return pr
**/
PeriodicRod define_periodic_rod(std::vector<Eigen::Vector3d> pts, RodMaterial material){
    std::size_t n_pts = pts.size();
    bool duplicated_0 = (pts[0] - pts[n_pts - 2]).norm()  < 1e-12;
    bool duplicated_1 = (pts[1] - pts[n_pts - 1]).norm()  < 1e-12;
    
    if (!duplicated_0 && !duplicated_1){
        pts.push_back(pts[0]);
        pts.push_back(pts[1]);
    }
    else if (duplicated_0 != duplicated_1){
        throw "Only one of the first two nodes was duplicated.";
    }
    PeriodicRod pr(pts,true);
    pr.setMaterial(material);
    pr.setTotalOpeningAngle(0);
    pr.rod.setBendingEnergyType(ElasticRod::BendingEnergyType::Bergou2010);

    return pr;
}

/**
    def read_nodes_from_file(file):
    """
    Supported extensions: obj, txt
    """
    nodes = []
    connectivity = []
    n_rods = 0
    if file.endswith('.obj'):
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith('v'):
                    pt = []
                    for coord in line.split(' ')[1:4]:
                        pt.append(float(coord))
                    nodes.append(np.array(pt))
                if line.startswith('l'):
                    edge = []
                    for index in line.split(' ')[1:3]:
                        edge.append(int(index))
                        if len(edge) == 2 and int(index) < edge[0]: # last edge of a rod 
                            n_rods += 1
                    connectivity.append(edge)
        if n_rods == 1:
            return np.array(nodes)
        elif n_rods > 1:
            indices_connections = [i for i in range(len(connectivity)) if connectivity[i][0] > connectivity[i][1]]
            ne_per_rod = np.append(indices_connections[0] + 1, np.diff(indices_connections))
            pts = np.array(nodes)
            pts_list = [pts[0:ne_per_rod[0], :]]
            for ri in range(0, n_rods-1):
                pts_list.append(pts[indices_connections[ri]:indices_connections[ri+1]])
            return pts_list
    
    elif file.endswith('.txt'):
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                pt = []
                for coord in line.split(' ')[0:3]:
                    pt.append(float(coord))
                nodes.append(np.array(pt))
        return np.array(nodes)
    
    elif not '.' in file.split('/')[0]: # no extension, assum same formatting as .txt
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                pt = []
                for coord in line.split(' ')[0:3]:
                    pt.append(float(coord))
                nodes.append(np.array(pt))
        return np.array(nodes)
 **/
std::vector<Eigen::Vector3d> read_nodes_from_file(std::string &filename){
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);
    
    std::vector<Eigen::Vector3d> nodes;
    std::string line;

    if (filename.find(".obj") != std::string::npos){ //it is a .obj file
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::istringstream stream(line);
            double x, y, z;
            char c;
            stream >> c;
            if(c == 'v' ){ //vertex
                stream >> x >> y >> z;
                nodes.emplace_back(x, y, z);
            }
            else if (c == 'l'){ // connectivity
                continue;   
            }
            else{
                continue; // skip invalid lines
            }
                
        }
        return nodes;
    }
    else if (filename.find(".txt") != std::string::npos){ //it is a .txt file
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::istringstream stream(line);
            double x, y, z;

            stream >> x >> y >> z;
            nodes.emplace_back(x, y, z);

                

        }
        return nodes;
    }
}
//extract the positions out of DoFs. Number of points is needed
std::vector<Eigen::Vector3d> DoFsToPos(Eigen::VectorXd dofs, uint n_pts){
    std::vector<Eigen::Vector3d> pts;
    for (uint i = 0; i < 3 * n_pts; i += 3 ){
        pts.emplace_back(dofs[i],dofs[i+1],dofs[i+2]);
    }
    return pts;
}
//extract the twist out of DoFs. Number of points is needed
Eigen::VectorXd DoFsToTwist(Eigen::VectorXd dofs, uint n_pts){

    Eigen::VectorXd twist = dofs.segment(3*n_pts,n_pts);
 
    return twist;
}

//reduce the knot resolution. Call before defining the knot
std::vector<Eigen::Vector3d> reduce_knot_resolution(std::vector<Eigen::Vector3d> pts, size_t factor){
    std::vector<Eigen::Vector3d> reduced_pts;
    for (uint i = 0; i < pts.size(); i += factor ){
        reduced_pts.emplace_back(pts[i]);
    }
    return reduced_pts;
}
//SuiteSparseMatrix to Eigen::MatrixXd (note wll make it semetrik since only uppertriangle is usaly safed)
Eigen::MatrixXd toEigenDense(SuiteSparseMatrix& ssm){
    TripletMatrix<Triplet<double>> triplets = ssm.getTripletMatrix();
    int n = std::max(triplets.m, triplets.n);
    Eigen::MatrixXd em = Eigen::MatrixXd::Zero(n, n);

    for(size_t k = 0; k< triplets.nnz(); ++k){
        auto &t = triplets.nz[k];
        //we need to sum up, because there can be multiple entries for (i,j)
        em(t.i,t.j) += t.v;
        em(t.j,t.i) += t.v;
    }
    return em;
}
//SuiteSparseMatrix to Eigen::SparseMatrix (note wll make it semetrik since only uppertriangle is usaly safed)
Eigen::SparseMatrix<double> toEigenSparse(SuiteSparseMatrix& ssm){
    TripletMatrix<Triplet<double>> triplets = ssm.getTripletMatrix();
    int n = std::max(triplets.m, triplets.n);
    std::vector<Eigen::Triplet<double>> eig_t;

    for(size_t k = 0; k < triplets.nnz(); ++k){
        auto &t = triplets.nz[k];
        eig_t.emplace_back(t.i,t.j,t.v);
        eig_t.emplace_back(t.j,t.i,t.v);
    }

    Eigen::SparseMatrix<double> em(n,n);
    em.setFromTriplets(eig_t.begin(), eig_t.end());
    return em;
}
//enlarge a sparse matrix, and force it to be symetric
Eigen::SparseMatrix<double> enlargeMatrix(Eigen::SparseMatrix<double>& small, size_t newSize){
    std::vector<Eigen::Triplet<double>> triplets;

    for (size_t k = 0; k < small.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(small, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
            triplets.emplace_back(it.col(), it.row(), it.value());
        }
    }
    Eigen::SparseMatrix<double> big(newSize,newSize);
    big.setFromTriplets(triplets.begin(), triplets.end());

    return big;
}
//use this instead of cp.hessian(). cp.hessian will fail, if contacts exists
Eigen::SparseMatrix<double> computeHessian(ContactProblem& cp){
    bool projectionMask = true;    
    SuiteSparseMatrix result = cp.hessianSparsityPattern();

    //hessian without contactEnergy
    cp.m_rods.hessian(result);

    Eigen::SparseMatrix<double> H_rod = toEigenSparse(result);

    //compute barrier hessian
    const bool projectIPCHessian = projectionMask && cp.m_options.projectContactHessianPSD;
    Eigen::SparseMatrix<double> IPCHessianEigen = cp.m_options.contactStiffness * compute_barrier_potential_hessian(cp.m_collisionMesh, cp.m_rods.deformedPointsMatrix(), cp.m_constraintSet, cp.m_options.dHat, projectIPCHessian);

    Eigen::SparseMatrix<double> resized_IPCHessianEigen = enlargeMatrix(IPCHessianEigen,H_rod.cols());

    return H_rod + resized_IPCHessianEigen;
}

//remove all twist entrys -> new Shape n_pts x n_pts, n_pts
HessianAndGradient removeTwist(Eigen::SparseMatrix<double, 0, int> H_sparse, Eigen::VectorXd g){
    int n = g.size() * 0.75;
    Eigen::SparseMatrix<double, 0, int> H_small(n,n);

    std::vector<Eigen::Triplet<double>> triplets;

    for (int k = 0; k < H_sparse.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H_sparse, k); it; ++it) {
            if (it.row() < n && it.col() < n) {
                triplets.emplace_back(it.row(), it.col(), it.value());
            }
        }
    }

    H_small.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd g_small = g.head(n);

    HessianAndGradient result;
    result.H = H_small;
    result.g = g_small;
    return result;

}

//removes the last entry -> new Shape n-1 x n-1
HessianAndGradient removeTheta(Eigen::SparseMatrix<double, 0, int> H_sparse, Eigen::VectorXd g){
    int n = g.size() -1;
    Eigen::SparseMatrix<double, 0, int> H_small(n,n);

    std::vector<Eigen::Triplet<double>> triplets;

    for (int k = 0; k < H_sparse.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H_sparse, k); it; ++it) {
            if (it.row() < n && it.col() < n) {
                triplets.emplace_back(it.row(), it.col(), it.value());
            }
        }
    }

    H_small.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd g_small = g.head(n);

    HessianAndGradient result;
    result.H = H_small;
    result.g = g_small;
    return result;

}
//TODO make removeTheta and remove Twist not repeat code

