#pragma once
#include "helpers.h"
#include "TreeVisualizer.h"
#include "projectToConstraintSpace.h"
#include <nanoflann.hpp> 

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RESET   "\033[0m"



struct RRTAdaptor {
    const std::vector<rrt_vertex>& pts;

    RRTAdaptor(const std::vector<rrt_vertex>& pts_) : pts(pts_) {}

    // number of points
    inline size_t kdtree_get_point_count() const {
        return pts.size();
    }

    // coordinate accessor
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return pts[idx].projection(dim);
    }

    // optional bounding box
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, RRTAdaptor>,
    RRTAdaptor,
    -1,
    size_t
>;

class RRT{
    public:
        RRT(const Eigen::VectorXd& start, const Eigen::VectorXd& goal,PCA _pca,bool all_permutations = true, size_t _every_k_permutation = 1);
        void add_all_permutations(const Eigen::VectorXd& start, const Eigen::VectorXd& goal);
        std::vector<rrt_vertex> pruneAllLeafNodes(std::vector<rrt_vertex>& tree);
        Eigen::VectorXd sampleRandConfig(ContactProblem& cp, const Eigen::VectorXd& goal, const std::vector<rrt_vertex>& tree);
        Eigen::VectorXd sampleRandDirection(ContactProblem& cp, const Eigen::VectorXd& current_config, const Eigen::VectorXd& goal);
        size_t nearestVertex(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree);
        void buildKDTrees();
        std::vector<size_t> kNearestNeighbor(const Eigen::VectorXd& query_projection, const KDTree& index, int k);
        std::vector<size_t> radiusNeighbor(const Eigen::VectorXd& query_projection, const KDTree& index, double r);
        void updateNeighboringWeights(const Eigen::VectorXd& config,std::vector<rrt_vertex>& tree,std::vector<double>& weights, const KDTree& kd_tree, double r);
        std::vector<double> updateAllWeights( const std::vector<rrt_vertex>& tree, KDTree& kd_tree, double r);
        Eigen::VectorXd steerTowardsConfig(ContactProblem& cp, Eigen::VectorXd& near, Eigen::VectorXd& rand);
        Eigen::VectorXd steerInDirection(ContactProblem& cp, Eigen::VectorXd& config, Eigen::VectorXd& direction);
        std::vector<Eigen::VectorXd> createPath(Eigen::VectorXd connection);
        void expand(ContactProblem& cp, size_t config_id, const Eigen::VectorXd& goal, std::vector<rrt_vertex>& tree,std::vector<double>& weights, KDTree& kd_tree, size_t k);
        bool connect(ContactProblem& cp, double l);
        void printTreeQuota();
        std::vector<Eigen::VectorXd> findConstrainedPath(ContactProblem& cp, size_t iterations, TreeVisualizer& Viewer);
        
        std::vector<rrt_vertex> start_tree;
        std::vector<rrt_vertex> goal_tree;

        std::vector<double> start_weight;
        std::vector<double> goal_weight;


        double step_length = 40;
        double max_energy = 10;
        double goal_bias = 0.2;
        double steer_step = 1;
        size_t pruning_interval = 1000;
        bool one_rand_direction_3d = false;
        double radius = 10;
        size_t k_neighbors = 10; //for goal bias 
        size_t k_expension = 1; //number of random directions to sample for each expansion
        double rotation_bias = 0.1;
        double gradient_bias = 0.1;
        size_t proj_dim = 15;
        size_t sample_proj_dim = 30;
        bool goal_bias_for_all_permutations = false;
        bool sample_in_projection_space = false;
        bool use_constraint_projection_for_sampling = false;
        bool reproject_direction = false;
        size_t every_k_permutation = 1;
        double max_rod_energy = 30000;
        bool start_with_rotation = false;
    private:
        double min_val;
        double max_val;
        size_t n_dofs;
        size_t n_permutations = 0;
        
        PCA pca;
        
        std::unique_ptr<RRTAdaptor> start_adaptor;
        std::unique_ptr<RRTAdaptor> goal_adaptor;

        std::unique_ptr<KDTree> start_kd_tree;
        std::unique_ptr<KDTree> goal_kd_tree;


};