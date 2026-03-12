#pragma once
#include "helpers.h"
#include "TreeVisualizer.h"
#include "projectToConstraintSpace.h"

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RESET   "\033[0m"

class RRT{
    public:
        RRT(const Eigen::VectorXd& start, const Eigen::VectorXd& goal,PCA _pca, double maxEnergy = 20, double stepLength = 2, double goalBias = 0.1, double steerStep = 0.1, size_t pruningInterval = 100, bool oneRandDirection = true,double _radius= 2);
        std::vector<rrt_vertex> pruneAllLeafNodes(std::vector<rrt_vertex>& tree);
        Eigen::VectorXd sampleRandConfig(ContactProblem& cp, const Eigen::VectorXd& goal, const std::vector<rrt_vertex>& tree);
        Eigen::VectorXd sampleRandDirection(ContactProblem& cp, const Eigen::VectorXd& current_config, const Eigen::VectorXd& goal);
        size_t nearestVertex(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree);
        std::vector<size_t> kNearestNeighbor(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree, int k);
        std::vector<size_t> radiusNeighbor(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree, double r);
        void updateNeighboringWeights(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree,std::vector<double>& weights, double r);
        std::vector<double> updateAllWeights( const std::vector<rrt_vertex>& tree, double r);
        Eigen::VectorXd steerTowardsConfig(ContactProblem& cp, Eigen::VectorXd& near, Eigen::VectorXd& rand);
        Eigen::VectorXd steerInDirection(ContactProblem& cp, Eigen::VectorXd& config, Eigen::VectorXd& direction);
        std::vector<Eigen::VectorXd> createPath(Eigen::VectorXd connection);
        void expand(ContactProblem& cp, size_t config_id, const Eigen::VectorXd& goal, std::vector<rrt_vertex>& tree,std::vector<double>& weights, size_t k);
        bool connect(ContactProblem& cp, std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree, double l);
        std::vector<Eigen::VectorXd> findPath(ContactProblem& cp, size_t iterations, TreeVisualizer& Viewer);
        std::vector<Eigen::VectorXd> findConstrainedPath(ContactProblem& cp, size_t iterations, TreeVisualizer& Viewer);
        
        std::vector<rrt_vertex> start_tree;
        std::vector<rrt_vertex> goal_tree;

        std::vector<double> start_weight;
        std::vector<double> goal_weight;
        
    private:
        double step_length;
        double max_energy;
        double min_val;
        double max_val;
        size_t n_dofs;
        double goal_bias;
        double steer_step;
        size_t pruning_interval;
        bool one_rand_direction_3d;
        double constraint_stiffness;
        double radius;

        PCA pca;
        


};