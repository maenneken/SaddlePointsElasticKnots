#pragma once
#include "helpers.h"
#include "KnotVisualizer.h"
#include "projectToConstraintSpace.h"

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define RESET   "\033[0m"

struct rrt_vertex{
    Eigen::VectorXd config;
    size_t parent;

    rrt_vertex(const Eigen::VectorXd& config_, size_t parent_)
        : config(config_), parent(parent_) {}
};
class RRT{
    public:
        RRT(const Eigen::VectorXd& start, const Eigen::VectorXd& goal, double maxEnergy = 20, double stepLength = 2, double goalBias = 0.1, double steerStep = 0.1, size_t pruningInterval = 100);
        std::vector<rrt_vertex> pruneAllLeafNodes(std::vector<rrt_vertex>& tree);
        Eigen::VectorXd sampleRandConfig(ContactProblem& cp, const Eigen::VectorXd& goal, const std::vector<rrt_vertex>& tree);
        Eigen::VectorXd sampleRandDirection(const Eigen::VectorXd& current_config, const Eigen::VectorXd& goal);
        size_t nearestVertex(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree);
        Eigen::VectorXd steerTowardsConfig(ContactProblem& cp, Eigen::VectorXd& near, Eigen::VectorXd& rand);
        Eigen::VectorXd steerInDirection(ContactProblem& cp, Eigen::VectorXd& config, Eigen::VectorXd& direction);
        std::vector<Eigen::VectorXd> createPath(Eigen::VectorXd connection);
        std::vector<Eigen::VectorXd> findPath(ContactProblem& cp, size_t iterations, KnotVisualizer& Viewer);
        std::vector<Eigen::VectorXd> findConstrainedPath(ContactProblem& cp, size_t iterations, KnotVisualizer& Viewer);
        
    private:
        double step_length;
        double max_energy;
        std::vector<rrt_vertex> start_tree;
        std::vector<rrt_vertex> goal_tree;
        double min_val;
        double max_val;
        size_t n_dofs;
        double goal_bias;
        double steer_step;
        size_t pruning_interval;


};