#pragma once
#include "helpers.h"
#include "KnotVisualizer.h"

struct rrt_vertex{
    Eigen::VectorXd config;
    int parent;

    rrt_vertex(const Eigen::VectorXd& config_, int parent_)
        : config(config_), parent(parent_) {}
};
class RRT{
    public:
        RRT(const Eigen::VectorXd& start, const Eigen::VectorXd& goal, double maxEnergy = 20, double stepLength = 2, double goalBias = 0.1, double steerStep = 0.1);
        Eigen::VectorXd sampleRandConfig(const Eigen::VectorXd& goal, const std::vector<rrt_vertex>& tree);
        size_t nearestVertex(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree);
        Eigen::VectorXd steerTowardsConfig(ContactProblem& cp, Eigen::VectorXd& near, Eigen::VectorXd& rand);
        std::vector<Eigen::VectorXd> createPath(Eigen::VectorXd connection);
        std::vector<Eigen::VectorXd> findPath(ContactProblem& cp, size_t iterations, KnotVisualizer& Viewer);
        
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


};