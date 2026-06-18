#include "KnotVisualizer.h"
#include "PCA.h"

struct rrt_vertex{
    Eigen::VectorXd config;
    int parent;
    Eigen::VectorXd projection;
    int unsuccessful_expansions = 0;
    int successful_expansions = 0;
    int depth = 0;

    rrt_vertex(const Eigen::VectorXd& config_, int parent_)
        : config(config_), parent(parent_), depth(0) {}
};
class TreeVisualizer : public KnotVisualizer{
    public:
        TreeVisualizer(const std::string PCA_file = "../data/PCA/25V_Dataset.txt");
        std::vector<std::array<size_t,2>> getEdges(std::vector<rrt_vertex> &tree);
        std::vector<Eigen::Vector3d> projectTree(std::vector<rrt_vertex> &tree, int dim);
        void setTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree);
        void updateTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree);

        PCA pca;
        bool update_tree_visualization = true;
    private:
        polyscope::CurveNetwork* start_tree_curve = nullptr;
        polyscope::CurveNetwork* goal_tree_curve = nullptr;

        std::vector<std::array<size_t,2>> start_edges;
        std::vector<Eigen::Vector3d> start_projection;

        std::vector<std::array<size_t,2>> goal_edges;
        std::vector<Eigen::Vector3d> goal_projection;

        size_t last_start_tree_size = 0;
        size_t last_goal_tree_size = 0;
        
};