#include "KnotVisualizer.h"
#include "PCA.h"

struct rrt_vertex{
    Eigen::VectorXd config;
    int parent;

    rrt_vertex(const Eigen::VectorXd& config_, int parent_)
        : config(config_), parent(parent_) {}
};
class TreeVisualizer : public KnotVisualizer{
    public:
        TreeVisualizer(std::string PCA_file = "../data/PCA/25V_Dataset.txt");

        void setTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree);
        void updateTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree);
    private:
        polyscope::CurveNetwork* start_tree_curve = nullptr;
        polyscope::CurveNetwork* goal_tree_curve = nullptr;
        PCA pca;
};