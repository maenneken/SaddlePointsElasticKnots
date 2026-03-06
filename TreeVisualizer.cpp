#include "TreeVisualizer.h"

TreeVisualizer::TreeVisualizer(std::string file){
    std::vector<Eigen::VectorXd> samples = removeTwist(loadPathTxt(file));
    pca.fit(samples);
}
std::vector<std::array<size_t,2>> getEdges(std::vector<rrt_vertex> &tree){
    std::vector<std::array<size_t,2>> edges;
    for (size_t i = 0; i < tree.size(); ++i){
        //it has a parent
        if(tree[i].parent > 0){
            edges.emplace_back({tree[i].parent, i});
        }
    }
    return edges;
}
std::vector<Eigen::Vector3d> projectTree(std::vector<rrt_vertex> &tree, int dim){
    std::vector<Eigen::Vector3d> projection;
    for (rrt_vertex v : tree){
        Eigen::VectorXd config = removeTwist(v.config);
        projection.emplace_back(pca.project(config,dim));
    }
    return projection;
    
}
void TreeVisualizer::setTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree){
    std::vector<std::array<size_t,2>> start_edges = getEdges(start_tree);
    std::vector<Eigen::Vector3d> start_projection = projectTree(start_tree);
    start_tree_curve = polyscope::registerCurveNetwork("start tree", start_projection, start_edges);
    start_tree_curve->setRadius(0.0001, true); 
    start_tree_curve->resetTransform();
}
void TreeVisualizer::updateTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree){
    return;
}