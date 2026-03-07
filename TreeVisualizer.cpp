#include "TreeVisualizer.h"

TreeVisualizer::TreeVisualizer(const std::string file){
    std::vector<Eigen::VectorXd> data = loadPathTxt(file);
    std::vector<Eigen::VectorXd> samples = removeTwist(data);
    pca.fit(samples);
}
std::vector<std::array<size_t,2>> TreeVisualizer::getEdges(std::vector<rrt_vertex> &tree){
    std::vector<std::array<size_t,2>> edges;
    for (size_t i = 0; i < tree.size(); ++i){
        //it has a parent
        if(tree[i].parent > 0){
            edges.push_back({tree[i].parent, i});
        }
    }
    return edges;
}
std::vector<Eigen::Vector3d> TreeVisualizer::projectTree(std::vector<rrt_vertex> &tree, int dim){
    std::vector<Eigen::Vector3d> projection;
    for (rrt_vertex v : tree){
        Eigen::VectorXd config = removeTwist({v.config})[0];
        projection.emplace_back(pca.project(config,dim));
    }
    return projection;
    
}
void TreeVisualizer::setTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree){
    std::vector<std::array<size_t,2>> start_edges = getEdges(start_tree);
    std::vector<Eigen::Vector3d> start_projection = projectTree(start_tree,3);
    start_tree_curve = polyscope::registerCurveNetwork("start tree", start_projection, start_edges);
    start_tree_curve->setRadius(0.0005, true); 
    start_tree_curve->resetTransform();

    std::vector<std::array<size_t,2>> goal_edges = getEdges(goal_tree);
    std::vector<Eigen::Vector3d> goal_projection = projectTree(goal_tree,3);
    goal_tree_curve = polyscope::registerCurveNetwork("goal tree", goal_projection, goal_edges);
    goal_tree_curve->setRadius(0.0005, true); 
    goal_tree_curve->resetTransform();
}
void TreeVisualizer::updateTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree){
    return;
}