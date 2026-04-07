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
        if(tree[i].parent >= 0){
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

    start_edges = getEdges(start_tree);
    start_projection = projectTree(start_tree,3);
    start_tree_curve = polyscope::registerCurveNetwork("start tree", start_projection, start_edges);
    start_tree_curve->setRadius(0.0005, true); 
    start_tree_curve->resetTransform();

    goal_edges = getEdges(goal_tree);
    goal_projection = projectTree(goal_tree,3);
    goal_tree_curve = polyscope::registerCurveNetwork("goal tree", goal_projection, goal_edges);
    goal_tree_curve->setRadius(0.0005, true); 
    goal_tree_curve->resetTransform();

    last_start_tree_size = start_tree.size();
    last_goal_tree_size = goal_tree.size();
}
void TreeVisualizer::updateTree(std::vector<rrt_vertex>& start_tree, std::vector<rrt_vertex>& goal_tree){
    // Update start tree incrementally
    if(!update_tree_visualization) return;
    size_t current_start_size = start_tree.size();
    for (size_t i = last_start_tree_size; i < current_start_size; ++i) {
        if (start_tree[i].parent >= 0) {
            start_edges.push_back({static_cast<size_t>(start_tree[i].parent), i});
        }
        Eigen::VectorXd config = removeTwist({start_tree[i].config})[0];
        start_projection.emplace_back(pca.project(config, 3));
    }
    if (current_start_size > last_start_tree_size) {
        start_tree_curve = polyscope::registerCurveNetwork("start tree", start_projection, start_edges);
        start_tree_curve->setRadius(0.0005, true);
        start_tree_curve->resetTransform();
    }
    last_start_tree_size = current_start_size;

    // Update goal tree incrementally
    size_t current_goal_size = goal_tree.size();
    for (size_t i = last_goal_tree_size; i < current_goal_size; ++i) {
        if (goal_tree[i].parent >= 0) {
            goal_edges.push_back({static_cast<size_t>(goal_tree[i].parent), i});
        }
        Eigen::VectorXd config = removeTwist({goal_tree[i].config})[0];
        goal_projection.emplace_back(pca.project(config, 3));
    }
    if (current_goal_size > last_goal_tree_size) {
        goal_tree_curve = polyscope::registerCurveNetwork("goal tree", goal_projection, goal_edges);
        goal_tree_curve->setRadius(0.0005, true);
        goal_tree_curve->resetTransform();
    }
    last_goal_tree_size = current_goal_size;
}