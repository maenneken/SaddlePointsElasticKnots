#include "rrt.h"
//create a rrt tree to find the path between to states of the knot
//use max Energie as constrained. Knots of higher energie do not exists, they are in "collision"
//create tree from start and from goal knot
//use a step size
//connection between nodes is linear movement
//try different sampling methods

//sample random Knot
//look for nerarest knot in tree
//go towards sample from it till max Energy, sample or step length is reached
//create new node and connect and add to tree
//try to connect new node with nearest node from other tree
//if possible end and look for path between start and goal
//else switch tree and sample again

double uniform01() {
    static std::mt19937 rng(std::random_device{}()); // random engine
    static std::uniform_real_distribution<double> dist(0.0, 1.0); // 0 ≤ x < 1
    return dist(rng);
}
Eigen::VectorXd gaussianVector(size_t n, double sigma) {
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> dist(0.0, sigma);

    Eigen::VectorXd v(n);
    for (size_t i = 0; i < n; ++i)
        v[i] = dist(rng);

    return v;
}


RRT::RRT(const Eigen::VectorXd& start,
         const Eigen::VectorXd& goal,
         double maxEnergy,
         double stepLength,
         double goalBias,
         double steerStep){
    assert(start.size() == goal.size() && "start and goal must have same dimension");
    assert(maxEnergy > 0 && "maxEnergy must be positive");
    assert(stepLength > 0 && "stepLength must be positive");
    assert(steerStep > 0 && "steerStep must be positive");

    rrt_vertex start_vertex(start, -1);
    rrt_vertex goal_vertex(goal, -1);
    start_tree.emplace_back(start_vertex);
    goal_tree.emplace_back(goal_vertex);

    max_energy = maxEnergy;
    step_length = stepLength;
    steer_step = steerStep;
    goal_bias = goalBias;

    min_val = std::min(start.minCoeff(), goal.minCoeff()) - 10;
    max_val = std::max(start.maxCoeff(), goal.maxCoeff()) + 10;

    n_dofs = start.size();
}
//sampels a rand configuration around existing knots in the tree or goal
//try sampling a random position for one random point in the vertex and have the neighboring vertecies go the same direction but less and less to simulate bending 
Eigen::VectorXd RRT::sampleRandConfig(const Eigen::VectorXd& goal, const std::vector<rrt_vertex>& tree){
    // goal bias
    if (uniform01() < goal_bias)
        return goal;

    // pick random node from the tree
    size_t idx = static_cast<size_t>(uniform01() * tree.size());
    const Eigen::VectorXd& center = tree[idx].config;

    double sigma = step_length;  // key parameter
    Eigen::VectorXd rand_config = center + gaussianVector(n_dofs, sigma);

    // keep twist fixed
    rand_config.tail(static_cast<int>(0.75 * n_dofs)).setZero();

    return rand_config;
}
size_t RRT::nearestVertex(const Eigen::VectorXd& config, const std::vector<rrt_vertex>& tree){
    size_t nearest = 0;
    double min_dist = std::numeric_limits<double>::infinity();
    for(size_t i = 0; i< tree.size(); ++i){
        double dist = (tree[i].config - config).squaredNorm();

        if (dist < min_dist) {
            min_dist = dist;
            nearest = i;
        }
    }
    return nearest;
}
//steer towards the rand configuration linear, till rand is reached, energy is to high, or max steer length is reached 
Eigen::VectorXd RRT::steerTowardsConfig(ContactProblem& cp, Eigen::VectorXd& near, Eigen::VectorXd& rand){
    Eigen::VectorXd direction = rand - near;
    
    if (direction.norm() < 1e-8)
        return near;  // nothing to do
    direction.normalize();
    Eigen::VectorXd current = near;
    cp.setVars(current);
    while(cp.energy() <= max_energy && (current - near).norm() <= step_length && ((current - rand).norm() >= 1e-8)){
        current += steer_step * direction;
        cp.setVars(current);
    }
    if((current - near).norm() < 1e-8){
        return near;
    }
    return current;
}
std::vector<Eigen::VectorXd> RRT::createPath(Eigen::VectorXd connection){
    auto current = nearestVertex(connection, start_tree);
    std::vector<Eigen::VectorXd> start_path;
    //go till you find the start
    while( start_tree[current].parent >=0 ){
        start_path.emplace_back(start_tree[current].config);
        current = start_tree[current].parent;
    }
    //the path is backwards
    std::reverse(start_path.begin(), start_path.end());
    
    current = nearestVertex(connection, goal_tree);
    std::vector<Eigen::VectorXd> goal_path;
    //go till you find goal
    while( goal_tree[current].parent >=0 ){
        goal_path.emplace_back(goal_tree[current].config);
        current = goal_tree[current].parent;
    }
    
    //remove duplicate
    if(!goal_path.empty()) goal_path.erase(goal_path.begin());

    std::vector<Eigen::VectorXd> full_path = start_path;
    full_path.insert(full_path.end(), goal_path.begin(), goal_path.end());

    return full_path;

}

std::vector<Eigen::VectorXd> RRT::findPath(ContactProblem& cp, size_t iterations, KnotVisualizer& Viewer){
    std::vector<Eigen::VectorXd> path;
    std::array<std::vector<rrt_vertex>*, 2> trees = { &start_tree, &goal_tree };

    for(size_t i = 0; i < iterations; ++i){
        //switch between the two trees to expand them
        for(size_t t = 0; t < 2; ++ t){
            
            Eigen::VectorXd goal = (*trees[1 - t])[0].config;
            
            //sample rand config
            Eigen::VectorXd rand_config = sampleRandConfig(goal,*trees[t]);

            //find the nearest vertex to the tree
            size_t nearest_id = nearestVertex(rand_config, *trees[t]);
            Eigen::VectorXd nearest = (*trees[t])[nearest_id].config;
           
            //try to gotowards rand_config
            Eigen::VectorXd new_config  = steerTowardsConfig(cp, nearest, rand_config);
            //could not step towards it
            if((new_config - nearest).norm() < 1e-8) continue;

            //add new to the tree
            rrt_vertex new_vertex(new_config,nearest_id);
            (*trees[t]).emplace_back(new_vertex);
            nearest = new_config;

                 

            //try to connect to the other tree
            nearest_id = nearestVertex(new_config, *trees[1-t]);
            nearest = (*trees[1-t])[nearest_id].config;
            Eigen::VectorXd other_config = steerTowardsConfig(cp, nearest, new_config);

            if((other_config - nearest).norm() < 1e-8)continue;

            //add new to the tree
            rrt_vertex other_new_vertex(other_config,nearest_id);
            (*trees[1-t]).emplace_back(other_new_vertex);
            nearest = other_config;
            
            //we can connect both trees
            if((other_config - new_config).norm() < 1e-8){
                return createPath(new_config);
            }
        }    
        
        
        if(i % 100 == 0){
            std::cout << std::endl;
            std::cout << "Iteration: " << i << "; Number of Vertcies in start_tree: " << start_tree.size() << "; Number of Vertcies in goal_tree: " << goal_tree.size() << std::endl;
            auto nearest_goal = nearestVertex(goal_tree[0].config, start_tree);
            cp.setVars(start_tree[nearest_goal].config);
            std::cout << "nearest Vertex to goal is: " << nearest_goal << " with distance " << (start_tree[nearest_goal].config - goal_tree[0].config).norm() << " and Energy: " <<cp.energy() << std::endl;
            auto nearest_start = nearestVertex(start_tree[0].config, goal_tree);
            std::cout << "nearest Vertex to start is: " << nearest_start << " with distance " << (goal_tree[nearest_start].config - start_tree[0].config).norm() <<  std::endl;

            auto pts = DoFsToPos(cp.getVars(),0.25*n_dofs);
            Viewer.updateKnot(pts);
            Viewer.frameTick();
        }   
        
    }
    return path;
}