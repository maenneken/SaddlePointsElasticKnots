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
         double steerStep,
         size_t pruningInterval,
         bool oneRandDirection){
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
    pruning_interval = pruningInterval;
    one_rand_direction_3d = oneRandDirection;

    min_val = std::min(start.minCoeff(), goal.minCoeff()) - 10;
    max_val = std::max(start.maxCoeff(), goal.maxCoeff()) + 10;

    n_dofs = start.size();
}
std::vector<rrt_vertex> RRT::pruneAllLeafNodes(std::vector<rrt_vertex>& tree){
    std::vector<rrt_vertex> pruned_tree;
    //root must be in the tree
    pruned_tree.emplace_back(tree[0]);
    //go through all knots in tree and count if they are a parent
    for( size_t i = 1; i < tree.size(); ++i){
        //look if i is a parent
        size_t parent_count = 0;
        for( size_t j = i; j < tree.size();++j){
            if(tree[j].parent == i){
                ++parent_count;
                //update parent id
                tree[j].parent = pruned_tree.size();
            }
        }
        if(parent_count > 0){
            pruned_tree.emplace_back(tree[i]);
        }
        
    }
    return pruned_tree;
    
}
//sampels a rand configuration around existing knots in the tree or goal
//try sampling a random position for one random point in the vertex and have the neighboring vertecies go the same direction but less and less to simulate bending 
Eigen::VectorXd RRT::sampleRandConfig(ContactProblem& cp, const Eigen::VectorXd& goal, const std::vector<rrt_vertex>& tree){

    size_t n_vertices = 0.25*n_dofs;
    size_t idx_tree; //id of Knot in the tree
    size_t idx_knot = static_cast<size_t>(uniform01() * n_vertices); //random Point in Knot
    //const Eigen::Vector3d center(tree[idx_tree].config[3*idx_knot], tree[idx_tree].config[3*idx_knot +1], tree[idx_tree].config[3*idx_knot+2]);


    //pick knot nearest to goal
    if(uniform01() < goal_bias){
        idx_tree = nearestVertex(goal,tree);
    } 
    else {
        idx_tree = static_cast<size_t>(uniform01() * tree.size()); //random Knot
    }

    Eigen::VectorXd randKnot = tree[idx_tree].config;

    Eigen::Vector3d randDirection;

    //goalBias
    if (uniform01() < goal_bias){
        randKnot.segment<3>(3*idx_knot) = goal.segment<3>(3*idx_knot);
    } 
    else {
        randDirection = Eigen::Vector3d::Random().normalized()*step_length;
        randKnot.segment<3>(3*idx_knot) += randDirection;
    }
    
    return randKnot;
}
//Samples a random direction and projects it to presever spring constraint
//if goal bias direction goes towards goal
Eigen::VectorXd RRT::sampleRandDirection(const Eigen::VectorXd& current_config, const Eigen::VectorXd& goal){
    
    size_t n_pts = 0.75*n_dofs;
    Eigen::VectorXd current_config_short = current_config.head(n_pts);
    Eigen::VectorXd goal_short = goal.head(n_pts);
    Eigen::VectorXd randDirection = Eigen::VectorXd::Zero(n_pts);;


    //move one vertex in 3d space
    if(one_rand_direction_3d){
        size_t idx_vertex = static_cast<size_t>(uniform01() * (n_pts/3)); //random Point in Knot
        //goalBias
        if (uniform01() < goal_bias){
        randDirection.segment<3>(3*idx_vertex) = goal.segment<3>(3*idx_vertex) - current_config.segment<3>(3*idx_vertex);
        } 
        else {
        Eigen::Vector3d rand3d = Eigen::Vector3d::Random().normalized();
        randDirection.segment<3>(3*idx_vertex) += rand3d;
        }
    }
    //move all vertecies in random directions
    else {
        //goalBias
        if (uniform01() < goal_bias){
            randDirection = goal_short - current_config_short;
        } 
        else {
            randDirection = Eigen::VectorXd::Random(n_pts);
        }
    }
    
    Eigen::SparseMatrix<double> J = stackJacobians(springJacobian(current_config_short),bendJacobian(current_config_short));
    Eigen::VectorXd projectedDirection = Eigen::VectorXd::Zero(n_dofs);
    projectedDirection.head(n_pts) = projectToTangentSpace(J,randDirection);

    return projectedDirection;

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
    //go tillnot possible anymore; removed && (current - near).norm() <= step_length
    while(cp.contactEnergy() <= max_energy  && ((current - rand).norm() >= steer_step)){
        current += steer_step * direction;
        cp.setVars(current);
    }
    if((current - near).norm() < 1e-8){
        return near;
    }
    //we reached rand position
    if((current - rand).norm() <= steer_step){
        return rand;
    }
    //we went one to far
    current -= steer_step * direction;
    return current;
}

//steer in direction linear, till energy is to high, or max steer length is reached 
Eigen::VectorXd RRT::steerInDirection(ContactProblem& cp, Eigen::VectorXd& config, Eigen::VectorXd& direction){
    
    if (direction.norm() < 1e-8)
        return config;  // nothing to do
    direction.normalize();
    Eigen::VectorXd current = config;

    cp.setVars(current);
    //go till not possible anymore or step_length is reached
    while(cp.contactEnergy() <= max_energy && (current - config).norm() <= step_length){
        current += steer_step * direction;
        cp.setVars(current);
    }
    if((current - config).norm() < 1e-8){
        return config;
    }
    //we went one to far
    current -= steer_step * direction;
    return current;
}

std::vector<Eigen::VectorXd> RRT::createPath(Eigen::VectorXd connection){
    auto current = nearestVertex(connection, start_tree);
    std::vector<Eigen::VectorXd> start_path;
    //go till you find the start
    while( start_tree[current].parent >=0 ){
        std::cout << start_tree[current].parent<< std::endl;
        start_path.emplace_back(start_tree[current].config);
        current = start_tree[current].parent;
    }
    //add start
    start_path.emplace_back(start_tree[0].config);
    //the path is backwards
    std::reverse(start_path.begin(), start_path.end());
    
    current = nearestVertex(connection, goal_tree);
    std::vector<Eigen::VectorXd> goal_path;
    //go till you find goal
    while( goal_tree[current].parent >=0 ){
        goal_path.emplace_back(goal_tree[current].config);
        current = goal_tree[current].parent;
    }
    goal_path.emplace_back(goal_tree[0].config);
    
    //remove duplicate
    if(!goal_path.empty()) goal_path.erase(goal_path.begin());

    std::vector<Eigen::VectorXd> full_path = start_path;
    full_path.insert(full_path.end(), goal_path.begin(), goal_path.end());

    return full_path;

}

std::vector<Eigen::VectorXd> RRT::findConstrainedPath(ContactProblem& cp, size_t iterations, KnotVisualizer& Viewer){
    std::vector<Eigen::VectorXd> path;
    std::array<std::vector<rrt_vertex>*, 2> trees = { &start_tree, &goal_tree };

    for(size_t i = 0; i < iterations; ++i){
        //switch between the two trees to expand them
        for(size_t t = 0; t < 2; ++ t){
            
            Eigen::VectorXd goal = (*trees[1 - t])[0].config;

            size_t config_id;
            //pick knot nearest to goal
            if(uniform01() < goal_bias){
                config_id = nearestVertex(goal,*trees[t]);
            } else {
                config_id = static_cast<size_t>(uniform01() * (*trees[t]).size()); //random Knot
            }
            Eigen::VectorXd config = (*trees[t])[config_id].config;
            Eigen::VectorXd rand_direction = sampleRandDirection(config, goal);
           
            //try to go in rand direction
            Eigen::VectorXd new_config  = steerInDirection(cp, config, rand_direction);


            //could not step towards it
            if((new_config - config).norm() < 1e-8) continue;

            //add new to the tree
            rrt_vertex new_vertex(new_config,config_id);
            (*trees[t]).emplace_back(new_vertex);

            Eigen::VectorXd nearest = new_config;

            //try to connect to the other tree
            size_t nearest_id = nearestVertex(new_config, *trees[1-t]);
            nearest = (*trees[1-t])[nearest_id].config;
            Eigen::VectorXd other_config = steerTowardsConfig(cp, nearest, new_config);

            if((other_config - nearest).norm() < 1e-8)continue;

            
            //we can connect both trees
            if((other_config - new_config).norm() < 1e-2){
                //add new to the tree
                std::cout <<"trees connect" << std::endl;
                rrt_vertex other_new_vertex(other_config,nearest_id);
                (*trees[1-t]).emplace_back(other_new_vertex);
                nearest = other_config;
                return createPath(new_config);
            }
        } 
        //pruning the tree   
        if(i % pruning_interval == 0){
            std::cout << std::endl;
            std::cout << RED <<"Trees are getting pruned" << RESET << std::endl;
            std::cout << "Number of Vertcies befor pruning: "<< BLUE << "start_tree: "  << start_tree.size() << "; goal_tree: " << goal_tree.size() << RESET << std::endl;
            auto nearest_goal = nearestVertex(goal_tree[0].config, start_tree);
            start_tree = pruneAllLeafNodes(start_tree);
            goal_tree = pruneAllLeafNodes(goal_tree);
            std::cout << "Number of Vertcies after pruning: "<< GREEN << "start_tree: " << start_tree.size() << "; goal_tree: " << goal_tree.size() << RESET <<std::endl;

        }
        
        if(i % 10 == 0){
            std::cout << std::endl;
            std::cout << "Iteration: " << i << "; Number of Vertcies in start_tree: " << start_tree.size() << "; Number of Vertcies in goal_tree: " << goal_tree.size() << std::endl;
            auto nearest_goal = nearestVertex(goal_tree[0].config, start_tree);
            cp.setVars(start_tree[nearest_goal].config);
            std::cout << "nearest Vertex to goal is: " << nearest_goal << " with distance " << (start_tree[nearest_goal].config - goal_tree[0].config).norm() << " and Energy: " <<cp.energy()<<" Contact Energy: "<< cp.contactEnergy() << std::endl;
            auto nearest_start = nearestVertex(start_tree[0].config, goal_tree);
            std::cout << "nearest Vertex to start is: " << nearest_start << " with distance " << (goal_tree[nearest_start].config - start_tree[0].config).norm() <<  std::endl;

            auto pts = DoFsToPos(start_tree[nearest_goal].config,0.25*n_dofs);
            Viewer.updateKnot(pts); 
            auto direction_to_goal = DoFsToPos(goal_tree[0].config - start_tree[nearest_goal].config,0.25*n_dofs);
            Viewer.showNodeGradient(direction_to_goal);
            
        } 
        Viewer.frameTick();  
        
    }
    return path;
}
std::vector<Eigen::VectorXd> RRT::findPath(ContactProblem& cp, size_t iterations, KnotVisualizer& Viewer){
    std::vector<Eigen::VectorXd> path;
    std::array<std::vector<rrt_vertex>*, 2> trees = { &start_tree, &goal_tree };

    for(size_t i = 0; i < iterations; ++i){
        //switch between the two trees to expand them
        for(size_t t = 0; t < 2; ++ t){
            
            Eigen::VectorXd goal = (*trees[1 - t])[0].config;
            
            //sample rand config
            Eigen::VectorXd rand_config = sampleRandConfig(cp, goal,*trees[t]);

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
            if((other_config - new_config).norm() < 1e-2){
                return createPath(new_config);
            }
        } 
        //pruning the tree   
        if(i % pruning_interval == 0){
            std::cout << std::endl;
            std::cout << RED <<"Trees are getting pruned" << RESET << std::endl;
            std::cout << "Number of Vertcies befor pruning: "<< BLUE << "start_tree: "  << start_tree.size() << "; goal_tree: " << goal_tree.size() << RESET << std::endl;
            auto nearest_goal = nearestVertex(goal_tree[0].config, start_tree);
            start_tree = pruneAllLeafNodes(start_tree);
            goal_tree = pruneAllLeafNodes(goal_tree);
            std::cout << "Number of Vertcies after pruning: "<< GREEN << "start_tree: " << start_tree.size() << "; goal_tree: " << goal_tree.size() << RESET <<std::endl;

        }
        
        if(i % 10 == 0){
            std::cout << std::endl;
            std::cout << "Iteration: " << i << "; Number of Vertcies in start_tree: " << start_tree.size() << "; Number of Vertcies in goal_tree: " << goal_tree.size() << std::endl;
            auto nearest_goal = nearestVertex(goal_tree[0].config, start_tree);
            cp.setVars(start_tree[nearest_goal].config);
            std::cout << "nearest Vertex to goal is: " << nearest_goal << " with distance " << (start_tree[nearest_goal].config - goal_tree[0].config).norm() << " and Energy: " <<cp.energy()<<" Contact Energy: "<< cp.contactEnergy() << std::endl;
            auto nearest_start = nearestVertex(start_tree[0].config, goal_tree);
            std::cout << "nearest Vertex to start is: " << nearest_start << " with distance " << (goal_tree[nearest_start].config - start_tree[0].config).norm() <<  std::endl;

            auto pts = DoFsToPos(cp.getVars(),0.25*n_dofs);
            Viewer.updateKnot(pts); 
        } 
        Viewer.frameTick();  
        
    }
    return path;
}