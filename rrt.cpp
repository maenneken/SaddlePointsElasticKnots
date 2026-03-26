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
Eigen::VectorXd reversePosIndicies(const Eigen::VectorXd& dofs){
    size_t n_pts = dofs.size()/4;
    Eigen::VectorXd reversed = Eigen::VectorXd::Zero(dofs.size());
    for ( size_t i = 0; i < n_pts; ++i){
        size_t j = n_pts - 1 - i;
        reversed.segment<3>(3*i) = dofs.segment<3>(3*j);
    }
    return reversed;
}
Eigen::VectorXd shiftPosIndcies(const Eigen::VectorXd& dofs){
    size_t n_pts = dofs.size()/4;
    Eigen::VectorXd shifted = Eigen::VectorXd::Zero(dofs.size());
    for ( size_t i = 0; i < n_pts; ++i){
        size_t i1 = (i+1) % n_pts;
        shifted.segment<3>(3*i) = dofs.segment<3>(3*i1);
    }
    return shifted;
}
void RRT::add_all_permutations(const Eigen::VectorXd& start, const Eigen::VectorXd& goal){
    size_t n_pts = start.size()/4;
    for(size_t i = 0; i < n_pts -1; i+= every_k_permutation){
        rrt_vertex current = start_tree.back();
        rrt_vertex start_shifted(shiftPosIndcies(current.config),-1);
        start_shifted.projection = pca.project(current.config,proj_dim);   
        start_tree.emplace_back(start_shifted);

        current = goal_tree.back();
        rrt_vertex goal_shifted(shiftPosIndcies(current.config),-1);
        goal_shifted.projection = pca.project(current.config,proj_dim);
        goal_tree.emplace_back(goal_shifted);

        //add weights so they can be selected
        start_weight.emplace_back(1.0);
        goal_weight.emplace_back(1.0);
    }
    //now all backwards
    rrt_vertex start_backwards (reversePosIndicies(start),-1);
    start_backwards.projection = pca.project(start_backwards.config,proj_dim);
    rrt_vertex goal_backwards (reversePosIndicies(goal),-1);
    goal_backwards.projection = pca.project(goal_backwards.config,proj_dim);
    start_tree.emplace_back(start_backwards);
    goal_tree.emplace_back(goal_backwards);

    for(size_t i = 0; i < n_pts -1; i+= every_k_permutation){
        rrt_vertex current = start_tree.back();
        rrt_vertex start_shifted(shiftPosIndcies(current.config),-1);
        start_shifted.projection = pca.project(current.config,proj_dim);
        start_tree.emplace_back(start_shifted);

        current = goal_tree.back();
        rrt_vertex goal_shifted(shiftPosIndcies(current.config),-1);
        goal_shifted.projection = pca.project(current.config,proj_dim);
        goal_tree.emplace_back(goal_shifted);

        //add weights so they can be selected
        start_weight.emplace_back(1.0);
        goal_weight.emplace_back(1.0);
    }
}

RRT::RRT(const Eigen::VectorXd& start,
         const Eigen::VectorXd& goal,
         PCA _pca, bool all_permutations,
         size_t _every_k_permutation
         ) {

    assert(start.size() == goal.size() && "start and goal must have same dimension");

    pca = _pca;

    rrt_vertex start_vertex(start, -1);
    start_vertex.projection = pca.project(start,proj_dim);
    rrt_vertex goal_vertex(goal, -1);
    goal_vertex.projection = pca.project(goal,proj_dim);
    

    start_tree.emplace_back(start_vertex);
    goal_tree.emplace_back(goal_vertex);
    start_weight.emplace_back(1.0);
    goal_weight.emplace_back(1.0);

    every_k_permutation = _every_k_permutation;

    //add all permutation of start and goal to the trees
    if(all_permutations){
        add_all_permutations(start,goal);
    }


    min_val = std::min(start.minCoeff(), goal.minCoeff()) - 10;
    max_val = std::max(start.maxCoeff(), goal.maxCoeff()) + 10;

    n_dofs = start.size();
    n_permutations = start_tree.size();
}
//prunes all leafs
//weights need to be recalulated
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
//TODO add goal Bias where we want to expand towards goal but use the knot that is most "effective"
Eigen::VectorXd RRT::sampleRandDirection(ContactProblem& cp,const Eigen::VectorXd& current_config, const Eigen::VectorXd& goal){
    
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

        //sample in PCA space
        } else if (sample_in_projection_space){
            auto projection = pca.project(current_config,sample_proj_dim);
            Eigen::VectorXd rand_proj = Eigen::VectorXd::Random(sample_proj_dim);
            projection += rand_proj;
            randDirection = pca.reconstruct(projection,sample_proj_dim) - current_config_short;
        }
        else {
            randDirection = Eigen::VectorXd::Random(n_pts);
        }
    }

    
    if(use_constraint_projection_for_sampling){
        //project to constraint space to preserve constraints
        //stretch and bend seem to be enough 
        Eigen::SparseMatrix<double> J = stackJacobians(stretchJacobian(current_config_short),bendJacobian(current_config_short));

        //force rotation
        if (uniform01() < rotation_bias){
            J = J = stackJacobians(J,twistJacobian(current_config_short));
        }
    

        Eigen::VectorXd projectedDirection = projectToTangentSpace(J,randDirection);
        randDirection = projectedDirection;
    }
    //proj to pca space and back to remove noise
    if(reproject_direction){
        Eigen::VectorXd  proj_step = current_config_short + randDirection.normalized() * step_length;
        Eigen::VectorXd  reprojected_proj_step = pca.reconstruct(pca.project(proj_step,sample_proj_dim),sample_proj_dim);
       randDirection = reprojected_proj_step - current_config_short;
    }

    Eigen::VectorXd fullRandomDirection = Eigen::VectorXd::Zero(n_dofs);
    fullRandomDirection.head(n_pts) = randDirection;


    return fullRandomDirection;    
}

void RRT::buildKDTrees(){
    std::cout << "RRT::buildKDTrees() called" << std::endl;
    std::cout << "start_tree.size() = " << start_tree.size() 
              << ", goal_tree.size() = " << goal_tree.size() << std::endl;
    std::cout << "proj_dim = " << proj_dim << std::endl;

    start_adaptor = std::make_unique<RRTAdaptor>(start_tree);
    goal_adaptor  = std::make_unique<RRTAdaptor>(goal_tree);

    std::cout << "Adaptors created" << std::endl;

    if(start_tree.empty() || goal_tree.empty()){
        std::cerr << "Error: one of the trees is empty!" << std::endl;
        return;
    }

    start_kd_tree = std::make_unique<KDTree>(
        proj_dim,
        *start_adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(10)
    );
    goal_kd_tree = std::make_unique<KDTree>(
        proj_dim,
        *goal_adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(10)
    );

    std::cout << "KDTree objects created" << std::endl;

    start_kd_tree->buildIndex();
    goal_kd_tree->buildIndex();

    std::cout << "KD-trees built successfully" << std::endl;
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
std::vector<size_t> RRT::kNearestNeighbor(const Eigen::VectorXd& query_projection, const KDTree& index, int k){

    assert((query_projection.size() == proj_dim ) && "knn wrong dim");

    std::vector<size_t> ret_indexes(k);
    std::vector<double> out_dists_sqr(k);

    index.knnSearch(
        query_projection.data(),
        k,
        ret_indexes.data(),
        out_dists_sqr.data()
    );

    return ret_indexes;
}
std::vector<size_t> RRT::radiusNeighbor(const Eigen::VectorXd& query_projection, const KDTree& index, double r){

    assert((query_projection.size() == proj_dim ) && "rnn wrong dim");

    std::vector<nanoflann::ResultItem<size_t,double>> matches;

    nanoflann::SearchParameters params;

    const double r2 = r * r;

    index.radiusSearch(
        query_projection.data(),
        r2,
        matches,
        params
    );

    std::vector<size_t> result;
    result.reserve(matches.size());

    for(const auto& m : matches)
        result.push_back(m.first);

    return result;
}

//steer towards the rand configuration linear, till rand is reached, energy is to high, or max steer length is reached 
Eigen::VectorXd RRT::steerTowardsConfig(ContactProblem& cp, Eigen::VectorXd& near, Eigen::VectorXd& rand){
    Eigen::VectorXd direction = rand - near;
    
    if (direction.squaredNorm() < 1e-8)
        return near;  // nothing to do
    direction.normalize();
    Eigen::VectorXd current = near;
    cp.setVars(current);
    //go tillnot possible anymore; removed && (current - near).norm() <= step_length
    while(cp.contactEnergy() <= max_energy  && ((current - rand).squaredNorm() >= (steer_step*steer_step))){
        current += steer_step * direction;
        cp.setVars(current);
    }
    if((current - near).squaredNorm() < 1e-8){
        return near;
    }
    //we reached rand position
    if((current - rand).squaredNorm() <= (steer_step*steer_step)){
        return rand;
    }
    //we went one to far
    current -= steer_step * direction;
    return current;
}

//steer in direction linear, till energy is to high, or max steer length is reached 
Eigen::VectorXd RRT::steerInDirection(ContactProblem& cp, Eigen::VectorXd& config, Eigen::VectorXd& direction){
    
    if (direction.squaredNorm() < 1e-8)
        return config;  // nothing to do
    direction.normalize();


    Eigen::VectorXd current = config;

    cp.setVars(current);

    if(uniform01() < gradient_bias){
        for (size_t i = 0; i <= 50; i++){
            cp.setVars(current);
            current -= 0.0001 * cp.gradient(true);
        }
        return current;
    }
    //go till not possible anymore or step_length is reached
    while(cp.contactEnergy() <= max_energy && (current - config).squaredNorm() <= (step_length*step_length)){
        current += steer_step * direction;
        cp.setVars(current);
    }
    if((current - config).squaredNorm() < 1e-8){
        return config;
    }
    //we went one to far
    current -= steer_step * direction;
    return current;
}

std::vector<Eigen::VectorXd> RRT::createPath(Eigen::VectorXd connection){
    //the connection is always the last knot in the tree list
    int current = start_tree.size()-1;
    std::vector<Eigen::VectorXd> start_path;
    //go till you find the start
    while( start_tree[current].parent >=0 ){
        std::cout << start_tree[current].parent<< std::endl;
        start_path.emplace_back(start_tree[current].config);
        current = start_tree[current].parent;
    }
    //add start
    start_path.emplace_back(start_tree[current].config);
    //the path is backwards
    std::reverse(start_path.begin(), start_path.end());
    
    //the connection is always the last knot in the tree list
    current = goal_tree.size()-1;
    std::vector<Eigen::VectorXd> goal_path;
    //go till you find goal
    while( goal_tree[current].parent >=0 ){
        goal_path.emplace_back(goal_tree[current].config);
        current = goal_tree[current].parent;
    }
    goal_path.emplace_back(goal_tree[current].config);
    
    //first two elements are the same
    if(!goal_path.empty()) goal_path.erase(goal_path.begin());
    if(!goal_path.empty()) goal_path.erase(goal_path.begin());

    std::vector<Eigen::VectorXd> full_path = start_path;
    full_path.insert(full_path.end(), goal_path.begin(), goal_path.end());

    return full_path;

}
std::vector<double> RRT::updateAllWeights( const std::vector<rrt_vertex>& tree, KDTree& kd_tree, double r){
    std::vector<double> weights;
    for (rrt_vertex v: tree){
        double weight = 1.0 / radiusNeighbor(v.projection, kd_tree,r).size();
        weights.emplace_back(weight);
    }
    return weights;
}

void RRT::updateNeighboringWeights(const Eigen::VectorXd& config,std::vector<rrt_vertex>& tree,std::vector<double>& weights, const KDTree& kd_tree, double r){
    std::vector<size_t> neighbors = radiusNeighbor(pca.project(config,proj_dim), kd_tree,r);
    for (size_t n: neighbors){
        double weight = 1.0 / (radiusNeighbor(tree[n].projection, kd_tree, r).size() + 1);
        weights[n] = weight;
    }
}

//todo remove viewer
void RRT::expand(ContactProblem& cp, size_t config_id, const Eigen::VectorXd& goal, std::vector<rrt_vertex>& tree,std::vector<double>& weights, KDTree& kd_tree, size_t k){
    Eigen::VectorXd config = tree[config_id].config;

    for(size_t i = 0; i < k; ++i){
        Eigen::VectorXd rand_direction = sampleRandDirection(cp, config, goal);
        
        Eigen::VectorXd desired = config + step_length * rand_direction.normalized();

        size_t nn = radiusNeighbor(pca.project(desired,proj_dim), kd_tree,radius).size();
        if(nn <=0){
            nn = 1;
        }
        double weight = 1.0 / nn;

        // if weight is low it is not that interesting. so skip
        if(uniform01() > weight){
            continue;
        } 
        
        //try to go in rand direction
        Eigen::VectorXd new_config  = steerInDirection(cp, config, rand_direction);

        //could not step towards it
        if((new_config - config).squaredNorm() < 1e-8) continue;

        //add new to the tree
        rrt_vertex new_vertex(new_config,config_id);
        new_vertex.projection = pca.project(new_config,proj_dim);
        tree.emplace_back(new_vertex);
        weights.emplace_back(1.0);

        // cp.setVars(new_config);
        // //relax if energy is to high
        // if(cp.energy() > max_rod_energy){

        //     auto optimizerOptions = NewtonOptimizerOptions();
        //     optimizerOptions.niter = 1000;
        //     optimizerOptions.gradTol =  1e-3;

        //     auto problemOptions = cp.m_options;

        //     compute_equilibrium(cp.m_rods,problemOptions,optimizerOptions);
        //     cp.updateCachedVars();
    
        //     Eigen::VectorXd relaxed =  cp.getVars();
        //     rrt_vertex relaxed_vertex(relaxed,tree.size()-1);
        //     relaxed_vertex.projection = pca.project(relaxed,proj_dim);
        //     tree.emplace_back(relaxed_vertex);
        //     weights.emplace_back(1.0);
        // }
    }
    updateNeighboringWeights(config, tree, weights,kd_tree, radius);
}
bool RRT::connect(ContactProblem& cp, double l){
    double l2 = l*l;
    for ( size_t i = start_tree.size() -1 ; i > start_tree.size() - 11; --i){
        rrt_vertex v_start = start_tree[i];
        
        auto nearest = kNearestNeighbor(v_start.projection, *goal_kd_tree, 1);
        for (size_t n : nearest){
            rrt_vertex v_goal = goal_tree[n];

            //if((v_goal.projection-v_start.projection).squaredNorm() > l2) continue;
                
            Eigen::VectorXd config = steerTowardsConfig(cp, v_start.config, v_goal.config);

            if((config - v_start.config).squaredNorm() < 1e-8) continue;


            //we can connect both trees
            if((config - v_goal.config).squaredNorm() < 1e-2){
                //add new to the tree
                std::cout <<"Start_Tree: trees connect" << std::endl;

                //add goal node to start Tree
                rrt_vertex start_new_vertex(v_goal.config,i);
                start_new_vertex.projection = v_goal.projection;
                start_tree.emplace_back(start_new_vertex);

                //add start node to goal Tree
                rrt_vertex goal_new_vertex(v_start.config,n);
                goal_new_vertex.projection = v_start.projection;
                goal_tree.emplace_back(goal_new_vertex);


                 return true;
            }            
        }
        
    }
    for ( size_t i = goal_tree.size() -1 ; i > goal_tree.size() - 11; --i){
        rrt_vertex v_goal = goal_tree[i];

        auto nearest = kNearestNeighbor(v_goal.projection, *start_kd_tree, 1);
        for (size_t n : nearest){
            rrt_vertex v_start = start_tree[n];

            //if((v_goal.projection-v_start.projection).squaredNorm() > l2) continue;
                
            Eigen::VectorXd config = steerTowardsConfig(cp, v_start.config, v_goal.config);

            if((config - v_start.config).squaredNorm() < 1e-8) continue;


            //we can connect both trees
            if((config - v_goal.config).squaredNorm() < 1e-2){
                //add new to the tree
                std::cout <<"Goal_Tree: trees connect" << std::endl;

                //add goal node to start Tree
                rrt_vertex start_new_vertex(v_goal.config,n);
                start_new_vertex.projection = v_goal.projection;
                start_tree.emplace_back(start_new_vertex);

                //add start node to goal Tree
                rrt_vertex goal_new_vertex(v_start.config,i);
                goal_new_vertex.projection = v_start.projection;
                goal_tree.emplace_back(goal_new_vertex);
                
                 return true;
            }            
        }
    }
    return false;
}

std::vector<Eigen::VectorXd> RRT::findConstrainedPath(ContactProblem& cp, size_t iterations, TreeVisualizer& Viewer){
    buildKDTrees();
    std::cout << "build kd trees"<< std::endl;
    std::vector<Eigen::VectorXd> path;
    std::array<std::vector<rrt_vertex>*, 2> trees = { &start_tree, &goal_tree };
    std::array<std::vector<double>*, 2> weights = {&start_weight,&goal_weight};
    std::array<KDTree*, 2> kd_trees = { start_kd_tree.get(), goal_kd_tree.get() };
    

    std::random_device rd;
    std::mt19937 gen(rd());
    std::cout << "start search" << std::endl;
    for(size_t i = 0; i < iterations; ++i){
        //switch between the two trees to expand them
        for(size_t t = 0; t < 2; ++ t){
            size_t goal_id = 0;
            if(goal_bias_for_all_permutations){
                size_t goal_id = static_cast<size_t>(uniform01() * n_permutations); //random goal Knot
            }
            
            Eigen::VectorXd goal = (*trees[1 - t])[goal_id].config;
            Eigen::VectorXd start = (*trees[t])[0].config;

            size_t config_id = 0;

            
            //goalbias = exploitation
            //pick knot nearest to goal
            if(uniform01() < goal_bias){
                auto nears = kNearestNeighbor(pca.project(goal,proj_dim),*kd_trees[t], k_neighbors);
                if(nears.empty()){
                    continue; // skip this expansion, KD-tree has no neighbors
                }
                size_t id = static_cast<size_t>(uniform01() * nears.size()); //random near Knot
                config_id = nears[id];
            } else {
                std::discrete_distribution<> dd(weights[t]->begin(), weights[t]->end());
                config_id = dd(gen);
            }
            //try to go towards the tree
            if(uniform01() < (goal_bias)){
                auto nears = kNearestNeighbor((*trees[t])[config_id].projection,*kd_trees[1-t], k_neighbors);
                if(nears.empty()){
                    continue; // skip this expansion, KD-tree has no neighbors
                }
                size_t id = static_cast<size_t>(uniform01() * nears.size()); //random near Knot
                goal = (*trees[1-t])[nears[id]].config;
            }
            expand(cp,config_id, goal, *trees[t], *weights[t], *kd_trees[t], 10);
        
        } 
        if(i % 10 == 0){
            start_kd_tree->buildIndex();
            goal_kd_tree->buildIndex();
        }
        if(connect(cp, 2*step_length)){
            Viewer.setTree(start_tree,goal_tree); 
            return createPath(start_tree.back().config);
        }

        // //pruning the tree   
        // if(i>0 && i % pruning_interval == 0){
        //     std::cout << std::endl;
        //     std::cout << RED <<"Trees are getting pruned" << RESET << std::endl;
        //     std::cout << "Number of Vertcies befor pruning: "<< BLUE << "start_tree: "  << start_tree.size() << "; goal_tree: " << goal_tree.size() << RESET << std::endl;
        //     auto nearest_goal = nearestVertex(goal_tree[0].config, start_tree);
        //     start_tree = pruneAllLeafNodes(start_tree);
        //     goal_tree = pruneAllLeafNodes(goal_tree);
        //     std::cout << "Number of Vertcies after pruning: "<< GREEN << "start_tree: " << start_tree.size() << "; goal_tree: " << goal_tree.size() << RESET <<std::endl;
        //     start_kd_tree->buildIndex();
        //     goal_kd_tree->buildIndex();
        //     start_weight = updateAllWeights(start_tree,*start_kd_tree,radius);
        //     goal_weight = updateAllWeights(goal_tree,*goal_kd_tree,radius);

        // }
        
        if(i % 10 == 0){
            std::cout << std::endl;
            std::cout << "Iteration: " << i << "; Number of Vertcies in start_tree: " << start_tree.size() << "; Number of Vertcies in goal_tree: " << goal_tree.size() << std::endl;
            auto nearest_goal = kNearestNeighbor(goal_tree[0].projection,*start_kd_tree, 1)[0];
            cp.setVars(start_tree[nearest_goal].config);
            std::cout << "nearest Vertex to goal is: " << nearest_goal << " with distance " << (start_tree[nearest_goal].config - goal_tree[0].config).norm() << " and Energy: " <<cp.energy()<<" Contact Energy: "<< cp.contactEnergy() << std::endl;
            auto nearest_start = kNearestNeighbor(start_tree[0].projection,*goal_kd_tree, 1)[0];
            // std::cout << "nearest Vertex to start is: " << nearest_start << " with distance " << (goal_tree[nearest_start].config - start_tree[0].config).norm() <<  std::endl;
            // std::cout << "highest weight in start tree: " << *std::max_element(start_weight.begin(), start_weight.end()) << std::endl;
            // std::cout << "lowest weight in start tree: " << *std::min_element(start_weight.begin(), start_weight.end()) << std::endl;
            // std::cout << "highest weight in goal tree: " << *std::max_element(goal_weight.begin(), goal_weight.end()) << std::endl;
            // std::cout << "lowest weight in goal tree: " << *std::min_element(goal_weight.begin(), goal_weight.end()) << std::endl;
            // auto pca_pts = pca.reconstruct(start_tree[nearest_goal].projection,proj_dim);
            // auto pts = DoFsToPos(pca_pts,0.25*n_dofs);
            auto pts = DoFsToPos(start_tree[nearest_goal].config,0.25*n_dofs);
            
            Viewer.updateKnot(pts); 
            auto direction_to_goal = DoFsToPos(goal_tree[0].config - start_tree[nearest_goal].config,0.25*n_dofs);
            Viewer.showNodeGradient(direction_to_goal);
            Viewer.setTree(start_tree,goal_tree); 
            
        } 
        Viewer.frameTick();  
        
    }
    return path;
}