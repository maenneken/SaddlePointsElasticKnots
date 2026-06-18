#include "rrt.h"
#include <algorithm>
#include <thread>   // for std::this_thread::sleep_for
#include <chrono>   // for std::chrono::milliseconds

void reflect_knot(Eigen::VectorXd& dofs, int axs){
    int n_pts = dofs.size() / 4 ;
    for(int i = 0; i < n_pts; ++i){
        dofs(3 * i + axs) *= -1;
    }
}
int main(int argc, char** argv) {
    //std::string start_file = "../data/L400-r0.2-UpTo9Crossings/4_1/0001.obj";
    //std::string goal_file = "../data/L400-r0.2-UpTo9Crossings/4_1/0033.obj";
    std::string start_file = "../data/NoCollision/reduced0001.obj";
    std::string goal_file = "../data/NoCollision/reduced0033.obj";
    double rod_radius = 0.2;
    int reductionFactor = 4;
    bool hasCollisions = true;
    int contactStiffness = 10000;


    // Parse command line arguments
    if (argc >= 3) {
        start_file = argv[1];
        goal_file = argv[2];
    }
    if (argc >= 4) {
        reductionFactor = std::stod(argv[3]);
    }
    if (reductionFactor < 1){
        reductionFactor = 1;
    }


    std::vector<double> params = {rod_radius, rod_radius};

    //Create material
    RodMaterial material(
        "ellipse",
        2000,     // Young's modulus
        0.3,      // Poisson's ratio
        params
    );

    // Read centerline nodes
    std::vector<Eigen::Vector3d> start_centerline = reduce_knot_resolution(read_nodes_from_file(start_file), reductionFactor);
    std::vector<Eigen::Vector3d> goal_centerline = reduce_knot_resolution(read_nodes_from_file(goal_file), reductionFactor);
    for (size_t i = 0; i < start_centerline.size(); ++i){
        start_centerline[i][0] *= 1;
        start_centerline[i][1] *= 1;
        start_centerline[i][2] *= 1;
    }
   

    int n_pts = start_centerline.size();  

    PeriodicRod start_pr = define_periodic_rod(start_centerline,material);
    PeriodicRod goal_pr = define_periodic_rod(goal_centerline,material);

    Eigen::VectorXd start_dofs = start_pr.getDoFs();
    Eigen::VectorXd goal_dofs = goal_pr.getDoFs();

    



    if(start_dofs.size() != goal_dofs.size()){
        throw std::runtime_error("start and goal are not the same size");
    }

    std::cout << "created pr"<< std::endl;
    
    // 4. Wrap into PeriodicRodList
    PeriodicRodList rod_list = PeriodicRodList(start_pr);

    std::cout << "created rodlist"<< std::endl;
    // 5. Setup problem options
    ContactProblemOptions problemOptions;
    problemOptions.hasCollisions = hasCollisions;
    problemOptions.contactStiffness = contactStiffness;
    problemOptions.dHat = 2 * rod_radius;

    std::cout << "set pr options"<< std::endl;
    
    // 6. Create contact problem
    ContactProblem cp(rod_list, problemOptions);

    
    

    double start_energy = cp.contactEnergy();
    std::cout << "Start contact Energy: " << cp.contactEnergy() << "Rod Energy: " << cp.energy() << std::endl;
    cp.setVars(goal_dofs);
    double goal_energy = cp.contactEnergy();
    std::cout << "Goal contact Energy: " << cp.contactEnergy()<< "Rod Energy: " << cp.energy() << std::endl;

    double edgeLength = (start_dofs.segment<3>(0) - start_dofs.segment<3>(3)).norm();
    std::cout << "Start Knot Edge length: " << edgeLength <<std::endl;

    edgeLength = (goal_dofs.segment<3>(0) - goal_dofs.segment<3>(3)).norm();
    std::cout << "Goal Knot Edge length: " << edgeLength <<std::endl;



    std::string pcaFile;
    int multiplier = 1;

    if(n_pts == 25){
        pcaFile = "../data/PCA/25V_Dataset.txt";
    }
    else if(n_pts == 50){
        pcaFile = "../data/PCA/50V_Dataset.txt";
        multiplier = 2;
    }
    else if(n_pts == 100){
        pcaFile = "../data/PCA/100V_Dataset.txt";
        multiplier = 4;
    }
    else{
        throw std::runtime_error("No PCA file for this number of vertices");
    }
    

    TreeVisualizer Viewer = TreeVisualizer(pcaFile);
    Viewer.setKnot(DoFsToPos(start_dofs,n_pts),0.01 * rod_radius);

    //show the goal state
    Viewer.setGoalKnot(DoFsToPos(goal_dofs,n_pts),0.01 * rod_radius);
    //Save first Knot Vars
    auto startKnot = cp.getVars();
     
    static int iterations = 5000;
    static int banchmark_runs = 10;
    static int pruningInterval = 1000;
    static double maxEnergy = std::max({start_energy,goal_energy})*2 +1;
    static double stepsize = 1 * multiplier;
    static double steplength = 40 * multiplier;
    static double goalBias = 0.2;
    static double near_goal_bias = 0.0;
    static double tree_bias = 0.2;
    static bool oneRandDirection = false;
    static bool allPermutations = true;
    static double neighbor_radius = 10;
    static bool goal_bias_for_all_permutations = false;
    static bool sample_in_projection_space = false;
    static bool use_constraint_projection_for_sampling = true;
    static bool reproject_direction = false;
    static int sample_proj_dim = 30;
    static int every_k_permutation = 1;
    static double max_rod_energy = 30000;
    static int k_neighbors = 10;
    static bool start_with_rotation = false;
    static double rotation_bias = 0.1;
    static double gradient_bias = 0.1;
    static bool pca_data_whitening = false;
    static int k_expension = 1;
    static bool update_tree_visualization = true;
    static int max_depth = 15;
    
   
    bool running = false;
    bool show_path = false;
    bool benchmarking = false;
    std::vector<Eigen::VectorXd> path;
    //set buttons
    Viewer.setUserCallback([&]() {
        //todo add controls to load a Knot and set up a contactproblem with all options
        ImGui::Begin("Controls");
        ImGui::InputInt("Iterations", &iterations,1000,10000);
        ImGui::InputInt("Benchmark runs", &banchmark_runs,1,100);
        //ImGui::InputInt("Pruning Interval", &pruningInterval,1000,10000);
        ImGui::InputDouble("Max Energy", &maxEnergy,(0.001),(1),"%.2f");
        ImGui::InputDouble("stepsize", &stepsize,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("stepLength", &steplength,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("goal Bias", &goalBias,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("near goal bias", &near_goal_bias,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("tree bias", &tree_bias,(0.001),(0.01),"%.4f");
        ImGui::InputDouble("rotation bias", &rotation_bias,(0.001),(1),"%.4f");
        ImGui::InputDouble("gradient bias", &gradient_bias,(0.001),(1),"%.4f");
        ImGui::InputDouble("neighbor radius", &neighbor_radius,(0.1),(0.1),"%.4f");
        ImGui::InputInt("sample projection dimension", &sample_proj_dim,1,10);
        ImGui::InputInt("every k permutation", &every_k_permutation,1,10);
        //ImGui::InputDouble("max rod energy", &max_rod_energy,(0.001),(1),"%.2f");
        ImGui::InputInt("k neighbors for goal bias", &k_neighbors,1,100);
        ImGui::InputInt("k expansion", &k_expension,1,10);
        ImGui::InputInt("max depth for weighting", &max_depth,1,100);


        ImGui::Checkbox("oneRandDirection", &oneRandDirection);
        ImGui::Checkbox("all Permutations of start and goal", &allPermutations);
        ImGui::Checkbox("goal bias for all permutations", &goal_bias_for_all_permutations);
        ImGui::Checkbox("sample in projection space", &sample_in_projection_space);
        ImGui::Checkbox("use constraint projection for sampling", &use_constraint_projection_for_sampling);
        ImGui::Checkbox("reproject direction", &reproject_direction);
        ImGui::Checkbox("start with rotation", &start_with_rotation);
        ImGui::Checkbox("PCA data whitening", &pca_data_whitening);

    
        

        if(ImGui::Button("Relax start and goal")){
            relax_start_goal(cp, start_dofs,goal_dofs);
            Viewer.setKnot(DoFsToPos(start_dofs,n_pts),0.01 * rod_radius);
            Viewer.setGoalKnot(DoFsToPos(goal_dofs,n_pts),0.01 * rod_radius);
        }
        if(ImGui::Button("reflect goal")){
            reflect_knot(goal_dofs, 2);
            Viewer.setGoalKnot(DoFsToPos(goal_dofs,n_pts),0.01 * rod_radius);
        }
        if (ImGui::Button("Find Path")) {
            running=true;
            
        }
        if (ImGui::Button("Benchmark")){ 
            running=true;
            benchmarking = true;     
        }

        if (ImGui::Button("Show Path")) {
            show_path=true;
        }
        if (ImGui::Button("Save Path")) {
            std::string header = "Stepsize: " + std::to_string(stepsize) + ", StepLength: " + std::to_string(steplength) + ", GoalBias: " + std::to_string(goalBias) + ", NeighborRadius: " + std::to_string(neighbor_radius) + ", SampleProjDim: " + std::to_string(sample_proj_dim) + ", OneRandDirection: " + std::to_string(oneRandDirection) + ", AllPermutations: " + std::to_string(allPermutations) + ", GoalBiasForAllPermutations: " + std::to_string(goal_bias_for_all_permutations) + ", SampleInProjectionSpace: " + std::to_string(sample_in_projection_space) + ", UseConstraintProjectionForSampling: " + std::to_string(use_constraint_projection_for_sampling) + ", ReprojectDirection: " + std::to_string(reproject_direction);
            savePathTxt("foundPath.txt",path, header);
        }
        if (ImGui::Button("toggle tree visualization")){
            update_tree_visualization = !update_tree_visualization;
            Viewer.update_tree_visualization = update_tree_visualization;
        }
  
        ImGui::End();   

    });
    //main loop
    while (!polyscope::windowRequestsClose()) { 
        Viewer.frameTick();
        if(running){
            Viewer.pca.whiten = pca_data_whitening;
            



            if (benchmarking) {
                std::vector<double> durations;
                std::vector<double> expension_quotas;
                std::vector<size_t> path_lengths;
                std::vector<size_t> tree_sizes;
                std::vector<size_t> iterations_quotas;
                size_t successful_runs = 0;

                for (int i = 0; i < banchmark_runs; ++i) {
                    RRT rrt(start_dofs, goal_dofs, Viewer.pca, allPermutations,every_k_permutation);
                    rrt.goal_bias = goalBias;
                    rrt.near_goal_bias = near_goal_bias;
                    rrt.tree_bias = tree_bias;
                    rrt.max_energy = maxEnergy;
                    rrt.pruning_interval = pruningInterval;
                    rrt.steer_step = stepsize;
                    rrt.step_length = steplength;
                    rrt.radius = neighbor_radius;
                    rrt.one_rand_direction_3d = oneRandDirection;
                    rrt.goal_bias_for_all_permutations = goal_bias_for_all_permutations;
                    rrt.sample_in_projection_space = sample_in_projection_space;
                    rrt.use_constraint_projection_for_sampling = use_constraint_projection_for_sampling;
                    rrt.reproject_direction = reproject_direction;
                    rrt.sample_proj_dim = sample_proj_dim;
                    rrt.max_rod_energy = max_rod_energy;
                    rrt.k_neighbors = k_neighbors;
                    rrt.start_with_rotation = start_with_rotation;
                    rrt.rotation_bias = rotation_bias;
                    rrt.gradient_bias = gradient_bias;
                    rrt.k_expension = k_expension;
                    rrt.max_depth = max_depth;

                    auto start_time = std::chrono::high_resolution_clock::now();
                    Viewer.setTree(rrt.start_tree, rrt.goal_tree);
                    rrt.findConstrainedPath(cp, iterations, Viewer);
                    auto end_time = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
                    if(rrt.path_length_quota > 0){
                        path_lengths.push_back(rrt.path_length_quota);
                        durations.push_back(duration);
                        tree_sizes.push_back(rrt.tree_size_quota);
                        expension_quotas.push_back(rrt.expension_quota);
                        iterations_quotas.push_back(rrt.iterations_quota);
                        successful_runs++;  
                    }
                    std::cout << "Run " << i + 1 << ": " << duration << " milliseconds, Path Length: " << rrt.path_length_quota << ", Tree Size: " << rrt.tree_size_quota << ", Expansion Quota: " << rrt.expension_quota << std::endl;
                }
                double average_duration = std::accumulate(durations.begin(), durations.end(), 0.0) / durations.size();
                double average_expansion_quota = std::accumulate(expension_quotas.begin(), expension_quotas.end(), 0.0) / expension_quotas.size();
                double average_path_length = std::accumulate(path_lengths.begin(), path_lengths.end(), 0.0) / path_lengths.size();
                double average_tree_size = std::accumulate(tree_sizes.begin(), tree_sizes.end(), 0.0) / tree_sizes.size();
                double average_iterations_quota = std::accumulate(iterations_quotas.begin(), iterations_quotas.end(), 0.0) / iterations_quotas.size();

                std::cout << std::endl << std::endl;
                std::cout << "Benchmarking completed: " << successful_runs << "/" << banchmark_runs << " successful runs." << std::endl;
                std::cout << "Average Duration: " << average_duration << " milliseconds" << std::endl;
                std::cout << "Average Iterations Quota: " << average_iterations_quota << std::endl;
                std::cout << "Average Tree Size: " << average_tree_size << std::endl;
                std::cout << "Average Expansion Quota: " << average_expansion_quota << std::endl;
                std::cout << "Average Path Length: " << average_path_length << std::endl;
                sort(iterations_quotas.begin(), iterations_quotas.end());
                std::cout << "Iteration Quotas: ";
                for (const auto& quota : iterations_quotas) {
                    std::cout << quota << " ";
                }
                std::cout << std::endl;

                auto now = std::chrono::system_clock::now();
                auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( now.time_since_epoch()).count();
                std::ofstream results_file("benchmark_results" + std::to_string(timestamp) + ".txt");
                results_file << "Knots: from " << start_file << " to " << goal_file << std::endl;
                results_file << "succsessrate  & average_duration & average_iterations_quota & average_tree_size & average_expansion_quota & average_path_length \\\\" << std::endl;
                results_file << successful_runs << "/" << banchmark_runs << " & " << average_duration << " & " << average_iterations_quota << " & " << average_tree_size << " & " << average_expansion_quota << " & " << average_path_length << " \\\\" << std::endl;
                results_file << std::endl << "Individual runs:" << std::endl;
                results_file << "Iterations Quota: ";
                for (const auto& quota : iterations_quotas) {
                    results_file << quota << ", ";
                }
                results_file << std::endl;
                results_file.close();
                benchmarking = false;
            }
            else{
                RRT rrt(start_dofs, goal_dofs, Viewer.pca, allPermutations,every_k_permutation);
                rrt.goal_bias = goalBias;
                rrt.near_goal_bias = near_goal_bias;
                rrt.tree_bias = tree_bias;
                rrt.max_energy = maxEnergy;
                rrt.pruning_interval = pruningInterval;
                rrt.steer_step = stepsize;
                rrt.step_length = steplength;
                rrt.radius = neighbor_radius;
                rrt.one_rand_direction_3d = oneRandDirection;
                rrt.goal_bias_for_all_permutations = goal_bias_for_all_permutations;
                rrt.sample_in_projection_space = sample_in_projection_space;
                rrt.use_constraint_projection_for_sampling = use_constraint_projection_for_sampling;
                rrt.reproject_direction = reproject_direction;
                rrt.sample_proj_dim = sample_proj_dim;
                rrt.max_rod_energy = max_rod_energy;
                rrt.k_neighbors = k_neighbors;
                rrt.start_with_rotation = start_with_rotation;
                rrt.rotation_bias = rotation_bias;
                rrt.gradient_bias = gradient_bias;
                rrt.k_expension = k_expension;
                rrt.max_depth = max_depth;

                std::cout << "Starting RRT with " << iterations << " iterations" << std::endl;
                auto start_time = std::chrono::high_resolution_clock::now();
                Viewer.setTree(rrt.start_tree, rrt.goal_tree);
                path = rrt.findConstrainedPath(cp, iterations, Viewer);
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
                std::cout << "RRT finished in " << duration << " milliseconds" << std::endl;
                std::cout << "Found a path of size: " << path.size() << std::endl; 
                showPath(path,cp,Viewer);
            }
            
            running=false;
        }
        if(show_path){
            showPath(path,cp,Viewer);
            show_path = false;
        }

    }


    return 0;
}