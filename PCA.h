#include "helpers.h"
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>

std::vector<Eigen::VectorXd> removeTwist(const std::vector<Eigen::VectorXd> &data);

class PCA
{
public:

    void fit(const std::vector<Eigen::VectorXd>& samples);
    

    Eigen::VectorXd project(const Eigen::VectorXd& v, int k = 3);
    
private:
    Eigen::RowVectorXd mean;
    Eigen::MatrixXd components;
};