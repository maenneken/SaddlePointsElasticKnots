#include "helpers.h"
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>

std::vector<Eigen::VectorXd> removeTwist(const std::vector<Eigen::VectorXd> &data);
std::vector<Eigen::VectorXd> createDataset(std::string folder,  int reductionFactor);
class PCA
{
public:

    void fit(const std::vector<Eigen::VectorXd>& samples);
    Eigen::VectorXd project(const Eigen::VectorXd& v, int k = 3);
    Eigen::VectorXd reconstruct(const Eigen::VectorXd& y, int k);
    double dist(const Eigen::VectorXd& x, const Eigen::VectorXd& y, int k = 3);
    
private:
    Eigen::RowVectorXd mean;
    Eigen::MatrixXd components;
    size_t data_dim;
};