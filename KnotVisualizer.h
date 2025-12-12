#pragma once
#include <polyscope/curve_network.h>
#include <Eigen/Dense>
#include <vector>

class KnotVisualizer {
public:
    KnotVisualizer();
    void setKnot(const std::vector<Eigen::Vector3d>& pts, const double radius);
    void colorTwist(Eigen::VectorXd& twist);
    void updateKnot(const std::vector<Eigen::Vector3d>& pts);
    void showNodeGradient(const std::vector<Eigen::Vector3d>& grad);
    void showContactForce(const std::vector<Eigen::Vector3d>& cForce);
    void show();
    void frameTick();
    void setUserCallback(std::function<void()> f);

private:
    polyscope::CurveNetwork* curve = nullptr;
    polyscope::CurveNetworkNodeVectorQuantity* node_grad = nullptr;
    polyscope::CurveNetworkEdgeScalarQuantity* edge_twist = nullptr;
    polyscope::CurveNetworkNodeVectorQuantity* contact_force = nullptr;
    std::function<void()> userCallback;
};