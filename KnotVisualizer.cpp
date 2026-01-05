#include "KnotVisualizer.h"

KnotVisualizer::KnotVisualizer() {
    polyscope::init();
    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
}

void KnotVisualizer::setKnot(const std::vector<Eigen::Vector3d>& pts, const double radius) {
    std::vector<std::array<size_t,2>> edges;
    for (size_t i = 0; i + 1 < pts.size(); ++i)
        edges.push_back({i, i+1});

    edges.push_back({0,pts.size()-1});

    curve = polyscope::registerCurveNetwork("knot", pts, edges);
    curve->setRadius(radius, true);

    std::vector<Eigen::Vector3d> grad(pts.size(), Eigen::Vector3d::Zero());
    Eigen::VectorXd twist = Eigen::VectorXd::Zero(edges.size());

    edge_twist = curve->addEdgeScalarQuantity("twist", twist);
    edge_twist->setEnabled(true);
    

    node_grad = curve->addNodeVectorQuantity("node gradient", grad);
    node_grad->setEnabled(true);
    node_grad->setVectorLengthScale(10);


    node_grad_mod = curve->addNodeVectorQuantity("node gradient modified", grad);
    node_grad_mod->setEnabled(true);
    node_grad_mod->setVectorLengthScale(10);
    

    contact_force = curve->addNodeVectorQuantity("contct force", grad);
    contact_force->setEnabled(true);
    contact_force->setVectorLengthScale(1);
    
}
void KnotVisualizer::setTheta(const std::vector<Eigen::Vector3d>& pts, const double radius) {

    std::vector<Eigen::Vector3d> grad_theta(pts.size(), Eigen::Vector3d::Zero());
    Eigen::VectorXd value_theta = Eigen::VectorXd::Zero(pts.size());

    theta = polyscope::registerPointCloud("theta", pts);
    theta->setPointRadius(radius,true);

    theta_value = theta->addScalarQuantity("theta value", value_theta);
    theta_value->setEnabled(true);

    theta_grad = theta->addVectorQuantity("theta gradient", grad_theta);
    theta_grad->setEnabled(true);
    theta_grad->setVectorLengthScale(100);

    theta_grad_mod = theta->addVectorQuantity("theta gradient modified", grad_theta);
    theta_grad_mod->setEnabled(true);
    theta_grad_mod->setVectorLengthScale(100);
}
void KnotVisualizer::colorTwist(Eigen::VectorXd& twist){
    edge_twist->updateData(twist);
    edge_twist->setMapRange({twist.minCoeff(),twist.maxCoeff()});
} 
void KnotVisualizer::showNodeGradient(const std::vector<Eigen::Vector3d>& grad){
    node_grad->updateData(grad);
} 
void KnotVisualizer::showNodeGradientModified(const std::vector<Eigen::Vector3d>& grad){
    node_grad_mod->updateData(grad);
} 
void KnotVisualizer::showContactForce(const std::vector<Eigen::Vector3d>& cForce){
    contact_force->updateData(cForce);
} 


void KnotVisualizer::updateKnot(const std::vector<Eigen::Vector3d>& pts) {
    if (!curve) return;
    curve->updateNodePositions(pts);
}

void KnotVisualizer::updateTheta(double value, double grad, double grad_mod) {
    if (!theta) return;
    theta_value->updateData(std::vector<double>{value,0});
    theta_value->setMapRange({-value,value});
    theta_grad->updateData(std::vector<Eigen::Vector3d>{ Eigen::Vector3d(0,0,grad),Eigen::Vector3d(0,0,grad) });
    theta_grad_mod->updateData(std::vector<Eigen::Vector3d>{ Eigen::Vector3d(0,0,grad_mod),Eigen::Vector3d(0,0,grad_mod) });
}



void KnotVisualizer::show() {
    polyscope::show();
}

void KnotVisualizer::frameTick() {
    polyscope::frameTick();
}

void KnotVisualizer::setUserCallback(std::function<void()> f) {
    userCallback = f;
    polyscope::state::userCallback = [this]() {
        if (userCallback) userCallback();
    };
}

