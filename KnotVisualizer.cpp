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

