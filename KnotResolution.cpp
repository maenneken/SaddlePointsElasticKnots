#include "helpers.h"

Eigen::VectorXd doubleKnotResolution(Eigen::VectorXd DoFs){
    Eigen::VectorXd double_DoFs =  Eigen::VectorXd::Zero(DoFs.size()-1);
    size_t n_pts = Dofs.size() / 4;
    //3d Points
    for ( size_t i = 0; i < n_pts; ++i){
        size_t i1 = (i+1) % n_pts; //wrap to start
        double_DoFs.segment<3>(6*i) = Dofs.segment<3>(3*i);
        double_DoFs.segment<3>(6*i1) = (Dofs.segment<3>(3*i1) - Dofs.segment<3>(3*i1))/2;
    }
    //twist
    for ( size_t i = 3*n_pts; i < 4*n_pts; ++i){
        size_t i1 = (i+1) % n_pts; //wrap to start
        double_DoFs(6*i) = Dofs.segment<3>(3*i); //same twist
        double_DoFs(6*i1) = (Dofs.segment<3>(3*i1) + Dofs.segment<3>(3*i1))/2; //mean twist
    }
    //theta
    double_DoFs(8*n_pts) = DoFs(4*n_pts);
    return double_DoFs;

}

Eigen::VectorXd halfKnotResolution(Eigen::VectorXd DoFs){
    Eigen::VectorXd double_DoFs =  Eigen::VectorXd::Zero(DoFs.size()-1);
    size_t n_pts = Dofs.size() / 4;
    //3d Points
    for ( size_t i = 0; i < n_pts; ++i){
        size_t i1 = (i+1) % n_pts; //wrap to start
        double_DoFs.segment<3>(6*i) = Dofs.segment<3>(3*i);
        double_DoFs.segment<3>(6*i1) = (Dofs.segment<3>(3*i1) - Dofs.segment<3>(3*i1))/2;
    }
    //twist
    for ( size_t i = 3*n_pts; i < 4*n_pts; ++i){
        size_t i1 = (i+1) % n_pts; //wrap to start
        double_DoFs(6*i) = Dofs.segment<3>(3*i); //same twist
        double_DoFs(6*i1) = (Dofs.segment<3>(3*i1) + Dofs.segment<3>(3*i1))/2; //mean twist
    }
    //theta
    double_DoFs(8*n_pts) = DoFs(4*n_pts);
    return double_DoFs;

}