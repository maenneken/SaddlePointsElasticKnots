#include "KnotResolution.h"

//double the resolution of the Knot
Eigen::VectorXd doubleKnotResolution(Eigen::VectorXd DoFs){
    Eigen::VectorXd double_DoFs =  Eigen::VectorXd::Zero(2*(DoFs.size())-1);
    size_t n_pts = DoFs.size() / 4;
    //3d Points
    for ( size_t i = 0; i < n_pts; ++i){
        size_t i1 = (i+1) % n_pts; //wrap to start
        double_DoFs.segment<3>(3*(2*i)) = DoFs.segment<3>(3*i);
        double_DoFs.segment<3>(3*(2*i +1)) = (DoFs.segment<3>(3*i) + DoFs.segment<3>(3*i1)) / 2.0;;
    }
    //twist
    size_t offset = 3*n_pts;
    for ( size_t i = 0; i < n_pts; ++i){
        size_t i1 = (i+1) % n_pts; //wrap to start
        double_DoFs(2*(i+offset)) = DoFs(i+offset); //same twist
        double_DoFs(2*(i +offset)+1) = (DoFs(i+offset) + DoFs(i1+offset))/2; //mean twist
    }
    //theta
    double_DoFs(8*n_pts) = DoFs(4*n_pts);
    return double_DoFs;

}

//half the resolution of the knot, whyle keeping topology
//not implemented yet
Eigen::VectorXd halfKnotResolution(Eigen::VectorXd DoFs){
    Eigen::VectorXd double_DoFs =  Eigen::VectorXd::Zero(DoFs.size()-1);
    return double_DoFs;

}