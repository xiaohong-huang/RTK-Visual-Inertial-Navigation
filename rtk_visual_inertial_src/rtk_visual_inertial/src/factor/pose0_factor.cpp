
#include "pose0_factor.h"
#include "../utility/utility.h"



bool InitPose0Factor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    Eigen::Vector3d Pwb(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qwb(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Vector3d Magw1 = Rwgw.transpose() * Qwb * Mag;
    double mag1_2_norm = Magw1.segment(0, 2).norm();
    memset(residuals, 0, sizeof(double) * 4);

    if (posefix) {
        residuals[0] = (Pwb.x() - InitPwb.x()) * sqrt_info;
        residuals[1] = (Pwb.y() - InitPwb.y()) * sqrt_info;
        residuals[2] = (Pwb.z() - InitPwb.z()) * sqrt_info;
    }
    if (orientationfix) {
        residuals[3] = Magw1(0) / mag1_2_norm * sqrt_info;
    }

    if (jacobians) {
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
            jacobian_pose.setZero();
            if (posefix) {
                jacobian_pose.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * sqrt_info;
            }
            if (orientationfix) {
                jacobian_pose.block<1, 3>(3, 3) = (Rwgw.transpose() * Qwb.toRotationMatrix() *
                                                   Utility::skewSymmetric(-Mag) / mag1_2_norm).block(0, 0, 1, 3) * sqrt_info;
            }
        }
    }



    return true;
}

