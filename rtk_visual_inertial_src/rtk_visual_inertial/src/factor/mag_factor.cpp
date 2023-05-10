
#include "mag_factor.h"
#include "../parameter/parameters.h"
#include "../utility/utility.h"
#define MAGSQRTINFO 20.0
bool MAGFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    Eigen::Vector3d Pwb(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qwb(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Vector3d Bmg(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Vector3d Magw1 = Rwgw.transpose() * Qwb * (Mag + Bmg);

    double mag1_2_norm = Magw1.segment(0, 2).norm();
    residuals[0] = Magw1(0) / mag1_2_norm * MAGSQRTINFO;
    if (jacobians) {
        Eigen::Matrix3d Rwb = Rwgw.transpose() * Qwb.toRotationMatrix();
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);//b.segment(idx_i, size_i)
            jacobian_pose.setZero();
            jacobian_pose.segment(3, 3) = (Rwb * Utility::skewSymmetric(-Mag - Bmg) / mag1_2_norm).block(0, 0, 1, 3) * MAGSQRTINFO;

        }
        if (jacobians[1]) {
            Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_bias(jacobians[1]);//b.segment(idx_i, size_i)
            jacobian_bias.setZero();
            jacobian_bias = Rwb.block(0, 0, 1, 3) / mag1_2_norm * MAGSQRTINFO;
        }

    }



    return true;
}
