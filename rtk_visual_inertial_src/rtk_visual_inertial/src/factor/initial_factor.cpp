#include "initial_factor.h"
#include "../utility/utility.h"



bool InitialPoseFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Q(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
    residual.block<3, 1>(0, 0) = P - init_P;
    residual.block<3, 1>(3, 0) = 2 * (init_Q.inverse() * Q).vec();
    residual = sqrt_info * residual;

    if (jacobians) {
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
            jacobian_pose.setZero();
            jacobian_pose.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
            jacobian_pose.block<3, 3>(3, 3) = Utility::Qleft(init_Q.inverse() * Q).bottomRightCorner<3, 3>();
            jacobian_pose = sqrt_info * jacobian_pose;
        }

    }
    return true;

}

InitialPoseFactor::InitialPoseFactor(const Eigen::Vector3d& _P, const Eigen::Quaterniond& _Q, Eigen::Matrix<double, 6, 6>sqrt_info_) {
    init_P = _P;
    init_Q = _Q;
    sqrt_info = sqrt_info_;
}




bool InitialFactor33::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    Eigen::Vector3d P(parameters[0][0], parameters[0][1], parameters[0][2]);

    Eigen::Map<Eigen::Matrix<double, 3, 1>> residual(residuals);
    residual.block<3, 1>(0, 0) = P - init_P;
    residual = sqrt_info * residual;

    if (jacobians) {
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
            jacobian_pose.setZero();
            jacobian_pose.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
            jacobian_pose = sqrt_info * jacobian_pose;
        }

    }
    return true;

}
InitialFactor33::InitialFactor33(const Eigen::Vector3d& _P, Eigen::Matrix<double, 3, 3>sqrt_info_) {
    init_P = _P;
    sqrt_info = sqrt_info_;
}


bool InitialFactor11::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {

    residuals[0] = (parameters[0][0] - init_P) * sqrt_info;
    if (jacobians) {
        if (jacobians[0]) {
            jacobians[0][0] = sqrt_info;
        }
    }
    return true;

}
InitialFactor11::InitialFactor11(double P_, double sqrt_info_) {
    init_P = P_;
    sqrt_info = sqrt_info_;
}

bool InitialBlackFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    residuals[0] = parameters[0][0] * istd;
    if (jacobians) {
        if (jacobians[0])jacobians[0][0] = 1 * istd;
    }
    return true;
}



bool InitialBiasFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    Eigen::Vector3d    V(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Vector3d Ba(parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Vector3d Bg(parameters[0][6], parameters[0][7], parameters[0][8]);

    Eigen::Map<Eigen::Matrix<double, 9, 1>> residual(residuals);
    residual.block<3, 1>(0, 0) = V - init_V;
    residual.block<3, 1>(3, 0) = Ba - init_Ba;
    residual.block<3, 1>(6, 0) = Bg - init_Bg;
    residual = sqrt_info * residual;

    if (jacobians) {
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 9, 9, Eigen::RowMajor>> jacobian(jacobians[0]);
            jacobian.setZero();
            jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
            jacobian.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
            jacobian.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity();
            jacobian = sqrt_info * jacobian;
        }
    }
    return true;
}

InitialBiasFactor::InitialBiasFactor(const Eigen::Vector3d& _V, const Eigen::Vector3d& _Ba, const Eigen::Vector3d& _Bg, const Eigen::Matrix<double, 9, 9>sqrt_info_) {
    init_Ba = _Ba;
    init_Bg = _Bg;
    init_V = _V;
    sqrt_info = sqrt_info_;
}