

#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>

class InitialPoseFactor : public ceres::SizedCostFunction<6, 7> {
  public:
    InitialPoseFactor(const Eigen::Vector3d& _P, const Eigen::Quaterniond& _Q, Eigen::Matrix<double, 6, 6>sqrt_info_);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    Eigen::Vector3d init_P;
    Eigen::Quaterniond init_Q;
    Eigen::Matrix<double, 6, 6> sqrt_info;
};

class InitialFactor33 : public ceres::SizedCostFunction<3, 3> {
  public:
    InitialFactor33(const Eigen::Vector3d& _P, Eigen::Matrix<double, 3, 3>sqrt_info_);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    Eigen::Vector3d init_P;
    Eigen::Matrix<double, 3, 3> sqrt_info;
};


class InitialFactor11 : public ceres::SizedCostFunction<1, 1> {
  public:
    InitialFactor11(double P_, double sqrt_info_);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    double init_P;
    double sqrt_info;
};



class InitialBlackFactor : public ceres::SizedCostFunction<1, 1> {
  public:
    InitialBlackFactor(double istd_): istd(istd_) {
    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;
    double istd;
};


class InitialBiasFactor : public ceres::SizedCostFunction<9, 9> {
  public:
    InitialBiasFactor(const Eigen::Vector3d& _V, const Eigen::Vector3d& _Ba, const Eigen::Vector3d& _Bg, const Eigen::Matrix<double, 9, 9>sqrt_info_);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d init_Ba, init_Bg, init_V;
    Eigen::Matrix<double, 9, 9> sqrt_info;
};
