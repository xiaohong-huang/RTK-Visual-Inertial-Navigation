

#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>



//feature poistion factor
class projection_factor : public ceres::SizedCostFunction<2, 7, 7, 3> {
  public:
    projection_factor( const Eigen::Vector3d& _pts);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d pts;
    static Eigen::Matrix2d sqrt_info;
};

//feature inverse depth factor
class ProjectionOneFrameTwoCamFactor : public ceres::SizedCostFunction<2, 7, 7, 1> {
  public:
    ProjectionOneFrameTwoCamFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d pts_i, pts_j;
    static Eigen::Matrix2d sqrt_info;
};


//feature inverse depth factor
class ProjectionTwoFrameOneCamFactor : public ceres::SizedCostFunction<2, 7, 7, 7, 1> {
  public:
    ProjectionTwoFrameOneCamFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d pts_i, pts_j;
    static Eigen::Matrix2d sqrt_info;
};


//feature inverse depth factor
class ProjectionTwoFrameTwoCamFactor : public ceres::SizedCostFunction<2, 7, 7, 7, 7, 1> {
  public:
    ProjectionTwoFrameTwoCamFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d pts_i, pts_j;
    static Eigen::Matrix2d sqrt_info;
};
