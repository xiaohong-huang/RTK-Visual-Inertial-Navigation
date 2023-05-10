#pragma once
#include <Eigen/Dense>
#include <ceres/ceres.h>
// magnetormeter factor
class MAGFactor : public ceres::SizedCostFunction<1, 7, 3> {
  public:
    MAGFactor() = delete;
    MAGFactor(const Eigen::Vector3d& Mag_): Mag(Mag_) {}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d Mag;
};

