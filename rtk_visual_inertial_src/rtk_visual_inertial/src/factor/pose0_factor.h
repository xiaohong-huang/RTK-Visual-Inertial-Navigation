
#pragma once
#include <Eigen/Dense>
#include <ceres/ceres.h>
#include "../parameter/parameters.h"

class InitPose0Factor : public ceres::SizedCostFunction<4, 7> {
  public:
    InitPose0Factor() = delete;
    InitPose0Factor(const Eigen::MatrixXd Rwb, Eigen::Vector3d Pwb, bool posefix_, bool orientationfix_,
                    double sqrt_info_): InitPwb(Pwb),
        posefix(posefix_), orientationfix(orientationfix_), sqrt_info(sqrt_info_) {
        Magw(0) = 0;
        Magw(1) = 1;
        Magw(2) = 0;
        Mag = Rwb.transpose() * Rwgw * Magw;
    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    Eigen::Vector3d InitPwb;
    Eigen::Vector3d Magw;
    Eigen::Vector3d Mag;
    bool posefix, orientationfix;
    double sqrt_info;

};

