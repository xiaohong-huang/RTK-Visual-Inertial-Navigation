
#pragma once

#include <eigen3/Eigen/Dense>
#include "integration_base.h"
#include <ceres/ceres.h>

class IMUFactor : public ceres::SizedCostFunction<15, 7, 9, 7, 9> {
  public:
    IMUFactor() = delete;
    IMUFactor(IntegrationBase* _pre_integration): pre_integration(_pre_integration) {
    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    bool Evaluate2(double const* const* parameters, double* residuals, double** jacobians) ;

    IntegrationBase* pre_integration;

};


