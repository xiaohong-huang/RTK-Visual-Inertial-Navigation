
#pragma once

#include <eigen3/Eigen/Dense>
#include "integration_base.h"
#include "../solver/cost_function.h"



class IMUFactor  {
  public:
    IMUFactor() = delete;
    IMUFactor(IntegrationBase* pre_integration_): pre_integration(pre_integration_) {
    }
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    IntegrationBase* pre_integration;


};


