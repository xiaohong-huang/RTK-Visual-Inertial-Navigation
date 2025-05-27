
#pragma once
#include <eigen3/Eigen/Dense>
#include "../solver/cost_function.h"
#include "../parameter/parameters.h"
#include "visual_inerial_base.h"


class VisualInertialFactor: public CostFunction {
  public:
    VisualInertialFactor(VisualInertialBase* );
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)const;
    bool  EvaluateLhsRhs(double const* const* parameters, double* rhs, double* lhs, double* gradient, double mu);
    bool  UpdateInsideGaussStep(double const* const* parameters, double* gradient, double* gauss,
                                double* gauss_newton_step_inside_square_norm, double* alpha1, double* gradient_squared_norm, double* gradient_dot_gauss_newton_inside);
    bool  UpdateInsideStateAndEvaluateCost(double const* const* parameters, double* cost, double* model_cost_change,
                                           double gradient_scale, double gauss_scale);
    bool SaveOrRestoreHidenStates(bool is_save);
    bool  EvaluateCost(double const* const* parameters, double* cost);
    VisualInertialBase* VI_info;
};





