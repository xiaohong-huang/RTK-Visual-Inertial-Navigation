

#pragma once

#include <cstdlib>
#include <pthread.h>
#include "../solver/cost_function.h"
#include <unordered_map>
#include "../parameter/parameters.h"
#include "../factor/loss_function.h"
class SolverResidualBlockInfo {
  public:
    ~SolverResidualBlockInfo() {
        delete cost_function;

    }
    SolverResidualBlockInfo(CostFunction* _cost_function, LossFunction* _loss_function, std::vector<double*> _parameter_blocks)
        : cost_function(_cost_function), loss_function(_loss_function), parameter_blocks(_parameter_blocks) {
        std::vector<int> block_sizes = cost_function->parameter_block_sizes();
    }


    void Evaluate(double mu);
    void FactorUpdateInsideGaussStep(double& gauss_newton_step_inside_square_norm0,
                                     double& gradient_square_norm0,
                                     double& alpha0,
                                     double& gradient_dot_gauss_newton_inside0);
    void FactorUpdateInsideStateAndEvaluateModelCostChange(double gradient_scale, double gauss_scale, double& model_cost_change0);
    void FactorEvaluateCost(double& cost0);

    CostFunction* cost_function;
    LossFunction* loss_function;
    std::vector<double*> parameter_blocks;

};


