

#include "solver.h"
#include "../parameter/parameters.h"
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "omp.h"







void SolverResidualBlockInfo::Evaluate(double mu) {

    cost_function->EvaluateLhsRhs(parameter_blocks.data(), 0,0,0, mu);
    ASSERT(loss_function == 0);


}



void SolverResidualBlockInfo::FactorUpdateInsideGaussStep(double& gauss_newton_step_inside_square_norm0,
                                                          double& gradient_square_norm0,
                                                          double& alpha0,
                                                          double& gradient_dot_gauss_newton_inside0) {

    double alpha1 = 0;
    double gradient_square_norm = 0;
    double gradient_dot_gauss_newton_inside = 0;
    double gauss_newton_step_inside_square_norm = 0;


    cost_function->UpdateInsideGaussStep(parameter_blocks.data(),
                                                  0,0, &gauss_newton_step_inside_square_norm,
                                                  &alpha1, &gradient_square_norm, &gradient_dot_gauss_newton_inside);

    gauss_newton_step_inside_square_norm0 += gauss_newton_step_inside_square_norm;
    gradient_square_norm0 += gradient_square_norm;
    gradient_dot_gauss_newton_inside0 += gradient_dot_gauss_newton_inside;
    alpha0 += alpha1;
}



void SolverResidualBlockInfo::FactorUpdateInsideStateAndEvaluateModelCostChange(double gradient_scale, double gauss_scale, double& model_cost_change0) {

    double model_cost_change = 0;
    cost_function->UpdateInsideStateAndEvaluateCost(parameter_blocks.data(),
                                                             0, &model_cost_change, gradient_scale, gauss_scale);

    model_cost_change0 += model_cost_change;

}



void SolverResidualBlockInfo::FactorEvaluateCost(double& cost0) {
    double cost = 0;
    cost_function->EvaluateCost(parameter_blocks.data(), &cost);

    cost0 += cost;

}

