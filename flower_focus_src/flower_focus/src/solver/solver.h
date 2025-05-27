

#pragma once

#include <cstdlib>
#include <pthread.h>
#include <unordered_map>
#include "solver_residualblock.h"
#include "../parameter/parameters.h"


class SolverInfo {
  public:
    SolverInfo() {};


    ~SolverInfo();
    void addResidualBlockInfo(SolverResidualBlockInfo* residual_block_info);
    void solve(int num_iteration);
    void init_solver();
    void ConstructHessianRhs();
    void UpdateGaussAndGradientStep();
    void StepAccepted(double step_quality);
    void StepRejected(double step_quality);
    void SolverEvaluateCost();
    void SaveOrRestoreStates(bool is_save);
    void DenseSolver();
    double StepQuality(const double cost, const double model_cost_change);
    void StepAccepted(const double cost, const double model_cost_change);
    bool IsStepSuccessful();
    void ComputeTraditionalDoglegStep();
    void StepIsInvalid();
    void UpdateParameters();
    bool FunctionToleranceReached();
    void marginalization_process();
    std::vector<SolverResidualBlockInfo*> factors;

    std::unordered_map<double*, int> parameter_block_size; //global size
    std::unordered_map<double*, double*> parameter_block_data;



    int num_rows;
    int pos = 0;
    int num_consecutive_nonmonotonic_steps_;
    int max_consecutive_nonmonotonic_steps_ = 0;


    double last_cost = 0;
    double decrease_factor;
    double cost = 0;
    double model_cost_change = 0;
    double minimum_cost_;
    double current_cost_;
    double reference_cost_;
    double candidate_cost_;
    double accumulated_reference_model_cost_change_;
    double accumulated_candidate_model_cost_change_;
    double relative_decrease = 0;
    double gauss_newton_square_norm = 0;
    double gauss_newton_norm = 0;
    double gradient_square_norm = 0;
    double gradient_norm = 0;
    double alpha1 = 0;
    double alpha = 0;
    double dogleg_step_norm = 0;
    double gradient_dot_gauss_newton_inside = 0;
    double gauss_newton_scale = 0;
    double gradient_scale = 0;
    double radius = 0;
    double mu = 0;
    double x_norm_ = 0;

    bool is_sucess = false;
    bool reuse_ = 0;


};

