

#include "solver.h"
#include "../parameter/parameters.h"
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "omp.h"
#include "../factor/visual_inerial_base.h"
const double max_radius = 1e32;
const double min_relative_decrease = 1e-3;
const double kMaxMu = 1e6;

const double function_tolerance = 1e-6;
const double parameter_tolerance = 1e-8;
const double decrease_threshold_ = 0.25;
const double increase_threshold_ = 0.75;
const double mu_increase_factor_ = 10.0;


SolverInfo::~SolverInfo() {


    for (int i = 0; i < (int)factors.size(); i++)
        delete factors[i];
    for (auto it = parameter_block_data.begin(); it != parameter_block_data.end(); it++)
        delete it->second;
}

void SolverInfo::addResidualBlockInfo(SolverResidualBlockInfo* residual_block_info) {
    factors.emplace_back(residual_block_info);

    std::vector<double*>& parameter_blocks = residual_block_info->parameter_blocks;
    std::vector<int> parameter_block_sizes = residual_block_info->cost_function->parameter_block_sizes();
    ASSERT(parameter_blocks.size() == parameter_block_sizes.size());

    for (int i = 0; i < static_cast<int>(residual_block_info->parameter_blocks.size()); i++) {
        double* addr = parameter_blocks[i];
        int size = parameter_block_sizes[i];
#if USE_ASSERT
        if (parameter_block_size.find(addr) != parameter_block_size.end())
            // if (parameter_block_size[addr] != size) {
            //     LOG_OUT << size << "," << parameter_block_size[addr] << addr << std::endl;
            ASSERT(parameter_block_size[addr] == size);
        // }
#endif
        parameter_block_size[addr] = size;

    }
}



void SolverInfo::ConstructHessianRhs() {
    VI_info_pointer = VI_info_pointer_head;
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) if (NUM_THREADS > 1)
    for (int i = 0; i < (int)factors.size(); i++)
        factors[i]->Evaluate(mu);
    ASSERT(VI_info_pointer == 0);
}


void SolverInfo::init_solver() {

    num_rows = pos;

    accumulated_reference_model_cost_change_ = 0;
    accumulated_candidate_model_cost_change_ = 0;
    radius = initial_radius;
    decrease_factor = 2;
    mu = kMinMu;



    for (auto it = parameter_block_size.begin(); it != parameter_block_size.end(); it++) {
        double* pointer = it->first;
        int size = it->second;
        double* data = new double[size];
        parameter_block_data[pointer] = data;
    }
}

void SolverInfo::SolverEvaluateCost() {

    std::vector<double>costs(factors.size(), 0);
    VI_info_pointer = VI_info_pointer_head;
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) if (NUM_THREADS > 1)
    for (int i = 0; i < (int)factors.size(); i++)
        factors[i]->FactorEvaluateCost(costs[i]);
    ASSERT(VI_info_pointer == 0);


    cost = 0;
    for (int i = 0; i < (int)costs.size(); i++)
        cost += costs[i];

    // ASSERT(cost<1e6&&cost>0);
}

void SolverInfo::UpdateGaussAndGradientStep() {



    std::vector<double>gauss_newton_step_inside_square_norms(factors.size(), 0);
    std::vector<double>gradient_square_norms(factors.size(), 0);
    std::vector<double>alphas(factors.size(), 0);
    std::vector<double>gradient_dot_gauss_newton_insides(factors.size(), 0);


    VI_info_pointer = VI_info_pointer_tail;
    ASSERT(!VI_info_pointer->next_vibase);
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)  if (NUM_THREADS > 1)
    for (int i = (int)factors.size() - 1; i >= 0 ; i--) {
        factors[i]->FactorUpdateInsideGaussStep(gauss_newton_step_inside_square_norms[i],
                                                gradient_square_norms[i], alphas[i],
                                                gradient_dot_gauss_newton_insides[i]);
    }
    ASSERT(VI_info_pointer == 0);

    gauss_newton_square_norm = 0;
    gradient_square_norm = 0;
    alpha1 = 0;
    gradient_dot_gauss_newton_inside = 0;

    for (int i = 0; i < (int)factors.size(); i++) {
        gauss_newton_square_norm += gauss_newton_step_inside_square_norms[i];
        gradient_square_norm += gradient_square_norms[i];
        alpha1 += alphas[i];
        gradient_dot_gauss_newton_inside += gradient_dot_gauss_newton_insides[i];
    }

    gauss_newton_norm = sqrt(gauss_newton_square_norm);
    gradient_norm = sqrt(gradient_square_norm);
    alpha = gradient_square_norm / alpha1;

    ASSERT(!isnan(gauss_newton_square_norm));
    ASSERT(!isnan(gradient_norm));
    ASSERT(!isinf(gauss_newton_square_norm));
    ASSERT(!isinf(gradient_norm));


}

void SolverInfo::UpdateParameters() {


    ComputeTraditionalDoglegStep();

    std::vector<double>model_cost_changes(factors.size(), 0);

    VI_info_pointer = VI_info_pointer_tail;
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)  if (NUM_THREADS > 1)
    for (int i = (int)factors.size() - 1; i >= 0 ; i--)
        factors[i]->FactorUpdateInsideStateAndEvaluateModelCostChange(gradient_scale, gauss_newton_scale, model_cost_changes[i]);
    ASSERT(VI_info_pointer == 0);

    model_cost_change = 0;
    for (int i = 0; i < (int)model_cost_changes.size(); i++)
        model_cost_change += model_cost_changes[i];

    SolverEvaluateCost();
    ASSERT(!isnan(model_cost_change));
    // ASSERT(model_cost_change > 0);




}

void SolverInfo::ComputeTraditionalDoglegStep() {


    if (gauss_newton_norm <= radius) {
        gauss_newton_scale = 1;
        gradient_scale = 0;
        dogleg_step_norm = gauss_newton_norm;
        return;
    }

    // Case 2. The Cauchy point and the Gauss-Newton steps lie outside
    // the trust region. Rescale the Cauchy point to the trust region
    // and return.
    if (gradient_norm * alpha >= radius) {
        gradient_scale = (-radius / gradient_norm);
        gauss_newton_scale = 0;
        dogleg_step_norm = gradient_norm;
        return;
    }

    // Case 3. The Cauchy point is inside the trust region and the
    // Gauss-Newton step is outside. Compute the line joining the two
    // points and the point on it which intersects the trust region
    // boundary.

    // a = alpha * -gradient
    // b = gauss_newton_step// gradient.dot(inc);
    double gradient_dot_gauss_newton = gradient_dot_gauss_newton_inside;
    const double b_dot_a = alpha * gradient_dot_gauss_newton;
    const double a_squared_norm = pow(alpha * gradient_norm, 2.0);
    const double b_minus_a_squared_norm =
        a_squared_norm - 2 * b_dot_a + pow(gauss_newton_norm, 2);

    // c = a' (b - a)
    //   = alpha * -gradient' gauss_newton_step - alpha^2 |gradient|^2
    const double c = b_dot_a - a_squared_norm;
    const double d = sqrt(c * c + b_minus_a_squared_norm *
                          (pow(radius, 2.0) - a_squared_norm));

    double beta = (c <= 0) ? (d - c) / b_minus_a_squared_norm
                  : (radius * radius - a_squared_norm) / (d + c);
    gradient_scale = (-alpha * (1.0 - beta));
    gauss_newton_scale = beta;
    // //a*A+b*B==aaAA+bbBB+2abAB
    // double a=(gradient_norm*gradient_scale);
    // double b=(gauss_newton_scale*gauss_newton_norm);
    // dogleg_step_norm=sqrt((a+b)*(a+b));
    dogleg_step_norm = sqrt(gradient_scale * gradient_scale * gradient_square_norm +
                            gauss_newton_scale * gauss_newton_scale * gauss_newton_square_norm -
                            2 * gradient_scale * gauss_newton_scale * gradient_dot_gauss_newton);
    // LOG_OUT<<dogleg_step_norm<<","<<gradient_scale<<","<<gauss_newton_scale<<std::endl;
}



void SolverInfo::SaveOrRestoreStates(bool is_save) {
    // if(is_save)LOG_OUT<<"SAVE\r\n";
    // else LOG_OUT<<"RESTORE\r\n";
    for (auto it = parameter_block_data.begin(); it != parameter_block_data.end(); it++) {
        double* src, *dst;
        if (is_save) {
            src = it->first;
            dst = it->second;
        } else {
            src = it->second;
            dst = it->first;
        }
        int size = parameter_block_size[it->first];
        memcpy(dst, src, size * sizeof(double));
    }
    for (int i = 0; i < (int)factors.size() ; i++)
        factors[i]->cost_function->SaveOrRestoreHidenStates(is_save);



}
void SolverInfo::DenseSolver() {
    is_sucess = true;
}


void SolverInfo::marginalization_process() {

    SolverEvaluateCost();
    SaveOrRestoreStates(Save);
    ConstructHessianRhs();//A,b

}

void SolverInfo::solve(int num_iteration) {



    static double ts[10] = {0};
    static int count = 0;
    count++;
    TicToc ta;
    SolverEvaluateCost();
    SaveOrRestoreStates(Save);
    ts[0] += ta.toc(); ta.tic();
    last_cost = cost;
    minimum_cost_ = last_cost;
    current_cost_ = last_cost;
    reference_cost_ = last_cost;
    candidate_cost_ = last_cost;
    LOG_OUT << "init cost:" << cost << std::endl;
    ts[1] += ta.toc(); ta.tic();
    static int count2 = 0;
    static int count3 = 0;
    for (int i = 0; i < num_iteration; i++) {
        if (!reuse_) {
            count2++;
            ConstructHessianRhs();
            ts[2] += ta.toc(); ta.tic();
            DenseSolver();
            ts[3] += ta.toc(); ta.tic();
            if (!is_sucess) {
                mu *= mu_increase_factor_;
                continue;
            }
            UpdateGaussAndGradientStep();
            ts[4] += ta.toc(); ta.tic();
        }


        UpdateParameters();
        ts[5] += ta.toc(); ta.tic();

        LOG_OUT << i << alpha << "," << gradient_norm << "," << gauss_newton_norm
                << "," << gauss_newton_scale << "," << gradient_scale << ",<<" << cost << ">>," << model_cost_change
                << "," << radius << "," << mu << "," << relative_decrease << std::endl;

        if (model_cost_change <= 0) { //invalid step
            StepIsInvalid();
            SaveOrRestoreStates(Restore);
            count3++;
            continue;
        }

        if (FunctionToleranceReached())
            break;


        if (IsStepSuccessful()) {
            SaveOrRestoreStates(Save);
            StepAccepted(relative_decrease);
            StepAccepted(cost, model_cost_change);
        } else {
            StepRejected(relative_decrease);
            SaveOrRestoreStates(Restore);
        }
        ts[6] += ta.toc(); ta.tic();
        last_cost = cost;
    }
    for (int i = 0; i < 7; i++)
        LOG_OUT << ts[i] / count << ",";
    LOG_OUT << std::endl;
    LOG_OUT << count2 << "," << count3;
    LOG_OUT << std::endl;


}



double SolverInfo::StepQuality(
    const double cost, const double model_cost_change)  {

    if (cost >= std::numeric_limits<double>::max())
        return std::numeric_limits<double>::lowest();

    const double relative_decrease = (current_cost_ - cost) / model_cost_change;
    const double historical_relative_decrease =
        (reference_cost_ - cost) /
        (accumulated_reference_model_cost_change_ + model_cost_change);
    return std::max(relative_decrease, historical_relative_decrease);
}

void SolverInfo::StepAccepted(const double cost,
                              const double model_cost_change) {


    current_cost_ = cost;
    accumulated_candidate_model_cost_change_ += model_cost_change;
    accumulated_reference_model_cost_change_ += model_cost_change;


    if (current_cost_ < minimum_cost_) {
        minimum_cost_ = current_cost_;
        num_consecutive_nonmonotonic_steps_ = 0;
        candidate_cost_ = current_cost_;
        accumulated_candidate_model_cost_change_ = 0.0;
    } else {
        ++num_consecutive_nonmonotonic_steps_;
        if (current_cost_ > candidate_cost_) {
            candidate_cost_ = current_cost_;
            accumulated_candidate_model_cost_change_ = 0.0;
        }
    }

    if (num_consecutive_nonmonotonic_steps_ ==
            max_consecutive_nonmonotonic_steps_) {
        reference_cost_ = candidate_cost_;
        accumulated_reference_model_cost_change_ =
            accumulated_candidate_model_cost_change_;
    }

}



bool SolverInfo::IsStepSuccessful() {
    relative_decrease = StepQuality(cost, model_cost_change);
    return  relative_decrease > min_relative_decrease;
}



void SolverInfo::StepAccepted(double step_quality) {
    //   CHECK_GT(step_quality, 0.0);

    if (step_quality < decrease_threshold_)
        radius *= 0.5;
    // LOG_OUT<<step_quality<<","<<decrease_threshold_<<","<<increase_threshold_<<","<<dogleg_step_norm<<","<<radius<<std::endl;
    if (step_quality > increase_threshold_)
        radius = std::max(radius, 3.0 * dogleg_step_norm);

    // Reduce the regularization multiplier, in the hope that whatever
    // was causing the rank deficiency has gone away and we can return
    // to doing a pure Gauss-Newton solve.
    mu = std::max(kMinMu, 2.0 * mu / mu_increase_factor_);
    reuse_ = false;
}

void SolverInfo::StepRejected(double step_quality) {
    radius *= 0.5;
    reuse_ = true;
}


void SolverInfo::StepIsInvalid() {
    mu *= mu_increase_factor_;
    reuse_ = false;
}



// Solver::Options::function_tolerance based convergence check.
bool SolverInfo::FunctionToleranceReached() {
    double cost_change = current_cost_ - cost;
    const double absolute_function_tolerance = function_tolerance * min(cost, current_cost_);

    if (fabs(cost_change) > absolute_function_tolerance)
        return false;

    LOG_OUT << "FunctionToleranceReached\r\n";
    return true;
}
