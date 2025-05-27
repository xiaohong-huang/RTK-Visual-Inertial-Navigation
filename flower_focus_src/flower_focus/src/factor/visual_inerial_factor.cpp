#include "visual_inerial_factor.h"
#include "../utility/tic_toc.h"



#define WAIT_BACKWARD_READY_FUNC(function){\
    if (!VI_info->next_vibase) {\
        while (1) {\
            function;\
            omp_set_lock(&omp_lock);\
            VI_info_pointer->ready_for_sequential_backward = true;\
            omp_unset_lock(&omp_lock);\
            if (VI_info_pointer->prev_vibase)\
                VI_info_pointer = VI_info_pointer->prev_vibase;\
            else {\
                VI_info_pointer = 0;\
                break;\
            }\
        }\
    }\
    while (1) {\
        omp_set_lock(&omp_lock);\
        if (VI_info->ready_for_sequential_backward) {\
            omp_unset_lock(&omp_lock);\
            break;\
        } else\
            omp_unset_lock(&omp_lock);\
    }\
    VI_info->ready_for_sequential_backward = false;\
}


#define WAIT_FORWARD_READY_FUNC(function){\
    if (!VI_info->prev_vibase) {\
        while (1) {\
            function;\
            omp_set_lock(&omp_lock);\
            VI_info_pointer->ready_for_sequential_forward = true;\
            omp_unset_lock(&omp_lock);\
            if (VI_info_pointer->next_vibase)\
                VI_info_pointer = VI_info_pointer->next_vibase;\
            else {\
                VI_info_pointer = 0;\
                break;\
            }\
        }\
    }\
    while (1) {\
        omp_set_lock(&omp_lock);\
        if (VI_info->ready_for_sequential_forward) {\
            omp_unset_lock(&omp_lock);\
            break;\
        } else\
            omp_unset_lock(&omp_lock);\
    }\
    VI_info->ready_for_sequential_forward = false;\
}



VisualInertialFactor::VisualInertialFactor(VisualInertialBase* VI_info_): VI_info(VI_info_) {
    int cnt = 0;
    for (int order = VI_info->O_POSR; order < VI_info->O_FULL; order++) {
        if (VI_info->order2pointers[order]) {
            ASSERT(VI_info->order2size[order]<10);
            mutable_parameter_block_sizes()->push_back(globalSize(VI_info->order2size[order]));
            cnt += VI_info->order2size[order];
        }
    }
    for (int i = 0; i < (int)VI_info->long_phase_bias.size(); i++) {
        mutable_parameter_block_sizes()->push_back(1);
        cnt += 1;
    }
    set_num_residuals(cnt);

}
bool  VisualInertialFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
    return false;
}

bool  VisualInertialFactor::EvaluateLhsRhs(double const* const* parameters, double* rhs, double* lhs, double* gradient, double mu)  {
    VI_info->EvaluateLhsRhs(parameters, rhs, lhs, gradient, mu);

    omp_set_lock(&omp_lock);
    VI_info->ready_for_sequential_forward = true;
    bool is_ok = false;
    ASSERT(VI_info_pointer);
    if (VI_info_pointer->ready_for_sequential_forward) {
        VI_info_pointer->ready_for_sequential_forward = false;
        is_ok = true;
    }
    omp_unset_lock(&omp_lock);

    while (is_ok) {
        is_ok = false;
        VI_info_pointer->UpdateLhsRhsGradientNew();
        omp_set_lock(&omp_lock);

        if (VI_info_pointer->next_vibase) {
            VI_info_pointer = VI_info_pointer->next_vibase;
            if (VI_info_pointer->ready_for_sequential_forward) {
                VI_info_pointer->ready_for_sequential_forward = false;
                is_ok = true;
            }
        } else
            VI_info_pointer = 0;

        omp_unset_lock(&omp_lock);
    }

    return true;
}


bool  VisualInertialFactor::UpdateInsideGaussStep(double const* const* parameters, double* gradient, double* gauss,
                                                  double* gauss_newton_step_inside_square_norm, double* alpha1, double* gradient_squared_norm, double* gradient_dot_gauss_newton_inside)  {

    WAIT_BACKWARD_READY_FUNC(VI_info_pointer->UpdateGaussAndGradientOutside())
    VI_info->SaveGaussStepOutside();
    VI_info->SaveGradientOutside();
    VI_info->UpdateInsideGaussStep();
    VI_info->SaveGaussStep();
    VI_info->EvaluateAlpha(alpha1);




    gauss_newton_step_inside_square_norm[0] = VI_info->gauss_newton_step_inside_square_norm;
    gradient_dot_gauss_newton_inside[0] = VI_info->gradient_dot_gauss_newton_inside;
    gradient_squared_norm[0] = VI_info->gradient_squared_norm_inside;

    return true;
}

bool  VisualInertialFactor::UpdateInsideStateAndEvaluateCost(double const* const* parameters, double* cost, double* model_cost_change,
                                                             double gradient_scale, double gauss_scale)  {

    WAIT_BACKWARD_READY_FUNC(VI_info_pointer->UpdateStateNew(gauss_scale, gradient_scale))
    VI_info->UpdateInsideStateUseGradientAndGauss(gradient_scale, gauss_scale);
    VI_info->EvaluateCost(0, model_cost_change);
    return true;
}

bool  VisualInertialFactor::EvaluateCost(double const* const* parameters, double* cost)  {
    WAIT_FORWARD_READY_FUNC(VI_info_pointer->ForwardDeliverIdepthValues())
    VI_info->EvaluateCost(cost, 0);
    return true;
}


bool  VisualInertialFactor::SaveOrRestoreHidenStates(bool is_save)  {
    VI_info->SaveOrRestoreHidenStates(is_save);
    VI_info->SaveLoadCandidateResidual(is_save);
    if (is_save)VI_info->EvaluateNonlinearQuantity();

    return true;
}


