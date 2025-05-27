#include "visual_inerial_base.h"

#include "../utility/utility.h"
#include "operation.h"

#include "../utility/tic_toc.h"






void VisualInertialBase::SaveLoadCandidateResidual(bool is_save) {


    if (is_save) {



        if (last_marg_info) {
            prior_rhs_save = prior_rhs;
            prior_residual_save = prior_residual;
            if (!have_hist) {
                yaw_constraint_residual_save = yaw_constraint_residual;
                yaw_constraint_model_residual_accum_save = yaw_constraint_model_residual_accum;
            }
            if (ESTIMATE_SCALE && fix_scale && !have_hist) {
                init_distance_constraint_residual_save = init_distance_constraint_residual;
                init_distance_constraint_model_residual_accum_save = init_distance_constraint_model_residual_accum;
            }
        }

        memcpy(visual_residuals_save, visual_residuals,  sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);

        memcpy(visual_model_residual_accum_save, visual_model_residual_accum,  sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);

    } else {



        if (last_marg_info) {
            prior_rhs = prior_rhs_save;
            prior_residual = prior_residual_save;
            if (!have_hist) {
                yaw_constraint_residual = yaw_constraint_residual_save;
                yaw_constraint_model_residual_accum = yaw_constraint_model_residual_accum_save;
            }
            if (ESTIMATE_SCALE && fix_scale && !have_hist) {
                init_distance_constraint_residual = init_distance_constraint_residual_save;
                init_distance_constraint_model_residual_accum = init_distance_constraint_model_residual_accum_save;
            }
        }

        memcpy(visual_residuals, visual_residuals_save, sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);
        memcpy(visual_model_residual_accum, visual_model_residual_accum_save, sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);

    }

}

// pindexi != pindexj;
void VisualInertialBase::VisualEvaluateRhs(int pindexi, int pindexj,
                                           double* rhs_idepth,
                                           double* jacobian_pointer,
                                           double* visual_residual_raw_in
#if USE_STEREO
                                           , double* lhs_scale_idepth, bool is_stereo
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                                           , double* lhs_extrinsic_idepth
#endif
#endif
                                          ) {

    int porderi = pindex2order[pindexi];
    int porderj = pindex2order[pindexj];

    bool need_update_pos = pindexi != pindexj;
    double tmp[6];
    if (need_update_pos) {
        MatrixTransposeVectorMultiply<2, 6, 0>
        (jacobian_pointer, 2, 6, visual_residual_raw_in, tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[porderi]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[porderi]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        MatrixTransposeVectorMultiply<2, 6, 0>
        (jacobian_pointer + 6 * 2, 2, 6, visual_residual_raw_in, tmp );
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[porderj]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[porderj]) += Eigen::Map<flower_focus::Vector6d>(tmp);
    }

    MatrixTransposeVectorMultiply<2, 1, 1>
    (jacobian_pointer + 6 * 4, 2, 1, visual_residual_raw_in, rhs_idepth);


#if USE_STEREO
    if (!need_update_pos)ASSERT(is_stereo);

    if (is_stereo) {
#if ESTIMATE_SCALE
        ASSERT(O_SCALE > porderj);
        ASSERT(O_SCALE > porderi);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        ASSERT(O_EXTRINSIC_RIGHT > porderj);
        ASSERT(O_EXTRINSIC_RIGHT > porderi);
#if ESTIMATE_SCALE
        ASSERT(O_EXTRINSIC_RIGHT > O_SCALE);
#endif
#endif

        if (need_update_pos) {
#if ESTIMATE_SCALE
            MatrixTransposeVectorMultiply<2, 1, 0>
            (jacobian_pointer + 6 * 4 + 2, 2, 1, visual_residual_raw_in, tmp );
            rhs_posbias_pointer[order2p_local[O_SCALE]] += tmp[0];
            rhs_posbias_save_pointer[order2p_local[O_SCALE]] += tmp[0];
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO

            MatrixTransposeVectorMultiply<2, 6, 0>
            (jacobian_pointer + 6 * 4 + 4, 2, 6, visual_residual_raw_in, tmp );
            Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT]) += Eigen::Map<flower_focus::Vector6d>(tmp);
            Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[O_EXTRINSIC_RIGHT]) += Eigen::Map<flower_focus::Vector6d>(tmp);

#endif

        }

    }
#endif
}



void VisualInertialBase::ComputeRhs() {


    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
        IMUGNSSPreFactor* imu_pre_factori = it->second;

        const auto& evaluate_idx = imu_pre_factori->evaluate_idx;
        const auto& orders = imu_pre_factori->orders;
        const auto& order_idx = imu_pre_factori->order_idx;
        const auto& parameter_size_sum = imu_pre_factori->parameter_size_sum;
        const auto& parameter_sizes = imu_pre_factori->parameter_sizes;

        Eigen::VectorXd imu_gradient(parameter_size_sum);
        flower_focus::VectorXd imu_rhs(parameter_size_sum);

        imu_pre_factori->factor->EvaluateRhs(imu_pre_factori->parameters.data(), imu_rhs.data(), imu_gradient.data());



        for (int i = 0; i < (int)orders.size(); ++i) {
            int orderi = orders[i];
            int sizei = parameter_sizes[i];

            Eigen::Map<flower_focus::VectorXd>(rhs_posbias_pointer + order2p_local[orderi] + order_idx[i], sizei) += imu_rhs.segment(evaluate_idx[i], sizei);
            Eigen::Map<flower_focus::VectorXd>(rhs_posbias_save_pointer + order2p_local[orderi] + order_idx[i], sizei) += imu_gradient.segment(evaluate_idx[i], sizei);

        }
    }

    int v_index = 0;
    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        VisualFactor* visual_factor = it1->second;
        visual_factor->rhs_idepth[0] = 0;
        auto& pindexs = visual_factor->pindexs;
        for (int i = 1; i < (int)pindexs.size(); ++i) {
            VisualEvaluateRhs(pindexs[0], pindexs[i],
                              visual_factor->rhs_idepth,
                              visual_jacobians + v_index * VRNUM * 2,
                              visual_residuals + v_index * 2
#if USE_STEREO
                              , & (visual_factor->lhs_scale_idepth), visual_factor->is_stereo[i]
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                              , visual_factor->lhs_extrinsic_idepth
#endif
#endif
                             );
            v_index++;
        }
        visual_factor->rhs_idepth_save[0] = visual_factor->rhs_idepth[0];
    }
    ASSERT(v_index == visual_residual_short_num);
    for (auto& it1 : long_feature_factors) {
        it1->rhs_idepth[0] = 0;
        auto& pindexs = it1->pindexs;
        for (int i = 1; i < (int)pindexs.size(); ++i) {
            VisualEvaluateRhs(pindexs[0], pindexs[i],
                              it1->rhs_idepth,
                              visual_jacobians + v_index * VRNUM * 2,
                              visual_residuals + v_index * 2
#if USE_STEREO
                              , & (it1->lhs_scale_idepth), it1->is_stereo[i]
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                              , it1->lhs_extrinsic_idepth
#endif
#endif
                             );
            v_index++;
        }
        it1->rhs_idepth_save[0] = it1->rhs_idepth[0];
    }
    ASSERT(v_index == visual_residual_short_num + visual_residual_long_num);


}


void VisualInertialBase::MargeRhs() {

    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        VisualFactor* visual_factor = it1->second;
        auto lhs_pos_idepth = visual_factor->lhs_pos_idepth;
        auto min_pindex = visual_factor->min_pindex;
        auto lhs_idepth = visual_factor->lhs_idepth;
        auto rhs_idepth = visual_factor->rhs_idepth;
        for (auto& pindexi : visual_factor->pindexs_set) {
            int porderi = pindex2order[pindexi];
            Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[porderi]) -= flower_focus::Vector6d(lhs_pos_idepth + 6 * (pindexi - min_pindex)) * (lhs_idepth[0] * rhs_idepth[0]);
        }
#if USE_STEREO
#if ESTIMATE_SCALE
        rhs_posbias_pointer[order2p_local[O_SCALE]] -= visual_factor->lhs_scale_idepth * lhs_idepth[0] * rhs_idepth[0];
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT]) -= Eigen::Map<flower_focus::Vector6d>(visual_factor->lhs_extrinsic_idepth) * (lhs_idepth[0] * rhs_idepth[0]);
#endif
#endif
    }

    for (int marg_order = 0; marg_order < E_O_MARG_BIAS; ++marg_order) {
        const std::vector<int>& connect_orders = connections_bias[marg_order];
        int lhs_posbias_indexm = marg_order * O_IDEPTHR;
        double* rhs_posbiasm = rhs_posbias_pointer + order2p_local[marg_order];
        for (int j = 0; j < (int)connect_orders.size(); ++j) {
            int orderj = connect_orders[j];
            int sizej = order2size[orderj];
            MatrixTransposeVectorMultiply <  9, Eigen::Dynamic, -1 > (
                lhs_posbias[lhs_posbias_indexm + orderj],  9, sizej,
                rhs_posbiasm,
                rhs_posbias_pointer + order2p_local[orderj]);
        }
    }

    if (short_phase_bias.size()) {
        int short_phase_bias_num = order2size[O_PHASE_BIAS_SHORT];
        int marg_order = O_PHASE_BIAS_SHORT;
        ASSERT((int)short_phase_bias.size() == order2size[O_PHASE_BIAS_SHORT]);
        std::vector<int>connect_orders;
        for (int order = O_MARG_POS; order < O_POSR_NEXT; order++)connect_orders.push_back(order);
        for (int order = O_SCALE; order < O_IDEPTHR; order++)connect_orders.push_back(order);
        int lhs_posbias_indexm = marg_order * O_IDEPTHR;
        double* rhs_posbiasm = rhs_posbias_pointer + order2p_local[marg_order];
        for (int j = 0; j < (int)connect_orders.size(); ++j) {
            int orderj = connect_orders[j];
            int sizej = order2size[orderj];
            MatrixTransposeVectorMultiply <  Eigen::Dynamic, Eigen::Dynamic, -1 > (
                lhs_posbias[lhs_posbias_indexm + orderj],  short_phase_bias_num, sizej,
                rhs_posbiasm,
                rhs_posbias_pointer + order2p_local[orderj]);
        }
    }



    int shift1 = order2p_local[O_MARG_POS];
    int shift2 = order2p_local[O_POSR];
    for (int marg_orderi = O_MARG_POS; marg_orderi < O_IDEPTHR; ++marg_orderi) {
        double* rhs_pos_biasi = rhs_posbias_pointer + order2p_local[marg_orderi];
        int sizei = order2size[marg_orderi];
        int marg_pindexi = order2p_local[marg_orderi] - shift1;
        int outside_pindexi = order2p_local[marg_orderi] - shift2;
        if (marg_orderi >= O_POSR) Eigen::Map<flower_focus::VectorXd>(rhs_outside.data() + outside_pindexi, sizei) = Eigen::Map<flower_focus::VectorXd>(rhs_pos_biasi, sizei);
        else Eigen::Map<flower_focus::VectorXd>(rhs_margpos_pointer + marg_pindexi, sizei) = Eigen::Map<flower_focus::VectorXd>(rhs_pos_biasi, sizei);
    }

    int findex = order2p_local[O_IDEPTHR] - order2p_local[O_POSR];
    for (auto& it : long_feature_factors) {
        double* rhs_idepth = it->rhs_idepth;
        rhs_outside(findex) = rhs_idepth[0];
        findex++;
    }


    for (int marg_orderi = O_POSR; marg_orderi < O_IDEPTHR; ++marg_orderi) {
        int sizei = order2size[marg_orderi];
        int outside_pindexi = order2p_local[marg_orderi] - shift2;
        gradient_outside.segment(outside_pindexi, sizei) = Eigen::Map<flower_focus::VectorXd>(rhs_posbias_save_pointer + order2p_local[marg_orderi], sizei);
    }

    findex = order2p_local[O_IDEPTHR] - shift2;
    for (auto& it : long_feature_factors) {
        double* rhs_idepth_save = it->rhs_idepth_save;
        gradient_outside(findex) = rhs_idepth_save[0];
        findex++;
    }



    if (last_marg_info) {

        int m = last_marg_info->m;

        const auto& keep_block_size = last_marg_info->keep_block_size;
        const auto& keep_block_idx = last_marg_info->keep_block_idx;

        int shift1 = order2p_local[O_POSR];
        for (int i = 0; i < static_cast<int>(keep_block_size.size()); ++i) {
            int sizei = localSize(keep_block_size[i]);
            int idxi = keep_block_idx[i] - m;
            int idxi2 = prior_idx[i];
            rhs_outside.segment(idxi2 - shift1, sizei) += prior_rhs.segment(idxi, sizei);
            gradient_outside.segment(idxi2 - shift1, sizei) += prior_rhs.segment(idxi, sizei);
        }

        if (!have_hist) {
            flower_focus::Vector6d yaw_constraint_rhs;
            yaw_constraint_rhs = yaw_constraint_jacobian.transpose() * yaw_constraint_residual;
            if (!ESTIMATE_R_WI_WC) {
                rhs_outside.segment(order2p_local[O_POS0] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
                gradient_outside.segment(order2p_local[O_POS0] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
            } else {
                rhs_outside.segment(order2p_local[O_POSGLOBAL] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
                gradient_outside.segment(order2p_local[O_POSGLOBAL] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
            }

        }

        if (ESTIMATE_SCALE && fix_scale && !have_hist) {

            for (int pindexi = 0; pindexi < pos_para_count - 1; pindexi++) {
                int orderi = pindex2order[pindexi];
                ASSERT(order2p_local[orderi] - order2p_local[O_POSR] >= 0);
                rhs_outside.segment(order2p_local[orderi] - order2p_local[O_POSR], 3) += init_distance_constraint_residual * init_distance_constraint_jacobian.segment(pindexi * 3, 3);
                gradient_outside.segment(order2p_local[orderi] - order2p_local[O_POSR], 3) +=  init_distance_constraint_residual * init_distance_constraint_jacobian.segment(pindexi * 3, 3);

            }
        }

    }

    int sizen = lhs_outside.rows();
    int sizem = order2p_local[O_POSR] - order2p_local[O_MARG_POS];
    if (sizem > 0) {
        MatrixTransposeVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
            lhs_margpos_outside_pointer,  sizem, sizen,
            rhs_margpos_pointer,
            rhs_outside.data());
    }

}

void VisualInertialBase::MargRhsNew() {


    ASSERT((int)long_feature_factors.size() - (int)deliver_idepth_pointer.size() == num_alone);

    if (prev_vibase) {
        ASSERT(cur_idxs.size());
        for (int i = 0; i < (int)cur_idxs.size(); ++i) {
            int sizei = cur_sizes[i];
            int cur_idxi = cur_idxs[i];
            int prev_idxi = prev_idxs[i];
            rhs_outside.segment(cur_idxi, sizei) += prev_vibase->rhs_outside.segment(prev_idxi, sizei);
            gradient_outside.segment(cur_idxi, sizei) += prev_vibase->gradient_outside.segment(prev_idxi, sizei);
        }
    }

    if (O_BIAS0 != O_POSR) {
        int m = order2p_local[O_BIAS0] - order2p_local[O_POSR];
        int n = lhs_outside.cols() - m;
        flower_focus::MatrixXd Amn = lhs_outside.block(0, m, m, n);
        MatrixTransposeVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
            Amn.data(),  m, n,
            rhs_outside.data() + 0,
            rhs_outside.data() + m);
    }

    if (num_alone) {
        //marginalize idepth alone
        int k = order2p_local[O_BIAS0] - order2p_local[O_POSR];
        int m = num_alone;
        int n = lhs_outside.cols() - m - k;
        flower_focus::MatrixXd Amn = lhs_outside.block(k, n + k, n, m).transpose();
        MatrixTransposeVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
            Amn.data(),  m, n,
            rhs_outside.data() + n + k,
            rhs_outside.data() + k);

    }

    if (deliver_idepth_pointer.size()) {

        int size = lhs_outside.cols() - num_alone - (order2p_local[O_BIAS0] - order2p_local[O_POSR]);
        int idx = order2p_local[O_BIAS0] - order2p_local[O_POSR];
        Eigen::VectorXd b = rhs_outside.segment(idx, size);
        Eigen::VectorXd b2 = b;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            vector_update(b, info.m, b2, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }

        rhs_outside.segment(idx, size) = b2;
        Eigen::VectorXd c = gradient_outside.segment(idx, size);
        Eigen::VectorXd c2 = c;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            vector_update(c, info.m, c2, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }

        gradient_outside.segment(idx, size) = c2;

    }

    {
        //marginalize pos0 bias0
        int k = order2p_local[O_BIAS0] - order2p_local[O_POSR];
        int m = 15 + phase_bias_alone_num;
        int n = lhs_outside.cols() - num_alone - 15 - k - phase_bias_alone_num;
        flower_focus::MatrixXd Amn = lhs_outside.block(k, k + m, m, n);
        MatrixTransposeVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
            Amn.data(),  m, n,
            rhs_outside.data() + k,
            rhs_outside.data() + k + m);
    }
}


void VisualInertialBase::EvaluateNonlinearQuantity() {
    if (history_flag) {
        nonlinear_quantity = 0;

        if (ESTIMATE_SCALE && last_marg_info && fix_scale && !have_hist) {
            double a = init_distance_constraint_residual + init_distance_constraint_model_residual_accum - init_distance_constraint_residual_linerized;
            nonlinear_quantity += a * a;
        }
        ASSERT(!isnan(nonlinear_quantity));
        if (last_marg_info && !have_hist)
            nonlinear_quantity += (yaw_constraint_residual + yaw_constraint_model_residual_accum - yaw_constraint_residual_linerized).squaredNorm();
        ASSERT(!isnan(nonlinear_quantity));

        nonlinear_quantity += (Eigen::Map<flower_focus::VectorXd>(visual_residuals, (visual_residual_long_num + visual_residual_short_num) * 2) +
                               Eigen::Map<flower_focus::VectorXd>(visual_model_residual_accum, (visual_residual_long_num + visual_residual_short_num) * 2) -
                               Eigen::Map<flower_focus::VectorXd>(visual_residuals_linerized, (visual_residual_long_num + visual_residual_short_num) * 2)).squaredNorm();

        ASSERT(!isnan(nonlinear_quantity));


        for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
            IMUGNSSPreFactor* imu_pre_factori = it->second;
            nonlinear_quantity += imu_pre_factori->factor->IMUGNSS_info->nonlinear_quantity;
        }
        ASSERT(!isnan(nonlinear_quantity));
    } else
        nonlinear_quantity = NLQ_THRESHOLD + 100;

    // LOG_OUT<<nonlinear_quantity<<std::endl;
    if (!prev_vibase) nonlinear_quantity_deliver_accum = nonlinear_quantity;
    else nonlinear_quantity_deliver_accum = (prev_vibase->nonlinear_quantity_deliver_accum > NLQ_THRESHOLD ? NLQ_THRESHOLD + 10 : nonlinear_quantity);

}
