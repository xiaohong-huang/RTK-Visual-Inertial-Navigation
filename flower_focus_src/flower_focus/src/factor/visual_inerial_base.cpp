#include "visual_inerial_base.h"
#include "../utility/utility.h"
#include "operation.h"
#include "gnss_imu_info.h"
#include "../utility/tic_toc.h"




#define YAWSQRTINFO 2e2
#define SQRT_INFO_SCALE (100000)

VisualInertialBase* VI_info_pointer;
VisualInertialBase* VI_info_pointer_head;
VisualInertialBase* VI_info_pointer_tail;

void UpdateTrustRegion(double* pointer, int size, double mu) {
    for (int i = 0; i < size; ++i)
        pointer[i * size + i] +=  mu;
}


VisualInertialBase::~VisualInertialBase() {
    for (auto it = idepth_map_factors_short.begin(); it != idepth_map_factors_short.end(); ++it)
        delete it->second;
    for (auto& it : long_feature_factors)
        delete it;
    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it)
        delete it->second;
    delete cauchy_loss_function;
    if (memory)
        delete memory;

}

void VisualInertialBase::RemoveFeature(int feature_id) {
    bool need_reset = false;
    if (idepth_map_factors_long.find(feature_id) != idepth_map_factors_long.end()) {
        need_reset = true;
        auto it = idepth_map_factors_long[feature_id];
        if (deliver_idepth_pointer.find(it->p_idepth) != deliver_idepth_pointer.end())deliver_idepth_pointer.erase(it->p_idepth);
        idepth_map_factors_long.erase(feature_id);
    }
    if (idepth_map_factors_short.find(feature_id) != idepth_map_factors_short.end()) {
        need_reset = true;
        idepth_map_factors_short.erase(feature_id);
    }
    if (!need_reset)return;

    ResetInit();
}

void VisualInertialBase::AddLastMargeInfo(MarginalizationInfo* last_marg_info_) {
    last_marg_info = last_marg_info_;
}


void VisualInertialBase::ConstructPriorIdx() {

#if USE_ASSERT
    {
        const auto& keep_block_size = last_marg_info->keep_block_size;
        const auto& keep_block_addr = last_marg_info->keep_block_addr;
        std::set<double*>feature_pointers;
        for (auto& it : long_feature_factors)
            feature_pointers.insert(it->p_idepth);

        for (int i = 0; i < (int)keep_block_size.size(); ++i) {
            if (keep_block_size[i] == 9) {
                bool found = false;
                for (auto it = bias_index2pointers.begin(); it != bias_index2pointers.end(); ++it) {
                    if (it->second == keep_block_addr[i])found = true;
                }
                ASSERT(found);
            }
            if (keep_block_size[i] == 7) {
                bool found = false;
                for (auto it = pos_index2pointers.begin(); it != pos_index2pointers.end(); ++it) {
                    if (it->second == keep_block_addr[i])found = true;
                }
                if (ESTIMATE_EXTRINSIC)
                    ASSERT(found || keep_block_addr[i] == global_pos_pointer || keep_block_addr[i] == extrinsic_pointer || keep_block_addr[i] == rwgw_pointer
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                           || keep_block_addr[i] == extrinsic2_pointer
#endif
                          );
                else
                    ASSERT(found || keep_block_addr[i] == global_pos_pointer || keep_block_addr[i] == rwgw_pointer
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                           || keep_block_addr[i] == extrinsic2_pointer
#endif
                    );
            }
            if (keep_block_size[i] == 1) {
                if (ESTIMATE_TD)
                    ASSERT(feature_pointers.find(keep_block_addr[i]) != feature_pointers.end() || keep_block_addr[i] == scale_pointer || phase_bias_pointers_set.find(keep_block_addr[i]) != phase_bias_pointers_set.end() || keep_block_addr[i] == td_pointer);
                else
                    ASSERT(feature_pointers.find(keep_block_addr[i]) != feature_pointers.end() || keep_block_addr[i] == scale_pointer || phase_bias_pointers_set.find(keep_block_addr[i]) != phase_bias_pointers_set.end());
            }
            if (ESTIMATE_ACC_SCALE || ESTIMATE_PBG) {
                if (keep_block_size[i] == 3)
                    ASSERT(keep_block_addr[i] == acc_scale_pointer || keep_block_addr[i] == pbg_pointer);
            }

        }
    }
#endif

    const auto& keep_block_size = last_marg_info->keep_block_size;
    const auto& keep_block_addr = last_marg_info->keep_block_addr;
    prior_idx.clear();
    for (int i = 0; i < static_cast<int>(keep_block_size.size()); ++i) {
        double* addri = keep_block_addr[i];
        if (long_phase_bias_idx_map.find(addri) != long_phase_bias_idx_map.end()) {
            prior_idx.push_back(order2p_local[O_PHASE_BIAS_LONG] + long_phase_bias_idx_map[addri]);
            continue;
        }
        ASSERT(pointer2orders.find(addri) != pointer2orders.end());
        int idxi2 = order2p_local[pointer2orders[addri]];
        prior_idx.push_back(idxi2);
    }
}

void VisualInertialBase::EvaluateLastMargInfo() {


    int m = last_marg_info->m;
    const auto& keep_block_size = last_marg_info->keep_block_size;
    const auto& keep_block_idx = last_marg_info->keep_block_idx;

    Eigen::MatrixXd prior_lhs;

    if (USE_REPARAM)
        last_marg_info->UpdateMatix(last_marg_info->FTAF, prior_lhs);
    else
        prior_lhs = last_marg_info->A;



    int shift1 = order2p_local[O_POSR];
    for (int i = 0; i < static_cast<int>(keep_block_size.size()); ++i) {

        int sizei = localSize(keep_block_size[i]);
        int idxi = keep_block_idx[i] - m;
        int idxi2 = prior_idx[i];
        rhs_outside.segment(idxi2 - shift1, sizei) += prior_rhs.segment(idxi, sizei);
        gradient_outside.segment(idxi2 - shift1, sizei) += prior_rhs.segment(idxi, sizei);

        for (int j = 0; j < static_cast<int>(keep_block_size.size()); ++j) {
            int sizej = localSize(keep_block_size[j]);
            int idxj = keep_block_idx[j] - m;
            int idxj2 = prior_idx[j];
            if (idxj2 >= idxi2)
                lhs_outside.block(idxi2 - shift1, idxj2 - shift1, sizei, sizej) += prior_lhs.block(idxi, idxj, sizei, sizej);
        }
    }

    if (!have_hist) {
        double* parameter0;
        if (!ESTIMATE_R_WI_WC)
            parameter0 = order2pointers[O_POS0];
        else
            parameter0 = order2pointers[O_POSGLOBAL];

        Eigen::Vector3d Pwb(parameter0[0], parameter0[1], parameter0[2]);
        Eigen::Quaterniond Qwb(parameter0[6], parameter0[3], parameter0[4], parameter0[5]);
        Eigen::Matrix3d Rwb = Qwb.toRotationMatrix();

        Eigen::Vector3d Magw1 = Rwb * InitMag;
        double mag1_2_norm = Magw1.segment(0, 2).norm();
        yaw_constraint_jacobian.setZero();
        yaw_constraint_jacobian.segment(3, 3) = (Rwb * Utility::skewSymmetric(-InitMag) / mag1_2_norm).block(0, 0, 1, 3) * (USE_GNSS ? 0 : YAWSQRTINFO);

        flower_focus::Matrix6d yaw_constraint_lhs;
        flower_focus::Vector6d yaw_constraint_rhs;

        yaw_constraint_lhs = yaw_constraint_jacobian.transpose() * yaw_constraint_jacobian;
        yaw_constraint_rhs = yaw_constraint_jacobian.transpose() * yaw_constraint_residual;

        if (!ESTIMATE_R_WI_WC) {
            lhs_outside.block(order2p_local[O_POS0] - order2p_local[O_POSR], order2p_local[O_POS0] - order2p_local[O_POSR], 6, 6) += yaw_constraint_lhs;
            rhs_outside.segment(order2p_local[O_POS0] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
            gradient_outside.segment(order2p_local[O_POS0] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
        } else {
            lhs_outside.block(order2p_local[O_POSGLOBAL] - order2p_local[O_POSR], order2p_local[O_POSGLOBAL] - order2p_local[O_POSR], 6, 6) += yaw_constraint_lhs;
            rhs_outside.segment(order2p_local[O_POSGLOBAL] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
            gradient_outside.segment(order2p_local[O_POSGLOBAL] - order2p_local[O_POSR], 6) += yaw_constraint_rhs;
        }


    }

    if (ESTIMATE_SCALE && fix_scale && !have_hist) {
        init_distance_constraint_jacobian.setZero();
        int order0 = pindex2order[0];
        Eigen::Vector3d P0(order2pointers[order0]);
        for (int pindex = 1; pindex < pos_para_count - 1; pindex++) {
            int orderi = pindex2order[pindex];
            Eigen::Vector3d Pk(order2pointers[orderi]);
            init_distance_constraint_jacobian.segment(pindex * 3, 3) += (Pk - P0) / (Pk - P0).norm() * SQRT_INFO_SCALE;
            init_distance_constraint_jacobian.segment(0 * 3, 3) += -(Pk - P0) / (Pk - P0).norm() * SQRT_INFO_SCALE;
        }

        for (int pindexi = 0; pindexi < pos_para_count - 1; pindexi++) {
            int orderi = pindex2order[pindexi];
            ASSERT(order2p_local[orderi] - order2p_local[O_POSR] >= 0);
            rhs_outside.segment(order2p_local[orderi] - order2p_local[O_POSR], 3) += init_distance_constraint_residual * init_distance_constraint_jacobian.segment(pindexi * 3, 3);
            gradient_outside.segment(order2p_local[orderi] - order2p_local[O_POSR], 3) +=  init_distance_constraint_residual * init_distance_constraint_jacobian.segment(pindexi * 3, 3);
            for (int pindexj = 0; pindexj < pos_para_count - 1; pindexj++) {
                int orderj = pindex2order[pindexj];
                if (orderi <= orderj) {
                    ASSERT(order2p_local[orderj] - order2p_local[O_POSR] >= 0);
                    lhs_outside.block(order2p_local[orderi] - order2p_local[O_POSR], order2p_local[orderj] - order2p_local[O_POSR], 3, 3) +=
                        init_distance_constraint_jacobian.segment(pindexi * 3, 3) * init_distance_constraint_jacobian.segment(pindexj * 3, 3).transpose();
                }
            }
        }
    }
}

void VisualInertialBase::EvaluatePriorAlpha(double* alpha1) {
    const auto& keep_block_size = last_marg_info->keep_block_size;
    const auto& keep_block_idx = last_marg_info->keep_block_idx;
    int n = last_marg_info->n;
    int m = last_marg_info->m;
    flower_focus::VectorXd gradient_prior(n);
    gradient_prior.setZero();
    int shift1 = order2p_local[O_POSR];
    for (int i = 0; i < static_cast<int>(keep_block_size.size()); ++i) {

        int sizei = localSize(keep_block_size[i]);
        int idxi = keep_block_idx[i] - m;
        int idxi2 = prior_idx[i];
        gradient_prior.segment(idxi, sizei) = Eigen::Map<flower_focus::VectorXd>(gradient_outside.data() + idxi2 - shift1, sizei);
    }

    if (USE_REPARAM) {
        Eigen::VectorXd tmp;
        last_marg_info->UpdateVector(gradient_prior, tmp, false);
        alpha1[0] += (last_marg_info->linearized_jacobiansF * tmp).squaredNorm();
    } else
        alpha1[0] += (last_marg_info->linearized_jacobians * gradient_prior).squaredNorm();
    if (!have_hist) {
        flower_focus::VectorXd model_residuals;
        if (!ESTIMATE_R_WI_WC)
            model_residuals = yaw_constraint_jacobian * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[O_POS0]);
        else
            model_residuals = yaw_constraint_jacobian * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[O_POSGLOBAL]);

        alpha1[0] += model_residuals.squaredNorm();
    }
    if (ESTIMATE_SCALE && fix_scale && !have_hist) {
        for (int pindexi = 0; pindexi < pos_para_count - 1; pindexi++) {
            int orderi = pindex2order[pindexi];
            ASSERT(order2p_local[orderi] - order2p_local[O_POSR] >= 0);
            alpha1[0] += (init_distance_constraint_jacobian.segment(pindexi * 3, 3).transpose() * gradient_outside.segment(order2p_local[orderi] - order2p_local[O_POSR], 3)).squaredNorm();
        }
    }

}


void VisualInertialBase::EvaluatePriorModelCostChange(double* model_cost_change) {
    const auto& keep_block_size = last_marg_info->keep_block_size;
    const auto& keep_block_idx = last_marg_info->keep_block_idx;
    int n = last_marg_info->n;
    int m = last_marg_info->m;
    flower_focus::VectorXd inc_prior(n);
    inc_prior.setZero();
    int shift1 = order2p_local[O_POSR];
    for (int i = 0; i < static_cast<int>(keep_block_size.size()); ++i) {

        int sizei = localSize(keep_block_size[i]);
        int idxi = keep_block_idx[i] - m;
        int idxi2 = prior_idx[i];
        inc_prior.segment(idxi, sizei) = Eigen::Map<flower_focus::VectorXd>(inc_outside_pointer + idxi2 - shift1, sizei);
    }


    ASSERT((prior_residual_save - prior_residual).norm() == 0);
    flower_focus::VectorXd model_residuals;
    if (USE_REPARAM) {
        Eigen::VectorXd tmp;
        last_marg_info->UpdateVector(inc_prior, tmp, false);
        model_residuals = last_marg_info->linearized_jacobiansF * tmp;
    } else
        model_residuals = last_marg_info->linearized_jacobians * inc_prior;

    model_cost_change[0] += model_residuals.transpose() * ( prior_residual_save - model_residuals / 2);

    if (!have_hist) {
        flower_focus::VectorXd model_residuals;
        if (!ESTIMATE_R_WI_WC)
            model_residuals = yaw_constraint_jacobian * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[O_POS0]);
        else
            model_residuals = yaw_constraint_jacobian * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[O_POSGLOBAL]);

        model_cost_change[0] += model_residuals.transpose() * ( yaw_constraint_residual_save - model_residuals / 2);
        yaw_constraint_model_residual_accum += model_residuals;
    }

    if (ESTIMATE_SCALE && fix_scale && !have_hist) {
        double model_residuals = 0;
        for (int pindexi = 0; pindexi < pos_para_count - 1; pindexi++) {
            int orderi = pindex2order[pindexi];
            ASSERT(order2p_local[orderi] - order2p_local[O_POSR] >= 0);
            model_residuals += (init_distance_constraint_jacobian.segment(pindexi * 3, 3).transpose() * flower_focus::Vector3d(inc_outside_pointer + order2p_local[orderi] - shift1))(0);
        }

        model_cost_change[0] += model_residuals * (init_distance_constraint_residual_save - model_residuals / 2);
        init_distance_constraint_model_residual_accum += model_residuals;
    }

}

void VisualInertialBase::EvaluatePriorCost(double* cost) {

    if (USE_REPARAM)
        last_marg_info->compute_residual_rhs(prior_residual, prior_rhs);
    else {
        last_marg_info->resetLinerizationPoint();
        prior_rhs = last_marg_info->b;
        prior_residual = last_marg_info->linearized_residuals;
    }


    cost[0] += 0.5 * prior_residual.squaredNorm();
    ASSERT(cost[0] > 0);


    if (!have_hist) {
        double* parameter0;
        if (!ESTIMATE_R_WI_WC)
            parameter0 = order2pointers[O_POS0];
        else
            parameter0 = order2pointers[O_POSGLOBAL];
        Eigen::Vector3d Pwb(parameter0[0], parameter0[1], parameter0[2]);
        Eigen::Quaterniond Qwb(parameter0[6], parameter0[3], parameter0[4], parameter0[5]);
        Eigen::Matrix3d Rwb = Qwb.toRotationMatrix();
        Eigen::Vector3d Magw1 = Rwb * InitMag;
        double mag1_2_norm = Magw1.segment(0, 2).norm();
        yaw_constraint_residual(0) = Magw1(0) / mag1_2_norm * (USE_GNSS ? 0 : YAWSQRTINFO);
        cost[0] +=  0.5 * yaw_constraint_residual.squaredNorm();
    }

    if (ESTIMATE_SCALE && fix_scale && !have_hist) {
        double sum_distance = 0;
        for (int pindex = 1; pindex < pos_para_count - 1; pindex++) {
            int orderi = pindex2order[pindex];
            int order0 = pindex2order[0];
            Eigen::Vector3d Pk(order2pointers[orderi]);
            Eigen::Vector3d P0(order2pointers[order0]);
            sum_distance += (Pk - P0).norm();
        }
        init_distance_constraint_residual = (sum_distance - distance0) * SQRT_INFO_SCALE;

        cost[0] += 0.5 * init_distance_constraint_residual * init_distance_constraint_residual;
    }

}




void VisualInertialBase::AddVisualFactorShort(int feature_id,
                                              double* p0, double* p1, double* p2,
                                              int pindexi, int pindexj,
                                              Eigen::Vector3d pts_i, Eigen::Vector3d pts_j
#if USE_STEREO
                                              , bool is_stereo_
#endif
                                             ) {
    if (!idepth_map_factors_short[feature_id]) {
        idepth_map_factors_short[feature_id] = new VisualFactor();
        idepth_map_factors_short[feature_id]->p_idepth = p2;
        idepth_map_factors_short[feature_id]->ptss.push_back(pts_i);
        idepth_map_factors_short[feature_id]->pindexs.push_back(pindexi);
        idepth_map_factors_short[feature_id]->pindexs_set.insert(pindexi);
#if USE_STEREO
        idepth_map_factors_short[feature_id]->is_stereo.push_back(0);
#endif
    }
    idepth_map_factors_short[feature_id]->ptss.push_back(pts_j);
    idepth_map_factors_short[feature_id]->pindexs.push_back(pindexj);
    idepth_map_factors_short[feature_id]->pindexs_set.insert(pindexj);
#if USE_STEREO
    idepth_map_factors_short[feature_id]->is_stereo.push_back(is_stereo_);
#endif

#if USE_ASSERT
    if (pos_index2pointers[pindexi])
        ASSERT(pos_index2pointers[pindexi] == p0);
    if (pos_index2pointers[pindexj])
        ASSERT(pos_index2pointers[pindexj] == p1);
#endif
    pos_index2pointers[pindexi] = p0;
    pos_index2pointers[pindexj] = p1;
    pos_pointer2indexs[p0] = pindexi;
    pos_pointer2indexs[p1] = pindexj;
}




void VisualInertialBase::AddVisualFactorLong(int feature_id,
                                             double* p0, double* p1, double* p2,
                                             int pindexi, int pindexj,
                                             Eigen::Vector3d pts_i, Eigen::Vector3d pts_j
#if USE_STEREO
                                             , bool is_stereo_
#endif
                                            ) {
    if (!idepth_map_factors_long[feature_id]) {
        idepth_map_factors_long[feature_id] = new VisualFactor();
        idepth_map_factors_long[feature_id]->p_idepth = p2;
        idepth_map_factors_long[feature_id]->ptss.push_back(pts_i);
        idepth_map_factors_long[feature_id]->pindexs.push_back(pindexi);
        idepth_map_factors_long[feature_id]->pindexs_set.insert(pindexi);
#if USE_STEREO
        idepth_map_factors_long[feature_id]->is_stereo.push_back(0);
#endif

    }
    idepth_map_factors_long[feature_id]->ptss.push_back(pts_j);
    idepth_map_factors_long[feature_id]->pindexs.push_back(pindexj);
    idepth_map_factors_long[feature_id]->pindexs_set.insert(pindexj);
#if USE_STEREO
    idepth_map_factors_long[feature_id]->is_stereo.push_back(is_stereo_);
#endif

#if USE_ASSERT
    if (pos_index2pointers[pindexi])
        ASSERT(pos_index2pointers[pindexi] == p0);
    if (pos_index2pointers[pindexj])
        ASSERT(pos_index2pointers[pindexj] == p1);
#endif

    pos_index2pointers[pindexi] = p0;
    pos_index2pointers[pindexj] = p1;
    pos_pointer2indexs[p0] = pindexi;
    pos_pointer2indexs[p1] = pindexj;

}


void VisualInertialBase::AddIMUFactor(IMUGNSSPreFactor* imu_pre_factor) {
    imu_map_factors[imu_pre_factor->pindexi] = imu_pre_factor;
    int j_idx = 0;
    if (ESTIMATE_SCALE) {
        scale_pointer = imu_pre_factor->parameters[j_idx]; j_idx++;
    }
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC) {
        global_pos_pointer = imu_pre_factor->parameters[j_idx]; j_idx++;
    }
    if (ESTIMATE_EXTRINSIC) {
        extrinsic_pointer = imu_pre_factor->parameters[j_idx]; j_idx++;
    }
#if USE_ASSERT
    if (pos_index2pointers[imu_pre_factor->pindexi])
        ASSERT(pos_index2pointers[imu_pre_factor->pindexi] == imu_pre_factor->parameters[j_idx]);
    if (pos_index2pointers[imu_pre_factor->pindexj])
        ASSERT(pos_index2pointers[imu_pre_factor->pindexj] == imu_pre_factor->parameters[j_idx + 2]);
    ASSERT(imu_pre_factor->pindexi < SWF_SIZE_IN + 1);
    ASSERT(imu_pre_factor->pindexj < SWF_SIZE_IN + 1);
#endif
    pos_index2pointers[imu_pre_factor->pindexi] = imu_pre_factor->parameters[j_idx]; j_idx++;
    bias_index2pointers[imu_pre_factor->pindexi] = imu_pre_factor->parameters[j_idx]; j_idx++;
    pos_index2pointers[imu_pre_factor->pindexj] = imu_pre_factor->parameters[j_idx]; j_idx++;
    bias_index2pointers[imu_pre_factor->pindexj] = imu_pre_factor->parameters[j_idx]; j_idx++;
    pos_pointer2indexs[pos_index2pointers[imu_pre_factor->pindexi]] = imu_pre_factor->pindexi;
    pos_pointer2indexs[pos_index2pointers[imu_pre_factor->pindexj]] = imu_pre_factor->pindexj;

    if (ESTIMATE_TD) {
        td_pointer = imu_pre_factor->parameters[j_idx]; j_idx++;
    }
    if (ESTIMATE_ACC_SCALE) {
        acc_scale_pointer = imu_pre_factor->parameters[j_idx]; j_idx++;
    }
    if (ESTIMATE_PBG) {
        pbg_pointer = imu_pre_factor->parameters[j_idx]; j_idx++;
    }


    for (int i = j_idx; i < (int)imu_pre_factor->parameters.size(); i++)
        phase_bias_pointers_set.insert(imu_pre_factor->parameters[i]);

}

int VisualInertialBase::GetBIndex(int index) {
    ASSERT(index < bias_para_count);
    return index;
}


int VisualInertialBase::GetPIndex(int index) {
    ASSERT(index < pos_para_count);
    return index + bias_para_count;
}





void VisualInertialBase::ConstructOrderingMap() {

    std::vector<int> order2index;
    //add margbias order
    std::vector<int> ordering_tmp;
    for (int index = 1; index < marg_count + 1; index++)
        ordering_tmp.push_back(GetBIndex(index));
    order2size.clear();
    while (ordering_tmp.size()) {
        int count = 0;
        for (int i = 0; i < (int)ordering_tmp.size(); ++i) {
            if (std::find(order2index.begin(), order2index.end(), ordering_tmp[i]) == order2index.end()) {
                if (count % 2 == 0) {
                    order2index.push_back(ordering_tmp[i]);
                    order2size.push_back(9);
                    ordering_tmp.erase(ordering_tmp.begin() + i);
                    i--;
                }
                count++;
            }
        }
    }

    //add marg pos order
    for (int index = 1 + outside_threshold; index < marg_count + 1; index++) {
        order2index.push_back(GetPIndex(index));
        order2size.push_back(6);
    }

    if (short_phase_bias.size()) {
        order2index.push_back(-1);
        order2size.push_back(short_phase_bias.size());
// LOG_OUT<<"short_phase_bias.size():"<<short_phase_bias.size()<<std::endl;
    }

    for (int index = 1; index < outside_threshold + 1; index++) {
        order2index.push_back(GetPIndex(index));
        order2size.push_back(6);
    }
    //add outside pos and bias
    // ASSERT((int)order2index.size() == marg_count * 2);

    order2index.push_back(GetBIndex(0)); order2size.push_back(9);
    order2index.push_back(GetPIndex(0)); order2size.push_back(6);
    if (long_phase_bias.size()) {
        order2index.push_back(-1);
        order2size.push_back(long_phase_bias.size());
    }

    order2index.push_back(GetPIndex(marg_count + 1)); order2size.push_back(6);
    order2index.push_back(GetBIndex(marg_count + 1)); order2size.push_back(9);


    for (int index = bias_para_count; index < pos_para_count; index++) {
        order2index.push_back(GetPIndex(index));
        order2size.push_back(6);
    }
    if (ESTIMATE_SCALE)
        order2size.push_back(1);//HXH
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
        order2size.push_back(6);
    if (ESTIMATE_TD)
        order2size.push_back(1);
    if (ESTIMATE_EXTRINSIC)
        order2size.push_back(6);
    if (ESTIMATE_ACC_SCALE)
        order2size.push_back(3);
    if (ESTIMATE_EXTRINSIC2 && USE_STEREO)
        order2size.push_back(6);
    if (ESTIMATE_PBG)
        order2size.push_back(3);
    for (int i = 0; i < idepth_long_count; ++i)
        order2size.push_back(1);









    pindex2order.resize(pos_para_count);
    bindex2order.resize(bias_para_count);
    //index2order,order2pointers
    order2pointers.clear();
    for (int order = 0; order < (int)order2index.size(); ++order) {

        for (int index = 0; index < bias_para_count; index++) {
            if (order2index[order] == index) {
                bindex2order[index] = order;
                order2pointers[order] = bias_index2pointers[index];
                break;
            }
        }
        for (int index = 0; index < pos_para_count; index++) {
            if (order2index[order] == index + bias_para_count) {
                pindex2order[index] = order;
                order2pointers[order] = pos_index2pointers[index];
                break;
            }
        }
    }
    if (ESTIMATE_SCALE)
        order2pointers[O_SCALE] = scale_pointer;
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
        order2pointers[O_POSGLOBAL] = global_pos_pointer;
    if (ESTIMATE_TD)
        order2pointers[O_TD] = td_pointer;
    if (ESTIMATE_EXTRINSIC)
        order2pointers[O_EXTRINSIC] = extrinsic_pointer;
    if (ESTIMATE_ACC_SCALE)
        order2pointers[O_ACC_S] = acc_scale_pointer;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    order2pointers[O_EXTRINSIC_RIGHT] = extrinsic2_pointer;
#endif
    if (ESTIMATE_PBG)
        order2pointers[O_PBG] = pbg_pointer;
    if (short_phase_bias.size())
        order2pointers[O_PHASE_BIAS_SHORT] = 0;
    if (long_phase_bias.size())
        order2pointers[O_PHASE_BIAS_LONG] = 0;

    int order = O_IDEPTHR;
    for (auto& it : long_feature_factors) {
        order2pointers[order] = it->p_idepth;
        ++order;
    }
    pointer2orders.clear();
    for (auto it = order2pointers.begin(); it != order2pointers.end(); ++it)
        pointer2orders[it->second] = it->first;




    order2p_local.resize(order2size.size() + 1);
    order2p_local[0] = 0;
    for (int i = 0; i < (int)order2size.size(); ++i)
        order2p_local[i + 1] = order2p_local[i] + order2size[i];

    order2p_global.resize(order2size.size() + 1);
    order2p_global[0] = 0;
    for (int i = 0; i < (int)order2size.size(); ++i)
        order2p_global[i + 1] = order2p_global[i] + globalSize(order2size[i]);
    pointer2p_local.clear();
    for (int i = 0; i < (int)order2size.size(); ++i)
        pointer2p_local[order2pointers[i]] = order2p_local[i];

}


void VisualInertialBase::ConstructBiasConnections() {
    connections_bias.clear();
    connections_bias.resize(marg_count);
    for (int index = 1; index < marg_count + 1; index++) {
        connections_bias[bindex2order[index]].push_back(pindex2order[index]);//posZ
        connections_bias[bindex2order[index]].push_back(bindex2order[index - 1]); //bias
        connections_bias[bindex2order[index]].push_back(pindex2order[index - 1]); //pos
        connections_bias[bindex2order[index]].push_back(bindex2order[index + 1]); //bias
        connections_bias[bindex2order[index]].push_back(pindex2order[index + 1]); //pos
        if (ESTIMATE_SCALE)
            connections_bias[bindex2order[index]].push_back(O_SCALE); //pos
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
            connections_bias[bindex2order[index]].push_back(O_POSGLOBAL); //pos
        if (ESTIMATE_TD) {
            connections_bias[bindex2order[index]].push_back(O_TD); //pos
        }
        if (ESTIMATE_EXTRINSIC) {
            connections_bias[bindex2order[index]].push_back(O_EXTRINSIC); //pos
        }
        if (ESTIMATE_ACC_SCALE) {
            connections_bias[bindex2order[index]].push_back(O_ACC_S); //pos
        }
        if (ESTIMATE_PBG) {
            connections_bias[bindex2order[index]].push_back(O_PBG); //pos
        }
        if (short_phase_bias.size())connections_bias[bindex2order[index]].push_back(O_PHASE_BIAS_SHORT); //pos
        if (long_phase_bias.size())connections_bias[bindex2order[index]].push_back(O_PHASE_BIAS_LONG); //pos
    }

    for (int i = 0; i < marg_count; ++i) { //marg_index
        const std::vector<int>& connect_orders = connections_bias[i];
        for (int j = 0; j < (int)connect_orders.size(); ++j) {
            int indexj = connect_orders[j]; //outside_index
            for (int k = 0; k < (int)connect_orders.size(); k++) {
                int indexk = connect_orders[k];
                if (indexk != indexj && indexj < marg_count && indexj > i &&
                        std::find(connections_bias[indexj].begin(), connections_bias[indexj].end(), indexk) == connections_bias[indexj].end())
                    connections_bias[indexj].push_back(indexk);
            }
        }
    }

    for (int i = 0; i < marg_count; ++i) {
        const std::vector<int>& connect_orders = connections_bias[i];
        for (int j = 0; j < (int)connect_orders.size(); ++j) {
            if (connect_orders[j] <= i)
                connections_bias[i].erase(connections_bias[i].begin() + j);
        }
    }

    for (int i = 0; i < marg_count; ++i)
        std::sort(connections_bias[i].begin(), connections_bias[i].end());

#if USE_ASSERT
    for (int i = 0; i < marg_count; ++i) {
        if (i < (marg_count + 1) / 2) {
            int p_idx = 5;
            if (ESTIMATE_TD)p_idx++;
            if (ESTIMATE_EXTRINSIC)p_idx++;
            if (ESTIMATE_ACC_SCALE)p_idx ++;
            if (ESTIMATE_PBG)p_idx ++;
            if (ESTIMATE_SCALE)p_idx++;
            if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)p_idx++;
            if (short_phase_bias.size())p_idx++;
            if (long_phase_bias.size())p_idx++;

            ASSERT((int)connections_bias[i].size() == p_idx);
        }

    }
    if (connections_bias.size()) {
        int p_idx = bias_para_count + 2;
        if (ESTIMATE_TD)p_idx++;
        if (ESTIMATE_EXTRINSIC)p_idx++;
        if (ESTIMATE_ACC_SCALE)p_idx ++;
        if (ESTIMATE_PBG)p_idx ++;
        if (ESTIMATE_SCALE)p_idx++;
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)p_idx++;
        if (short_phase_bias.size())p_idx++;
        if (long_phase_bias.size())p_idx++;
        ASSERT((int)connections_bias[connections_bias.size() - 1].size() == p_idx);
    }
    for (int i = 0; i < marg_count; ++i) {
        for (int j = 1; j < (int)connections_bias[i].size(); ++j)
            ASSERT(connections_bias[i][j] > connections_bias[i][j - 1]);
        ASSERT(connections_bias[i][1] > connections_bias[i][0]);
    }
#endif

}







void VisualInertialBase::PrepareMatrix() {


    visual_residual_short_num = 0;
    for (auto it = idepth_map_factors_short.begin(); it != idepth_map_factors_short.end(); ++it)
        visual_residual_short_num += it->second->pindexs.size() - 1;

    visual_residual_long_num = 0;
    for (auto& it : long_feature_factors)
        visual_residual_long_num += it->pindexs.size() - 1;


    int mpos_size = order2p_local[O_POSR] - order2p_local[O_MARG_POS];
    int r_size = order2p_local[O_FULL] - order2p_local[O_POSR];
    visual_lhs_count = 0;
    for (auto it = idepth_map_factors_short.begin(); it != idepth_map_factors_short.end(); ++it) {
        it->second->min_pindex = *std::min_element(it->second->pindexs_set.begin(), it->second->pindexs_set.end());
        int max_pindex = *std::max_element(it->second->pindexs_set.begin(), it->second->pindexs_set.end());
        visual_lhs_count += max_pindex - it->second->min_pindex + 1;
    }
    for (auto& it : long_feature_factors) {
        it->min_pindex = *std::min_element(it->pindexs_set.begin(), it->pindexs_set.end());
        int max_pindex = *std::max_element(it->pindexs_set.begin(), it->pindexs_set.end());
        visual_lhs_count += max_pindex - it->min_pindex + 1;
    }

    uint32_t cnt1 = 0;
    uint32_t cnt2 = 0;
    {
        cnt1 += order2p_global[O_FULL] - order2p_global[O_POSR];
        cnt1 += (visual_residual_long_num + visual_residual_short_num) * 2 * VRNUM;
        cnt1 += (visual_residual_long_num + visual_residual_short_num) * 2;
        cnt1 += (visual_residual_long_num + visual_residual_short_num) * 2;
        cnt1 += (visual_residual_long_num + visual_residual_short_num) * 2;
        cnt1 += (visual_residual_long_num + visual_residual_short_num) * 2;
        cnt1 += (visual_residual_long_num + visual_residual_short_num) * 2;
        cnt1 += (visual_residual_long_num + visual_residual_short_num);
        cnt1 += visual_lhs_count * 6;
        for (int orderi = 0; orderi < O_IDEPTHR; ++orderi) {
            for (int orderj = orderi; orderj < O_IDEPTHR; ++orderj)
                cnt1 += order2size[orderi] * order2size[orderj];
        }
        cnt1 += order2p_local[O_IDEPTHR];
        cnt1 += order2p_local[O_IDEPTHR];
        cnt1 += mpos_size * r_size;
        cnt1 += mpos_size * mpos_size;
        cnt1 += mpos_size;
        cnt1 += order2p_local[O_IDEPTHR];
        cnt1 += idepth_map_factors_short.size();
        cnt1 += r_size;
        cnt1 += order2p_local[O_IDEPTHR];
        cnt1 += idepth_map_factors_short.size();
        cnt1 += idepth_map_factors_short.size();
        cnt1 += order2p_global[E_O_MARG_POS];
        memory = new double [cnt1];
    }
    {
        old_estimations_pointer = memory + cnt2; cnt2 += order2p_global[O_FULL] - order2p_global[O_POSR];
        visual_jacobians = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num) * 2 * VRNUM;
        visual_residuals = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num) * 2;
        visual_residuals_save = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num) * 2;
        visual_model_residual_accum = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num) * 2;
        visual_model_residual_accum_save = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num) * 2;
        visual_residuals_linerized = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num) * 2;
        visual_costs = memory + cnt2; cnt2 += (visual_residual_long_num + visual_residual_short_num);
        if (last_marg_info)init_distance_constraint_jacobian = Eigen::VectorXd (pos_para_count * 3);



        lhs_pos_idepth_global = memory + cnt2;

        for (auto it = idepth_map_factors_short.begin(); it != idepth_map_factors_short.end(); ++it) {
            it->second->lhs_pos_idepth = memory + cnt2;
            int max_pindex = *std::max_element(it->second->pindexs_set.begin(), it->second->pindexs_set.end());
            cnt2 += (max_pindex - it->second->min_pindex + 1) * 6;
        }
        for (auto& it : long_feature_factors) {
            it->lhs_pos_idepth = memory + cnt2;
            int max_pindex = *std::max_element(it->pindexs_set.begin(), it->pindexs_set.end());
            cnt2 += (max_pindex - it->min_pindex + 1) * 6;
        }
        lhs_posbias.resize(O_IDEPTHR * O_IDEPTHR, 0);
        for (int orderi = 0; orderi < O_IDEPTHR; ++orderi) {
            for (int orderj = orderi; orderj < O_IDEPTHR; ++orderj) {
                lhs_posbias[orderi * O_IDEPTHR + orderj] = memory + cnt2;
                cnt2 += order2size[orderi] * order2size[orderj];
            }
        }
        rhs_posbias_pointer = memory + cnt2; cnt2 += order2p_local[O_IDEPTHR];
        rhs_posbias_save_pointer = memory + cnt2; cnt2 += order2p_local[O_IDEPTHR];
        lhs_margpos_outside_pointer = memory + cnt2; cnt2 += mpos_size * r_size;
        lhs_margpos_pointer = memory + cnt2; cnt2 += mpos_size * mpos_size;
        rhs_margpos_pointer = memory + cnt2; cnt2 += mpos_size;
        inc_posbias_pointer = memory + cnt2; cnt2 += order2p_local[O_IDEPTHR];
        inc_idepth_short_pointer = memory + cnt2; cnt2 += idepth_map_factors_short.size();
        inc_outside_pointer = memory + cnt2; cnt2 += r_size;
        gauss_posbias_pointer = memory + cnt2; cnt2 += order2p_local[O_IDEPTHR];
        gauss_idepth_short_pointer = memory + cnt2; cnt2 += idepth_map_factors_short.size();
        idepth_short_save_pointer = memory + cnt2; cnt2 += idepth_map_factors_short.size();
        para_posebias_inside_save_pointer = memory + cnt2; cnt2 += order2p_global[E_O_MARG_POS];
    }
    ASSERT(cnt2 == cnt1);




    {
        if (ESTIMATE_SCALE && last_marg_info && fix_scale && !have_hist)
            init_distance_constraint_model_residual_accum = 0;
        if (last_marg_info && !have_hist)
            yaw_constraint_model_residual_accum.setZero();
        memset(visual_model_residual_accum, 0, sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);
    }


    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
        IMUGNSSPreFactor* imu_pre_factori = it->second;
        auto& evaluate_idx = imu_pre_factori->evaluate_idx;
        auto& orders = imu_pre_factori->orders;
        auto& order_idx = imu_pre_factori->order_idx;
        auto& parameter_size_sum = imu_pre_factori->parameter_size_sum;
        auto& parameter_sizes = imu_pre_factori->parameter_sizes;
        auto& phase_biases = imu_pre_factori->phase_biases;
        const auto& pindexi = imu_pre_factori->pindexi;
        const auto& pindexj = imu_pre_factori->pindexj;


        evaluate_idx.clear();
        orders.clear();
        order_idx.clear();
        parameter_sizes.clear();
        parameter_size_sum = 0;


        if (ESTIMATE_SCALE) {
            parameter_sizes.push_back(1); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 1;
            orders.push_back(O_SCALE);
            order_idx.push_back(0);
        }
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC) {
            parameter_sizes.push_back(6); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 6;
            orders.push_back(O_POSGLOBAL);
            order_idx.push_back(0);
        }

        if (ESTIMATE_EXTRINSIC) {
            parameter_sizes.push_back(6); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 6;
            orders.push_back(O_EXTRINSIC);
            order_idx.push_back(0);
        }

        parameter_sizes.push_back(6); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 6; orders.push_back(pindex2order[pindexi]); order_idx.push_back(0);
        parameter_sizes.push_back(9); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 9; orders.push_back(bindex2order[pindexi]); order_idx.push_back(0);
        parameter_sizes.push_back(6); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 6; orders.push_back(pindex2order[pindexj]); order_idx.push_back(0);
        parameter_sizes.push_back(9); evaluate_idx.push_back(parameter_size_sum); parameter_size_sum += 9; orders.push_back(bindex2order[pindexj]); order_idx.push_back(0);
        if (ESTIMATE_TD) {
            parameter_sizes.push_back(1); evaluate_idx.push_back(parameter_size_sum);
            parameter_size_sum += 1;
            orders.push_back(O_TD);
            order_idx.push_back(0);
        }
        if (ESTIMATE_ACC_SCALE) {
            parameter_sizes.push_back(3); evaluate_idx.push_back(parameter_size_sum);
            parameter_size_sum += 3;
            orders.push_back(O_ACC_S);
            order_idx.push_back(0);
        }
        if (ESTIMATE_PBG) {
            parameter_sizes.push_back(3); evaluate_idx.push_back(parameter_size_sum);
            parameter_size_sum += 3;
            orders.push_back(O_PBG);
            order_idx.push_back(0);
        }
        for (int i = 0; i < (int)phase_biases.size(); i++) {
            parameter_sizes.push_back(1); evaluate_idx.push_back(parameter_size_sum);
            if (short_phase_bias_idx_map.find(phase_biases[i]) != short_phase_bias_idx_map.end()) {
                orders.push_back(O_PHASE_BIAS_SHORT);
                order_idx.push_back(short_phase_bias_idx_map[phase_biases[i]]);
            } else {
                ASSERT(long_phase_bias_idx_map.find(phase_biases[i]) != long_phase_bias_idx_map.end());
                orders.push_back(O_PHASE_BIAS_LONG);
                order_idx.push_back(long_phase_bias_idx_map[phase_biases[i]]);
            }
            parameter_size_sum += 1;
        }



    }


    // imu_gradient = Eigen::VectorXd(imu_parameter_size);
    // imu_gauss = Eigen::VectorXd(imu_parameter_size);
    // imu_lhs = flower_focus::MatrixXd(imu_parameter_size, imu_parameter_size);
    // imu_rhs = flower_focus::VectorXd(imu_parameter_size);
}


void VisualInertialBase::SaveOutsidePointerRaw() {
    outside_pointer_raw.clear();
    for (int order = O_POSR; order < O_FULL; ++order) {
        if (order2pointers[order])outside_pointer_raw.push_back(order2pointers[order]);
    }
    for (int i = 0; i < (int)long_phase_bias.size(); i++)
        outside_pointer_raw.push_back(long_phase_bias[i]);


    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
        ASSERT(order2pointers[O_POSGLOBAL] == global_pos_pointer);
    if (ESTIMATE_SCALE)
        ASSERT(order2pointers[O_SCALE] == scale_pointer);
    if (ESTIMATE_TD)
        ASSERT(order2pointers[O_TD] == td_pointer);
    if (ESTIMATE_EXTRINSIC)
        ASSERT(order2pointers[O_EXTRINSIC] == extrinsic_pointer);
    if (ESTIMATE_ACC_SCALE)
        ASSERT(order2pointers[O_ACC_S] == acc_scale_pointer);
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    ASSERT(order2pointers[O_EXTRINSIC_RIGHT] == extrinsic2_pointer && extrinsic2_pointer);
#endif
    if (ESTIMATE_PBG)
        ASSERT(order2pointers[O_PBG] == pbg_pointer);
}


void VisualInertialBase::BuildPhaseBias() {
    std::set<double*>prev_phase_bias_set;
    std::set<double*>next_phase_bias_set;

    long_phase_bias.clear();
    short_phase_bias.clear();
    long_phase_bias_idx_map.clear();
    short_phase_bias_idx_map.clear();
    phase_bias_alone_num = 0;



    if (!prev_vibase) {
        for (auto it = phase_bias_pointers_set.begin(); it != phase_bias_pointers_set.end(); it++) {
            if (last_marg_info->keep_block_addr_set.find(*it) != last_marg_info->keep_block_addr_set.end())
                prev_phase_bias_set.insert(*it);
        }
    } else {
        for (auto it = phase_bias_pointers_set.begin(); it != phase_bias_pointers_set.end(); it++) {
            if (prev_vibase->phase_bias_pointers_set.find(*it) != prev_vibase->phase_bias_pointers_set.end())
                prev_phase_bias_set.insert(*it);
        }
    }

    if (next_vibase) {
        for (auto it = phase_bias_pointers_set.begin(); it != phase_bias_pointers_set.end(); it++) {
            if (next_vibase->phase_bias_pointers_set.find(*it) != next_vibase->phase_bias_pointers_set.end())
                next_phase_bias_set.insert(*it);
        }
    }
    if (!next_vibase && USE_RTK) {
        phase_bias_alone_num = 0;
        for (auto it = phase_bias_pointers_set.begin(); it != phase_bias_pointers_set.end(); it++)
            long_phase_bias.push_back(*it);//for integer ambiguity resolution

    } else {
        //phase biases that exists in previous block or next block are considered as long  phase biases.
        for (auto it = prev_phase_bias_set.begin(); it != prev_phase_bias_set.end(); it++) {
            if (next_phase_bias_set.find(*it) == next_phase_bias_set.end())
                long_phase_bias.push_back(*it);//phase biases that exists in previous block and not in the next block are considered as long short phase biases.
        }
        phase_bias_alone_num = long_phase_bias.size();
        for (auto it = next_phase_bias_set.begin(); it != next_phase_bias_set.end(); it++)
            long_phase_bias.push_back(*it);//phase biases that exists in next block are considered as long long  phase biases.


        for (auto it = phase_bias_pointers_set.begin(); it != phase_bias_pointers_set.end(); it++) {
            if (prev_phase_bias_set.find(*it) == prev_phase_bias_set.end() && next_phase_bias_set.find(*it) == next_phase_bias_set.end())
                short_phase_bias.push_back(*it);
        }


    }

    for (int i = 0; i < (int)long_phase_bias.size(); i++)
        long_phase_bias_idx_map[long_phase_bias[i]] = i;
    for (int i = 0; i < (int)short_phase_bias.size(); i++)
        short_phase_bias_idx_map[short_phase_bias[i]] = i;

    ASSERT(long_phase_bias.size() + short_phase_bias.size() == phase_bias_pointers_set.size());

}


void VisualInertialBase::FactorUpdates() {
    if (init) return;
    Reset();
    history_flag = false;
    nonlinear_quantity = NLQ_THRESHOLD + 100;
    long_feature_factors.clear();
    for (auto it = idepth_map_factors_long.begin(); it != idepth_map_factors_long.end(); ++it) {
        if (deliver_idepth_pointer.find(it->second->p_idepth) != deliver_idepth_pointer.end())long_feature_factors.push_back(it->second);
    }
    for (auto it = idepth_map_factors_long.begin(); it != idepth_map_factors_long.end(); ++it) {
        if (deliver_idepth_pointer.find(it->second->p_idepth) == deliver_idepth_pointer.end())long_feature_factors.push_back(it->second);
    }
    BuildPhaseBias();



    bias_para_count = bias_index2pointers.size();
    pos_para_count = pos_index2pointers.size();
    para_count = bias_para_count + pos_para_count;
    marg_count = bias_para_count - 2;
    idepth_long_count = long_feature_factors.size();

    outside_threshold = LEAK_NUM;
    outside_threshold = std::min(outside_threshold, marg_count);

    Rwcs.resize(pos_para_count);
    Pwcs.resize(pos_para_count);

    int j_idx = 0;
    O_MARG_BIAS = j_idx; j_idx += marg_count;
    E_O_MARG_BIAS = j_idx;
    O_PHASE_BIAS_SHORT = j_idx;
    if (short_phase_bias.size())j_idx++;
    O_MARG_POS = j_idx; j_idx += marg_count - outside_threshold;
    E_O_MARG_POS = j_idx;
    O_POSR = j_idx; j_idx += outside_threshold;
    O_BIAS0 = j_idx; j_idx++;
    O_POS0 = j_idx; j_idx++;
    O_PHASE_BIAS_LONG = j_idx;
    if (long_phase_bias.size())j_idx++;
    O_POSK = j_idx; j_idx++;
    O_BIASK = j_idx; j_idx++;
    O_POSR_NEXT = j_idx; j_idx += pos_para_count - bias_para_count;

    O_SCALE = j_idx;
    if (ESTIMATE_SCALE)j_idx++;
    O_POSGLOBAL = j_idx;
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC) j_idx++;
    O_TD = j_idx;
    if (ESTIMATE_TD) j_idx++;
    O_EXTRINSIC = j_idx;
    if (ESTIMATE_EXTRINSIC) j_idx++;
    O_ACC_S = j_idx;
    if (ESTIMATE_ACC_SCALE)j_idx++;
    O_EXTRINSIC_RIGHT = j_idx;
    if (ESTIMATE_EXTRINSIC2 && USE_STEREO)j_idx++;
    O_PBG = j_idx;
    if (ESTIMATE_PBG) j_idx++;

    O_IDEPTHR = j_idx; j_idx += idepth_long_count;
    O_FULL = j_idx;







    ConstructOrderingMap();
    if (last_marg_info)ConstructPriorIdx();
    ConstructBiasConnections();
    PrepareMatrix();
    SaveOutsidePointerRaw();
    double* parameter0;
    if (!ESTIMATE_R_WI_WC)
        parameter0 = order2pointers[O_POS0];
    else
        parameter0 = order2pointers[O_POSGLOBAL];

    Eigen::Vector3d Pwb(parameter0[0], parameter0[1], parameter0[2]);
    Eigen::Quaterniond Qwb(parameter0[6], parameter0[3], parameter0[4], parameter0[5]);
    Eigen::Matrix3d Rwb = Qwb.toRotationMatrix();
    InitMag = Rwb.transpose() * Eigen::Vector3d({0, 1, 0});
    if (last_marg_info) {
        if (!have_hist) {
            distance0 = 0;
            for (int pindex = 1; pindex < pos_para_count - 1; pindex++) {
                int orderi = pindex2order[pindex];
                int order0 = pindex2order[0];
                Eigen::Vector3d Pk(order2pointers[orderi]);
                Eigen::Vector3d P0(order2pointers[order0]);
                distance0 += (Pk - P0).norm();
            }
        }
    }

    num_alone = (int)long_feature_factors.size() - (int)deliver_idepth_pointer.size();

    int r_size = order2p_local[O_FULL] - order2p_local[O_POSR];
    gauss_outside = flower_focus::VectorXd(r_size);
    lhs_outside = flower_focus::MatrixXd(r_size, r_size);
    rhs_outside = flower_focus::VectorXd(r_size);
    gradient_outside = flower_focus::VectorXd(r_size);

    BuildCurPrevIdx();
    init = true;

// LOG_OUT<<"short_phase_bias.size():"<<short_phase_bias.size()<<std::endl;
}


void VisualInertialBase::ResetMem() {
    for (int orderi = 0; orderi < O_IDEPTHR; ++orderi) {
        for (int orderj = orderi; orderj < O_IDEPTHR; ++orderj)
            memset(lhs_posbias[orderi * O_IDEPTHR + orderj], 0, sizeof(double [order2size[orderi]*order2size[orderj]]));
    }
    memset(lhs_pos_idepth_global, 0, sizeof(double)*visual_lhs_count * 6);
    memset(rhs_posbias_pointer, 0, sizeof(double)*order2p_local[O_IDEPTHR]);
    memset(rhs_posbias_save_pointer, 0, sizeof(double)*order2p_local[O_IDEPTHR]);
    memset(inc_posbias_pointer, 0, sizeof(double)*order2p_local[O_IDEPTHR]);
}


void VisualInertialBase::setRPwc() {
    for (int i = 0; i < pos_para_count; ++i) {
        Rwcs[i] = Eigen::Quaterniond(pos_index2pointers[i][6],
                                     pos_index2pointers[i][3],
                                     pos_index2pointers[i][4],
                                     pos_index2pointers[i][5]).toRotationMatrix();
        Pwcs[i] = Eigen::Vector3d(pos_index2pointers[i][0],
                                  pos_index2pointers[i][1],
                                  pos_index2pointers[i][2]);
    }
#if USE_STEREO
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    tC2C1 = Eigen::Vector3d(extrinsic2_pointer[0], extrinsic2_pointer[1], extrinsic2_pointer[2] ) ;
    rC2C1 = Eigen::Quaterniond(extrinsic2_pointer[6], extrinsic2_pointer[3], extrinsic2_pointer[4], extrinsic2_pointer[5]).normalized().toRotationMatrix();
#else
    tC2C1 = TC2C1;
    rC2C1 = RC2C1;
#endif
#endif
}


void VisualInertialBase::VisualEvaluate(int pindexi, int pindexj,
                                        double* p_idepth,
                                        double* residuals, double** jacobians,
                                        const Eigen::Vector3d pts_i, const Eigen::Vector3d pts_j
#if USE_STEREO
                                        , bool is_stereo
#endif
                                       )  {
#define LEN 6
    double inv_dep_i = p_idepth[0];
    Eigen::Matrix3d Ri = Rwcs[pindexi];
    Eigen::Vector3d Pi = Pwcs[pindexi];

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_camera_j0 = Rwcs[pindexj].transpose() * ( Ri * pts_camera_i + Pi - Pwcs[pindexj]);

    Eigen::Vector3d pts_camera_j;
#if USE_STEREO
    if (is_stereo) pts_camera_j = rC2C1 * pts_camera_j0 * (scale_pointer ? scale_pointer[0] : 1) + tC2C1;
    else pts_camera_j = pts_camera_j0;
#else
    pts_camera_j = pts_camera_j0;
#endif
    double dep_j = pts_camera_j.z();
    double raw_dep_j = scale_pointer ? dep_j * scale_pointer[0] : dep_j;
    if (raw_dep_j < 0.2) {
        if (residuals)Eigen::Map<Eigen::Vector2d> (residuals).setZero();
        if (jacobians) {
            if (jacobians[0]) Eigen::Map<Eigen::Matrix<double, 2, LEN, Eigen::RowMajor>> (jacobians[0]).setZero();
            if (jacobians[1]) Eigen::Map<Eigen::Matrix<double, 2, LEN, Eigen::RowMajor>> (jacobians[1]).setZero();
            if (jacobians[2]) Eigen::Map<Eigen::Vector2d> (jacobians[2]).setZero();
#if USE_STEREO
            if (jacobians[3])Eigen::Map<Eigen::Vector2d> (jacobians[3]).setZero();
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
            if (jacobians[4]) Eigen::Map<Eigen::Matrix<double, 2, LEN, Eigen::RowMajor>> (jacobians[4]).setZero();
#endif
#endif
        }
        return;
    }

    if (residuals) {
        Eigen::Map<Eigen::Vector2d> residual(residuals);
        residual = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
        residual = SQRT_INFO_FEATURE * residual;
#if USE_STEREO
        if (is_stereo) residual *= STEREO_SCALE;
#endif
    }



    if (jacobians) {

        Eigen::Matrix<double, 2, 3> reduce(2, 3);
        reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),
               0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j);


        reduce = SQRT_INFO_FEATURE * reduce;
#if USE_STEREO
#if ESTIMATE_SCALE
        if (is_stereo)reduce = reduce * rC2C1 * STEREO_SCALE * scale_pointer[0];
#else
        if (is_stereo)reduce = reduce * rC2C1 * STEREO_SCALE;
#endif
#endif
        Eigen::Matrix<double, 2, 3> reduce1 = reduce * Rwcs[pindexj].transpose();
        Eigen::Matrix<double, 2, 3> reduce2 = reduce1 * Ri;

        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 2, LEN, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
            jacobian_pose_i.setZero();
            jacobian_pose_i.block(0, 0, 2, 3) = reduce1;
            jacobian_pose_i.block(0, 3, 2, 3) = reduce2 * -Utility::skewSymmetric(pts_camera_i);
        }
        if (jacobians[1]) {
            Eigen::Map<Eigen::Matrix<double, 2, LEN, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
            jacobian_pose_j.setZero();
            jacobian_pose_j.block(0, 0, 2, 3) = -reduce1;
            jacobian_pose_j.block(0, 3, 2, 3) = reduce * Utility::skewSymmetric(pts_camera_j0);
        }
        if (jacobians[2]) {
            Eigen::Map<Eigen::Vector2d> jacobian_feature(jacobians[2]);
            jacobian_feature = reduce2 *  pts_camera_i * -1.0 /  inv_dep_i;
        }
#if USE_STEREO
        if (jacobians[3]) {
            Eigen::Map<Eigen::Vector2d> jacobian_scale(jacobians[3]);
            jacobian_scale.setZero();
            if (is_stereo && fix_scale)jacobian_scale = reduce / (scale_pointer ? scale_pointer[0] : 1) * pts_camera_j0;
        }
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        if (jacobians[4]) {
            Eigen::Map<Eigen::Matrix<double, 2, LEN, Eigen::RowMajor>> jacobian_extrinsic(jacobians[4]);
            jacobian_extrinsic.setZero();
            if (is_stereo) {
                jacobian_extrinsic.block(0, 0, 2, 3) = reduce * rC2C1.transpose() / (scale_pointer ? scale_pointer[0] : 1);
                jacobian_extrinsic.block(0, 3, 2, 3) = -reduce * Utility::skewSymmetric(pts_camera_j0);
            }
        }
#endif
#endif
    }
#undef LEN
}


void VisualInertialBase::VisualJacobianResidualUpdatelhsRhs(int pindexi, int pindexj, Eigen::Vector3d pts_i, Eigen::Vector3d pts_j,
                                                            double* lhs_pos_idepth, double* rhs_idepth, double* lhs_idepth,
                                                            LossFunction* loss_function0,
                                                            double* jacobian_pointer,
                                                            double* visual_residual_raw_in,
                                                            double visual_cost,
                                                            int min_pindex,
                                                            double* p_idepth
#if USE_STEREO
                                                            , double* lhs_scale_idepth, bool is_stereo
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                                                            , double* lhs_extrinsic_idepth
#endif
#endif
                                                           ) {

#if !USE_STEREO
    double* visual_jacobians_raw_in[3];
#else
#if !ESTIMATE_EXTRINSIC2&&USE_STEREO
    double* visual_jacobians_raw_in[4];
#else
    double* visual_jacobians_raw_in[5];
#endif
#endif
    visual_jacobians_raw_in[0] = jacobian_pointer;
    visual_jacobians_raw_in[1] = jacobian_pointer + 6 * 2;
    visual_jacobians_raw_in[2] = jacobian_pointer + 6 * 4;
#if USE_STEREO
    visual_jacobians_raw_in[3] = jacobian_pointer + 6 * 4 + 2;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    visual_jacobians_raw_in[4] = jacobian_pointer + 6 * 4 + 2 + 2;
#endif
#endif
    int porderi = pindex2order[pindexi];
    int porderj = pindex2order[pindexj];

    VisualEvaluate(pindexi, pindexj, p_idepth, 0, visual_jacobians_raw_in, pts_i, pts_j
#if USE_STEREO
                   , is_stereo
#endif
                  );

    Eigen::Map<flower_focus::Matrix2_6d>visual_jacobians_in1(jacobian_pointer);
    Eigen::Map<flower_focus::Matrix2_6d>visual_jacobians_in2(jacobian_pointer + 6 * 2);
    Eigen::Map<flower_focus::Vector2d>visual_jacobians_in3(jacobian_pointer + 6 * 4);
#if USE_STEREO
    Eigen::Map<flower_focus::Vector2d>visual_jacobians_in4(visual_jacobians_raw_in[3]);
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    Eigen::Map<flower_focus::Matrix2_6d>visual_jacobians_in5(visual_jacobians_raw_in[4]);
#endif
#endif
    Eigen::Map<flower_focus::Vector2d>visual_residual_in(visual_residual_raw_in);

    {
        double  alpha_sq_norm_, rho[3], sqrt_rho1_;

        loss_function0->Evaluate(visual_cost, rho);
        sqrt_rho1_ = sqrt(rho[1]);
        if ((visual_cost == 0.0) || (rho[2] <= 0.0))
            alpha_sq_norm_ = 0.0;
        else {
            const double D = 1.0 + 2.0 * visual_cost * rho[2] / rho[1];
            const double alpha = 1.0 - sqrt(D);
            alpha_sq_norm_ = alpha / visual_cost;
        }
        flower_focus::Matrix2d visual_scale =
            sqrt_rho1_ * (flower_focus::Matrix2d::Identity() - alpha_sq_norm_ * visual_residual_in * visual_residual_in.transpose());

        visual_jacobians_in1 = visual_scale * visual_jacobians_in1;
        visual_jacobians_in2 = visual_scale * visual_jacobians_in2;
        visual_jacobians_in3 = visual_scale * visual_jacobians_in3;
#if USE_STEREO
        if (is_stereo) {
            visual_jacobians_in4 = visual_scale * visual_jacobians_in4;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
            visual_jacobians_in5 = visual_scale * visual_jacobians_in5;
#endif
        }
#endif
    }

    bool need_update_pos = pindexi != pindexj;
    double tmp[6];
    if (need_update_pos) {
        MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
        (visual_jacobians_raw_in[0], 2, 6, visual_jacobians_raw_in[0], 2, 6, lhs_posbias[porderi * O_IDEPTHR + porderi], 0, 0,  6, 6);
        if (porderi < porderj)
            MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
            (visual_jacobians_raw_in[0], 2, 6, visual_jacobians_raw_in[1], 2, 6, lhs_posbias[porderi * O_IDEPTHR + porderj], 0, 0,  6, 6);
        else
            MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
            (visual_jacobians_raw_in[1], 2, 6, visual_jacobians_raw_in[0], 2, 6, lhs_posbias[porderj * O_IDEPTHR + porderi], 0, 0,  6, 6);
        double tmp[6];
        MatrixTransposeVectorMultiply<2, 6, 0>
        (visual_jacobians_raw_in[0], 2, 6, visual_residual_raw_in, tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[porderi]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[porderi]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        MatrixTransposeMatrixMultiply<2, 6, 2, 1, 1>
        (visual_jacobians_raw_in[0], 2, 6, visual_jacobians_raw_in[2], 2, 1, lhs_pos_idepth + (pindexi - min_pindex) * 6, 0, 0, 6, 1);
    }

    if (need_update_pos) {
        MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
        (visual_jacobians_raw_in[1], 2, 6, visual_jacobians_raw_in[1], 2, 6, lhs_posbias[porderj * O_IDEPTHR + porderj], 0, 0, 6, 6);
        MatrixTransposeMatrixMultiply<2, 6, 2, 1, 1>
        (visual_jacobians_raw_in[1], 2, 6, visual_jacobians_raw_in[2], 2, 1, lhs_pos_idepth + (pindexj - min_pindex) * 6, 0, 0, 6, 1);
        MatrixTransposeVectorMultiply<2, 6, 0>
        (visual_jacobians_raw_in[1], 2, 6, visual_residual_raw_in, tmp );

        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[porderj]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[porderj]) += Eigen::Map<flower_focus::Vector6d>(tmp);
    }

    MatrixTransposeMatrixMultiply<2, 1, 2, 1, 1>
    (visual_jacobians_raw_in[2], 2, 1, visual_jacobians_raw_in[2], 2, 1, lhs_idepth, 0, 0, 1, 1);

    MatrixTransposeVectorMultiply<2, 1, 1>
    (visual_jacobians_raw_in[2], 2, 1, visual_residual_raw_in, rhs_idepth);



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
            MatrixTransposeMatrixMultiply<2, 6, 2, 1, 1>
            (visual_jacobians_raw_in[0], 2, 6, visual_jacobians_raw_in[3], 2, 1, lhs_posbias[porderi * O_IDEPTHR + O_SCALE], 0, 0,  6, 1);

            MatrixTransposeMatrixMultiply<2, 6, 2, 1, 1>
            (visual_jacobians_raw_in[1], 2, 6, visual_jacobians_raw_in[3], 2, 1, lhs_posbias[porderj * O_IDEPTHR + O_SCALE], 0, 0,  6, 1);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
            MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
            (visual_jacobians_raw_in[0], 2, 6, visual_jacobians_raw_in[4], 2, 6, lhs_posbias[porderi * O_IDEPTHR + O_EXTRINSIC_RIGHT], 0, 0,  6, 6);
            MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
            (visual_jacobians_raw_in[1], 2, 6, visual_jacobians_raw_in[4], 2, 6, lhs_posbias[porderj * O_IDEPTHR + O_EXTRINSIC_RIGHT], 0, 0,  6, 6);
#endif
        }

#if ESTIMATE_SCALE
        MatrixTransposeMatrixMultiply<2, 1, 2, 1, 1>
        (visual_jacobians_raw_in[3], 2, 1, visual_jacobians_raw_in[3], 2, 1, lhs_posbias[O_SCALE * O_IDEPTHR + O_SCALE], 0, 0,  1, 1);

        MatrixTransposeMatrixMultiply<2, 1, 2, 1, 1>
        (visual_jacobians_raw_in[2], 2, 1, visual_jacobians_raw_in[3], 2, 1, lhs_scale_idepth, 0, 0, 1, 1);

        MatrixTransposeVectorMultiply<2, 1, 0>
        (visual_jacobians_raw_in[3], 2, 1, visual_residual_raw_in, tmp );
        rhs_posbias_pointer[order2p_local[O_SCALE]] += tmp[0];
        rhs_posbias_save_pointer[order2p_local[O_SCALE]] += tmp[0];
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        MatrixTransposeMatrixMultiply<2, 1, 2, 6, 1>
        (visual_jacobians_raw_in[2], 2, 1, visual_jacobians_raw_in[4], 2, 6, lhs_extrinsic_idepth, 0, 0, 1, 6);
#if ESTIMATE_SCALE
        MatrixTransposeMatrixMultiply<2, 1, 2, 6, 1>
        (visual_jacobians_raw_in[3], 2, 1, visual_jacobians_raw_in[4], 2, 6, lhs_posbias[O_SCALE * O_IDEPTHR + O_EXTRINSIC_RIGHT], 0, 0,  1, 6);
#endif
        MatrixTransposeMatrixMultiply<2, 6, 2, 6, 1>
        (visual_jacobians_raw_in[4], 2, 6, visual_jacobians_raw_in[4], 2, 6, lhs_posbias[O_EXTRINSIC_RIGHT * O_IDEPTHR + O_EXTRINSIC_RIGHT], 0, 0,  6, 6);


        MatrixTransposeVectorMultiply<2, 6, 0>
        (visual_jacobians_raw_in[4], 2, 6, visual_residual_raw_in, tmp );
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT]) += Eigen::Map<flower_focus::Vector6d>(tmp);
        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[O_EXTRINSIC_RIGHT]) += Eigen::Map<flower_focus::Vector6d>(tmp);

#endif

    }
#endif

}



void VisualInertialBase::IMUJacobianResidualUpdatelhsRhs() {


    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
        IMUGNSSPreFactor* imu_pre_factori = it->second;

        const auto& evaluate_idx = imu_pre_factori->evaluate_idx;
        const auto& orders = imu_pre_factori->orders;
        const auto& order_idx = imu_pre_factori->order_idx;
        const auto& parameter_size_sum = imu_pre_factori->parameter_size_sum;
        const auto& parameter_sizes = imu_pre_factori->parameter_sizes;

        flower_focus::VectorXd imu_gradient(parameter_size_sum);
        flower_focus::MatrixXd imu_lhs(parameter_size_sum, parameter_size_sum);
        flower_focus::VectorXd imu_rhs(parameter_size_sum);

        ASSERT(imu_pre_factori->pindexi + 1 == imu_pre_factori->pindexj);
        imu_pre_factori->factor->EvaluateLhsRhs(imu_pre_factori->parameters.data(), imu_rhs.data(), imu_lhs.data(), imu_gradient.data(), mu);

        for (int i = 0; i < (int)orders.size(); ++i) {
            int orderi = orders[i];
            int order_idxi = order_idx[i];
            int sizei = parameter_sizes[i];
            int evaluate_idxi = evaluate_idx[i];
            int order_sizei = order2size[orderi];
            int lhs_posbias_indexi = orderi * O_IDEPTHR;

            Eigen::Map<flower_focus::VectorXd>(rhs_posbias_pointer + order2p_local[orderi] + order_idxi, sizei) += imu_rhs.segment(evaluate_idxi, sizei);
            Eigen::Map<flower_focus::VectorXd>(rhs_posbias_save_pointer + order2p_local[orderi] + order_idxi, sizei) += imu_gradient.segment(evaluate_idxi, sizei);
            // LOG_OUT<<"short_phase_bias.size():"<<short_phase_bias.size()<<std::endl;
            for (int j = 0; j <  (int)orders.size(); ++j) {
                int orderj = orders[j];
                int order_sizej = order2size[orderj];
                int order_idxj = order_idx[j];
                int evaluate_idxj = evaluate_idx[j];
                int sizej = parameter_sizes[j];
                if (orderi <= orderj)
                    Eigen::Map<flower_focus::MatrixXd> (lhs_posbias[lhs_posbias_indexi + orderj], order_sizei, order_sizej).block(order_idxi, order_idxj, sizei, sizej) += imu_lhs.block(evaluate_idxi, evaluate_idxj, sizei, sizej);
            }
        }
    }
}






void VisualInertialBase::EvaluateIdepthsShort() {
    int v_index = 0;
    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        VisualFactor* visual_factor = it1->second;
        visual_factor->rhs_idepth[0] = 0;
        visual_factor->lhs_idepth[0] = 0;
        auto& pindexs = visual_factor->pindexs;
        auto& ptss = visual_factor->ptss;
#if USE_STEREO
#if ESTIMATE_SCALE
        visual_factor->lhs_scale_idepth = 0;
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        memset(visual_factor->lhs_extrinsic_idepth, 0, sizeof(visual_factor->lhs_extrinsic_idepth));
#endif
#endif
        for (int i = 1; i < (int)pindexs.size(); ++i) {
            VisualJacobianResidualUpdatelhsRhs(pindexs[0], pindexs[i], ptss[0], ptss[i],
                                               visual_factor->lhs_pos_idepth,
                                               visual_factor->rhs_idepth,
                                               visual_factor->lhs_idepth, cauchy_loss_function,
                                               visual_jacobians + v_index * VRNUM * 2,
                                               visual_residuals + v_index * 2,
                                               visual_costs[v_index],
                                               visual_factor->min_pindex,
                                               visual_factor->p_idepth
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
}


void VisualInertialBase::EvaluateIdepthsLong() {
    int v_index = visual_residual_short_num;
    for (auto& it1 : long_feature_factors) {
        it1->rhs_idepth[0] = 0;
        it1->lhs_idepth[0] = 0;
        auto& pindexs = it1->pindexs;
        auto& ptss = it1->ptss;
#if USE_STEREO
#if ESTIMATE_SCALE
        it1->lhs_scale_idepth = 0;
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        memset(it1->lhs_extrinsic_idepth, 0, sizeof(it1->lhs_extrinsic_idepth));
#endif
#endif
        for (int i = 1; i < (int)pindexs.size(); ++i) {
            VisualJacobianResidualUpdatelhsRhs(pindexs[0], pindexs[i], ptss[0], ptss[i],
                                               it1->lhs_pos_idepth,
                                               it1->rhs_idepth,
                                               it1->lhs_idepth, cauchy_loss_function,
                                               visual_jacobians + v_index * VRNUM * 2,
                                               visual_residuals + v_index * 2,
                                               visual_costs[v_index],
                                               it1->min_pindex,
                                               it1->p_idepth
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


void VisualInertialBase::MargeIdepthsShort() {

    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        VisualFactor* visual_factor = it1->second;
        UpdateTrustRegion(visual_factor->lhs_idepth, 1, IDEPTH_REGION);//mu
        MargIdepth(visual_factor->pindexs_set,
                   visual_factor->lhs_pos_idepth,
                   visual_factor->rhs_idepth,
                   visual_factor->lhs_idepth,
                   visual_factor->min_pindex
#if USE_STEREO
                   , & (visual_factor->lhs_scale_idepth)
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                   , visual_factor->lhs_extrinsic_idepth
#endif
#endif
                  );
    }
}


void VisualInertialBase::MargIdepth(const std::set<int>& pindexs_set,
                                    double* lhs_pos_idepth, double* rhs_idepth, double* lhs_idepth,
                                    int min_pindex
#if USE_STEREO
                                    , double* lhs_scale_idepth
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                                    , double* lhs_extrinsic_idepth
#endif
#endif
                                   ) {

    lhs_idepth[0] = 1 / lhs_idepth[0];
    for (auto& pindexi : pindexs_set) {
        int porderi = pindex2order[pindexi];
        int lhs_posbias_indexi = porderi * O_IDEPTHR;
        flower_focus::Vector6d Anm_Amm_inverse = flower_focus::Vector6d(lhs_pos_idepth + 6 * (pindexi - min_pindex)) * lhs_idepth[0];

        Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[porderi]) -= Anm_Amm_inverse * rhs_idepth[0];

        for (auto& pindexj : pindexs_set) {
            int porderj = pindex2order[pindexj];
            if (porderi > porderj)
                continue;
            MatrixTransposeMatrixMultiply < 1, 6, 1, 6, -1 >
            (Anm_Amm_inverse.data(), 1, 6,
             lhs_pos_idepth + (pindexj - min_pindex) * 6, 1, 6,
             lhs_posbias[lhs_posbias_indexi + porderj], 0, 0,  6, 6);
        }
#if USE_STEREO
#if ESTIMATE_SCALE
        MatrixTransposeMatrixMultiply < 1, 6, 1, 1, -1 >
        (Anm_Amm_inverse.data(), 1, 6,
         lhs_scale_idepth, 1, 1,
         lhs_posbias[lhs_posbias_indexi + O_SCALE], 0, 0,  6, 1);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        MatrixTransposeMatrixMultiply < 1, 6, 1, 6, -1 >
        (Anm_Amm_inverse.data(), 1, 6,
         lhs_extrinsic_idepth, 1, 6,
         lhs_posbias[lhs_posbias_indexi + O_EXTRINSIC_RIGHT], 0, 0,  6, 6);
#endif

#endif
    }
#if USE_STEREO
#if ESTIMATE_SCALE
    rhs_posbias_pointer[order2p_local[O_SCALE]] -= lhs_scale_idepth[0] * lhs_idepth[0] * rhs_idepth[0];
    lhs_posbias[O_SCALE * O_IDEPTHR + O_SCALE][0] -= lhs_scale_idepth[0] * lhs_idepth[0] * lhs_scale_idepth[0];
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    Eigen::Map<flower_focus::Vector6d>(rhs_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT]) -= Eigen::Map<flower_focus::Vector6d>(lhs_extrinsic_idepth) * (lhs_idepth[0] * rhs_idepth[0]);
    Eigen::Map<flower_focus::Matrix6d>(lhs_posbias[O_EXTRINSIC_RIGHT * O_IDEPTHR + O_EXTRINSIC_RIGHT]) -= Eigen::Map<flower_focus::Vector6d>(lhs_extrinsic_idepth) * lhs_idepth[0] * Eigen::Map<flower_focus::Matrix1_6d>(lhs_extrinsic_idepth);
#if ESTIMATE_SCALE
    Eigen::Map<flower_focus::Matrix1_6d>(lhs_posbias[O_SCALE * O_IDEPTHR + O_EXTRINSIC_RIGHT]) -= lhs_scale_idepth[0] * lhs_idepth[0] * Eigen::Map<flower_focus::Matrix1_6d>(lhs_extrinsic_idepth);
#endif
#endif
#endif
}

void VisualInertialBase::MargeShortPhaseBias() {

    if (short_phase_bias.size()) {
        int short_phase_bias_num = order2size[O_PHASE_BIAS_SHORT];
        int marg_order = O_PHASE_BIAS_SHORT;
        ASSERT((int)short_phase_bias.size() == order2size[O_PHASE_BIAS_SHORT]);
        std::vector<int>connect_orders;
        for (int order = O_MARG_POS; order < O_POSR_NEXT; order++)connect_orders.push_back(order);
        for (int order = O_SCALE; order < O_IDEPTHR; order++)connect_orders.push_back(order);


        UpdateTrustRegion(lhs_posbias[marg_order * O_IDEPTHR + marg_order], order2size[marg_order], mu);
        int lhs_posbias_indexm = marg_order * O_IDEPTHR;
        double* rhs_posbiasm = rhs_posbias_pointer + order2p_local[marg_order];
        flower_focus::MatrixXd inverse_lhs =
            InvertPSDMatrix(Eigen::Map<flower_focus::MatrixXd>(lhs_posbias[lhs_posbias_indexm + marg_order], short_phase_bias_num, short_phase_bias_num) );

        memcpy(lhs_posbias[lhs_posbias_indexm + marg_order], inverse_lhs.data(), sizeof(double) * short_phase_bias_num * short_phase_bias_num);

        for (int j = 0; j < (int)connect_orders.size(); ++j) {
            int orderj = connect_orders[j];
            int sizej = order2size[orderj];
            flower_focus::MatrixXd Anm_Amm_inverse;
            int lhs_posbias_indexj = orderj * O_IDEPTHR;

            Anm_Amm_inverse = flower_focus::MatrixXd(sizej, short_phase_bias_num);
            MatrixTransposeMatrixMultiply<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, 0>(
                lhs_posbias[lhs_posbias_indexm + orderj], short_phase_bias_num, sizej,
                inverse_lhs.data(), short_phase_bias_num, short_phase_bias_num,
                Anm_Amm_inverse.data(), 0, 0, sizej, short_phase_bias_num);

            MatrixVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
                Anm_Amm_inverse.data(), sizej, short_phase_bias_num,
                rhs_posbiasm,
                rhs_posbias_pointer + order2p_local[orderj]);

            for (int k = j; k < (int)connect_orders.size(); k++) {
                int orderk = connect_orders[k];
                int sizek = order2size[orderk];
                ASSERT(orderk >= orderj);
                if (orderk != orderj) {
                    MatrixMatrixMultiply
                    < Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, -1 > (
                        Anm_Amm_inverse.data(), sizej, short_phase_bias_num,
                        lhs_posbias[lhs_posbias_indexm + orderk], short_phase_bias_num, sizek,
                        lhs_posbias[lhs_posbias_indexj + orderk], 0, 0, sizej, sizek, false);
                } else {
                    MatrixMatrixMultiply
                    < Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, -1 > (
                        Anm_Amm_inverse.data(), sizej, short_phase_bias_num,
                        lhs_posbias[lhs_posbias_indexm + orderk], short_phase_bias_num, sizek,
                        lhs_posbias[lhs_posbias_indexj + orderk], 0, 0, sizej, sizek, true);
                }
            }
            Eigen::Map<flower_focus::MatrixXd>(lhs_posbias[lhs_posbias_indexm + orderj], short_phase_bias_num, sizej) = Anm_Amm_inverse.transpose();
        }

    }
}

void VisualInertialBase::MargInsideBias() {


    for (int marg_order = 0; marg_order < E_O_MARG_BIAS; ++marg_order) {
        UpdateTrustRegion(lhs_posbias[marg_order * O_IDEPTHR + marg_order], order2size[marg_order], mu);
        const std::vector<int>& connect_orders = connections_bias[marg_order];
        int lhs_posbias_indexm = marg_order * O_IDEPTHR;
        double* rhs_posbiasm = rhs_posbias_pointer + order2p_local[marg_order];
        flower_focus::Matrix9d inverse_lhs =
            InvertPSDMatrix(Eigen::Map<flower_focus::Matrix9d>(lhs_posbias[lhs_posbias_indexm + marg_order]) );

        memcpy(lhs_posbias[lhs_posbias_indexm + marg_order], inverse_lhs.data(), sizeof(double) * 9 * 9);

        for (int j = 0; j < (int)connect_orders.size(); ++j) {
            int orderj = connect_orders[j];
            int sizej = order2size[orderj];
            flower_focus::MatrixX9d Anm_Amm_inverse;
            int lhs_posbias_indexj = orderj * O_IDEPTHR;

            Anm_Amm_inverse = flower_focus::MatrixX9d(sizej, 9);
            MatrixTransposeMatrixMultiply<9, Eigen::Dynamic, 9, 9, 0>(
                lhs_posbias[lhs_posbias_indexm + orderj], 9, sizej,
                inverse_lhs.data(), 9, 9,
                Anm_Amm_inverse.data(), 0, 0, sizej, 9);

            MatrixVectorMultiply < Eigen::Dynamic, 9, -1 > (
                Anm_Amm_inverse.data(), sizej, 9,
                rhs_posbiasm,
                rhs_posbias_pointer + order2p_local[orderj]);

            for (int k = j; k < (int)connect_orders.size(); k++) {
                int orderk = connect_orders[k];
                int sizek = order2size[orderk];
                ASSERT(orderk >= orderj);
                if (orderk != orderj) {
                    MatrixMatrixMultiply
                    < Eigen::Dynamic, 9, 9, Eigen::Dynamic, -1 > (
                        Anm_Amm_inverse.data(), sizej, 9,
                        lhs_posbias[lhs_posbias_indexm + orderk], 9, sizek,
                        lhs_posbias[lhs_posbias_indexj + orderk], 0, 0, sizej, sizek, false);
                } else {
                    MatrixMatrixMultiply
                    < Eigen::Dynamic, 9, 9, Eigen::Dynamic, -1 > (
                        Anm_Amm_inverse.data(), sizej, 9,
                        lhs_posbias[lhs_posbias_indexm + orderk], 9, sizek,
                        lhs_posbias[lhs_posbias_indexj + orderk], 0, 0, sizej, sizek, true);
                }
            }
            Eigen::Map<flower_focus::MatrixXd>(lhs_posbias[lhs_posbias_indexm + orderj], 9, sizej) = Anm_Amm_inverse.transpose();
        }
    }




}


void VisualInertialBase::MargInsidePos() {

    int r_size = order2p_local[O_FULL] - order2p_local[O_POSR];
    int mpos_size = order2p_local[O_POSR] - order2p_local[O_MARG_POS];
    Eigen::Map<flower_focus::MatrixXd>lhs_margpos_outside(lhs_margpos_outside_pointer, mpos_size, r_size);
    Eigen::Map<flower_focus::MatrixXd>lhs_margpos(lhs_margpos_pointer, mpos_size, mpos_size);

    lhs_outside.setZero();
    memset(lhs_margpos_outside_pointer, 0, sizeof(double)*mpos_size * r_size);


    int shift1 = order2p_local[O_MARG_POS];
    int shift2 = order2p_local[O_POSR];
    for (int marg_orderi = O_MARG_POS; marg_orderi < O_IDEPTHR; ++marg_orderi) {
        double* rhs_pos_biasi = rhs_posbias_pointer + order2p_local[marg_orderi];
        int sizei = order2size[marg_orderi];
        int marg_pindexi = order2p_local[marg_orderi] - shift1;
        int outside_pindexi = order2p_local[marg_orderi] - shift2;

        if (marg_orderi >= O_POSR)
            Eigen::Map<flower_focus::VectorXd>(rhs_outside.data() + outside_pindexi, sizei) = Eigen::Map<flower_focus::VectorXd>(rhs_pos_biasi, sizei);
        else
            Eigen::Map<flower_focus::VectorXd>(rhs_margpos_pointer + marg_pindexi, sizei) = Eigen::Map<flower_focus::VectorXd>(rhs_pos_biasi, sizei);
        int lhs_posbias_indexi = marg_orderi * O_IDEPTHR;
        for (int marg_orderj = marg_orderi; marg_orderj < O_IDEPTHR; ++marg_orderj) {
            double* lhs_pos_biasij = lhs_posbias[lhs_posbias_indexi + marg_orderj];
            int sizej = order2size[marg_orderj];
            int marg_indexj = order2p_local[marg_orderj] - shift1;
            int outside_indexj = order2p_local[marg_orderj] - shift2;

            if (marg_orderi >= O_POSR && marg_orderj >= O_POSR)
                lhs_outside.block(outside_pindexi, outside_indexj, sizei, sizej) = Eigen::Map<flower_focus::MatrixXd>(lhs_pos_biasij, sizei, sizej);
            else if (marg_orderi < O_POSR && marg_orderj < O_POSR)
                lhs_margpos.block(marg_pindexi, marg_indexj, sizei, sizej) = Eigen::Map<flower_focus::MatrixXd>(lhs_pos_biasij, sizei, sizej);
            else if (marg_orderi < O_POSR && marg_orderj >= O_POSR)
                lhs_margpos_outside.block(marg_pindexi, outside_indexj, sizei, sizej) = Eigen::Map<flower_focus::MatrixXd>(lhs_pos_biasij, sizei, sizej);
        }
    }

    int findex = order2p_local[O_IDEPTHR] - order2p_local[O_POSR];
    for (auto& it : long_feature_factors) {
        double* lhs_pos_idepth = it->lhs_pos_idepth;
        double* rhs_idepth = it->rhs_idepth;
        double* lhs_idepth = it->lhs_idepth;

        lhs_outside(findex, findex) = lhs_idepth[0];
        rhs_outside(findex) = rhs_idepth[0];
        std::set<int>& pindexs_set = it->pindexs_set;
        int min_pindex = it->min_pindex;
        for (auto& pindexi : pindexs_set) {
            int porderi = pindex2order[pindexi];
            int marg_pindexi = order2p_local[porderi] - shift1;
            int outside_pindexi = order2p_local[porderi] - shift2;
            if (outside_pindexi >= 0)
                lhs_outside.block(outside_pindexi, findex, 6, 1) = Eigen::Map<flower_focus::Vector6d>(lhs_pos_idepth + 6 * (pindexi - min_pindex));
            else
                lhs_margpos_outside.block(marg_pindexi, findex, 6, 1) = Eigen::Map<flower_focus::Vector6d>(lhs_pos_idepth + 6 * (pindexi - min_pindex));
        }
#if USE_STEREO
#if ESTIMATE_SCALE
        {
            int outside_pindexi = order2p_local[O_SCALE] - shift2;
            ASSERT(outside_pindexi > 0);
            lhs_outside(outside_pindexi, findex) = it->lhs_scale_idepth;
        }
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        {
            int outside_pindexi = order2p_local[O_EXTRINSIC_RIGHT] - shift2;
            ASSERT(outside_pindexi > 0);
            lhs_outside.block(outside_pindexi, findex, 6, 1) = Eigen::Map<flower_focus::Vector6d>(it->lhs_extrinsic_idepth);
        }
#endif
#endif
        findex++;
    }
    // LOG_OUT << lhs_outside << std::endl;
    // LOG_OUT << rhs_outside.transpose() << std::endl;
    int sizen = lhs_outside.rows();
    int sizem = order2p_local[O_POSR] - order2p_local[O_MARG_POS];
    if (sizem > 0) {
        lhs_margpos.diagonal().array() += mu;
        lhs_margpos = InvertPSDMatrix(lhs_margpos );
        flower_focus::MatrixXd Anm_Amm_inverse(sizen, sizem);
        MatrixTransposeMatrixMultiply<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, 0>(
            lhs_margpos_outside_pointer, sizem, sizen,
            lhs_margpos_pointer, sizem, sizem,
            Anm_Amm_inverse.data(), 0, 0, sizen, sizem);
        MatrixMatrixMultiply < Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, -1 > (
            Anm_Amm_inverse.data(),  sizen, sizem,
            lhs_margpos_outside_pointer, sizem, sizen,
            lhs_outside.data(), 0, 0, sizen, sizen, true);
        MatrixVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
            Anm_Amm_inverse.data(),  sizen, sizem,
            rhs_margpos_pointer,
            rhs_outside.data());
        Eigen::Map<flower_focus::MatrixXd>(lhs_margpos_outside_pointer, sizem, sizen) = Anm_Amm_inverse.transpose();
    }



    {
        for (int marg_orderi = O_POSR; marg_orderi < O_IDEPTHR; ++marg_orderi) {
            int sizei = order2size[marg_orderi];
            int outside_pindexi = order2p_local[marg_orderi] - shift2;
            gradient_outside.segment(outside_pindexi, sizei) = Eigen::Map<flower_focus::VectorXd>(rhs_posbias_save_pointer + order2p_local[marg_orderi], sizei);
        }

        int findex = order2p_local[O_IDEPTHR] - shift2;
        for (auto& it : long_feature_factors) {
            double* rhs_idepth_save = it->rhs_idepth_save;
            gradient_outside(findex) = rhs_idepth_save[0];
            findex++;
        }
    }


    if (last_marg_info)
        EvaluateLastMargInfo();

    lhs_outside = lhs_outside.selfadjointView<Eigen::Upper>();


}




void VisualInertialBase::UpdateMargposUseGradientAndGauss(double gradient_scale, double gauss_scale) {
    if (order2p_local[O_POSR] - order2p_local[O_MARG_POS] == 0)
        return;
    for (int order = O_MARG_POS; order < O_POSR; ++order) {
        Eigen::Map<flower_focus::Vector3d> P(order2pointers[order]);
        Eigen::Map<Eigen::Quaterniond> Q(order2pointers[order] + 3);
        double* pointer1 = gauss_posbias_pointer + order2p_local[order];
        double* pointer2 = rhs_posbias_save_pointer + order2p_local[order];
        double* pointer3 = inc_posbias_pointer + order2p_local[order];
        for (int i = 0; i < 6; ++i)
            *pointer3++ = gauss_scale * (*pointer1++) - gradient_scale * (*pointer2++);
        P -= Eigen::Map<flower_focus::Vector3d>(inc_posbias_pointer + order2p_local[order]);
        Q = (  Q *  Utility::deltaQ(-Eigen::Map<flower_focus::Vector3d>(inc_posbias_pointer + order2p_local[order] + 3))  ).normalized();
    }
}


void VisualInertialBase::UpdateMargBiasUseGradientAndGauss(double gradient_scale, double gauss_scale) {
    for (int marg_order = E_O_MARG_BIAS - 1; marg_order >= 0; marg_order--) {
        Eigen::Map<flower_focus::Vector9d> B(order2pointers[marg_order]);
        double* pointer1 = gauss_posbias_pointer + order2p_local[marg_order];
        double* pointer2 = rhs_posbias_save_pointer + order2p_local[marg_order];
        double* pointer3 = inc_posbias_pointer + order2p_local[marg_order];
        for (int i = 0; i < 9; ++i)
            *pointer3++ = gauss_scale * (*pointer1++) - gradient_scale * (*pointer2++);
        B -=  Eigen::Map<flower_focus::Vector9d>(inc_posbias_pointer + order2p_local[marg_order]);
    }
}

void VisualInertialBase::UpdateMargPhaseBiasUseGradientAndGauss(double gradient_scale, double gauss_scale) {
    {
        int marg_order = O_PHASE_BIAS_SHORT;
        Eigen::Map<flower_focus::Vector9d> B(order2pointers[marg_order]);
        double* pointer1 = gauss_posbias_pointer + order2p_local[marg_order];
        double* pointer2 = rhs_posbias_save_pointer + order2p_local[marg_order];
        double* pointer3 = inc_posbias_pointer + order2p_local[marg_order];
        for (int i = 0; i < (int)short_phase_bias.size(); ++i) {
            *pointer3 = gauss_scale * (*pointer1++) - gradient_scale * (*pointer2++);
            short_phase_bias[i][0] -= pointer3[0];
            pointer3++;
        }
    }
}


void VisualInertialBase::UpdateMargFeatureUseGradientAndGauss(double gradient_scale, double gauss_scale) {
    int index = 0;
    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        inc_idepth_short_pointer[index] =
            gauss_idepth_short_pointer[index] * gauss_scale
            - gradient_scale * it1->second->rhs_idepth_save[0];
        it1->second->p_idepth[0] -= inc_idepth_short_pointer[index];
        index++;
    }
}


void VisualInertialBase::UpdateInsideStateUseGradientAndGauss(double gradient_scale, double gauss_scale) {
    UpdateMargposUseGradientAndGauss(gradient_scale, gauss_scale);
    UpdateMargPhaseBiasUseGradientAndGauss(gradient_scale, gauss_scale);
    UpdateMargBiasUseGradientAndGauss(gradient_scale, gauss_scale);
    UpdateMargFeatureUseGradientAndGauss(gradient_scale, gauss_scale);
    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it)
        it->second->factor->UpdateInsideStateUseGaussAndGradient(it->second->parameters.data(), gradient_scale, gauss_scale);



}


void VisualInertialBase::SaveGaussStep() {
    for (int order = 0; order < O_IDEPTHR; ++order) {
        memcpy(gauss_posbias_pointer + order2p_local[order],
               inc_posbias_pointer + order2p_local[order],
               order2size[order]*sizeof(double));
    }
    if (idepth_map_factors_short.size())
        Eigen::Map<flower_focus::VectorXd>(gauss_idepth_short_pointer, idepth_map_factors_short.size()) =
            Eigen::Map<flower_focus::VectorXd>(inc_idepth_short_pointer, idepth_map_factors_short.size());
}


void VisualInertialBase::UpdateMargposGaussStep() {
    int mpos_size = order2p_local[O_POSR] - order2p_local[O_MARG_POS];
    int r_size = order2p_local[O_FULL] - order2p_local[O_POSR];
    if (mpos_size == 0)
        return;
    flower_focus::VectorXd margpos_inc(mpos_size);



    MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 0>(
        lhs_margpos_pointer, mpos_size, mpos_size,
        rhs_margpos_pointer,
        margpos_inc.data());

    MatrixVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
        lhs_margpos_outside_pointer, mpos_size, r_size,
        inc_outside_pointer,
        margpos_inc.data());

    for (int order = O_MARG_POS; order < O_POSR; ++order) {
        memcpy(inc_posbias_pointer + order2p_local[order],
               margpos_inc.data() + (order - O_MARG_POS) * 6, sizeof(double) * 6);

        double* pointer1 = inc_posbias_pointer + order2p_local[order];
        double* pointer3 = rhs_posbias_save_pointer + order2p_local[order];
        for (int k = 0; k < 6; k++) {
            gauss_newton_step_inside_square_norm += (*pointer1) * (*pointer1);
            ASSERT(!std::isnan(gauss_newton_step_inside_square_norm));
            gradient_dot_gauss_newton_inside += (*pointer1++) * (*pointer3++);
        }
    }
}


void VisualInertialBase::UpdateMargBiasGaussStep() {
    for (int marg_order = E_O_MARG_BIAS - 1; marg_order >= 0; marg_order--) {
        const std::vector<int>& connect_orders = connections_bias[marg_order];
        int lhs_posbias_indexm = marg_order * O_IDEPTHR;

        MatrixVectorMultiply<9, 9, 0>(
            lhs_posbias[lhs_posbias_indexm + marg_order], 9, 9,
            rhs_posbias_pointer + order2p_local[marg_order],
            inc_posbias_pointer + order2p_local[marg_order]);
        for (int j = (int)connect_orders.size() - 1; j >= 0; j--) {
            int indexj = connect_orders[j];
            int sizen = order2size[indexj];
            MatrixVectorMultiply < 9, Eigen::Dynamic, -1 > (
                lhs_posbias[lhs_posbias_indexm + indexj], 9, sizen,
                inc_posbias_pointer + order2p_local[indexj],
                inc_posbias_pointer + order2p_local[marg_order]);
        }
        double* pointer1 = inc_posbias_pointer + order2p_local[marg_order];
        double* pointer3 = rhs_posbias_save_pointer + order2p_local[marg_order];
        for (int k = 0; k < 9; k++) {
            gauss_newton_step_inside_square_norm += (*pointer1) * (*pointer1);
            ASSERT(!std::isnan(gauss_newton_step_inside_square_norm));
            gradient_dot_gauss_newton_inside += (*pointer1++) * (*pointer3++);
        }
    }
}

void VisualInertialBase::UpdateMargPhaseBiasGaussStep() {
    if (short_phase_bias.size()) {
        int short_phase_bias_num = order2size[O_PHASE_BIAS_SHORT];
        int marg_order = O_PHASE_BIAS_SHORT;
        ASSERT((int)short_phase_bias.size() == order2size[O_PHASE_BIAS_SHORT]);
        std::vector<int>connect_orders;
        for (int order = O_MARG_POS; order < O_POSR_NEXT; order++)connect_orders.push_back(order);
        for (int order = O_SCALE; order < O_IDEPTHR; order++)connect_orders.push_back(order);

        int lhs_posbias_indexm = marg_order * O_IDEPTHR;

        MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 0>(
            lhs_posbias[lhs_posbias_indexm + marg_order], short_phase_bias_num, short_phase_bias_num,
            rhs_posbias_pointer + order2p_local[marg_order],
            inc_posbias_pointer + order2p_local[marg_order]);
        for (int j = (int)connect_orders.size() - 1; j >= 0; j--) {
            int indexj = connect_orders[j];
            int sizen = order2size[indexj];
            MatrixVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
                lhs_posbias[lhs_posbias_indexm + indexj], short_phase_bias_num, sizen,
                inc_posbias_pointer + order2p_local[indexj],
                inc_posbias_pointer + order2p_local[marg_order]);
        }
        double* pointer1 = inc_posbias_pointer + order2p_local[marg_order];
        double* pointer3 = rhs_posbias_save_pointer + order2p_local[marg_order];
        for (int k = 0; k < short_phase_bias_num; k++) {
            gauss_newton_step_inside_square_norm += (*pointer1) * (*pointer1);
            ASSERT(!std::isnan(gauss_newton_step_inside_square_norm));
            gradient_dot_gauss_newton_inside += (*pointer1++) * (*pointer3++);
        }
    }
}


void VisualInertialBase::UpdateMargFeatureGaussStep() {
    int index = 0;
    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {

        VisualFactor* visual_factor = it1->second;
        double* lhs_pos_idepth = visual_factor->lhs_pos_idepth;
        double* rhs_idepth = visual_factor->rhs_idepth;
        double* lhs_idepth = visual_factor->lhs_idepth;
        const std::set<int>& pindexs_set = it1->second->pindexs_set;
        int min_pindex = it1->second->min_pindex;
        for (auto& pindexi : pindexs_set) {
            MatrixVectorMultiply < 1, 6, -1 > (
                lhs_pos_idepth + (pindexi - min_pindex) * 6, 1, 6,
                inc_posbias_pointer + order2p_local[pindex2order[pindexi]],
                rhs_idepth);
        }
#if USE_STEREO
#if ESTIMATE_SCALE
        rhs_idepth[0] -= visual_factor->lhs_scale_idepth * inc_posbias_pointer[order2p_local[O_SCALE]];
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        MatrixVectorMultiply < 1, 6, -1 > (
            visual_factor->lhs_extrinsic_idepth, 1, 6,
            inc_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT],
            rhs_idepth);
#endif
#endif
        double inc = lhs_idepth[0] * rhs_idepth[0];
        gauss_newton_step_inside_square_norm += inc * inc;
        ASSERT(!std::isnan(gauss_newton_step_inside_square_norm));
        gradient_dot_gauss_newton_inside += inc * visual_factor->rhs_idepth_save[0];
        inc_idepth_short_pointer[index] = inc;
        index++;

    }
}


void VisualInertialBase::UpdateInsideGaussStep() {

    UpdateMargposGaussStep();
    UpdateMargPhaseBiasGaussStep();
    UpdateMargBiasGaussStep();
    UpdateMargFeatureGaussStep();

}





void VisualInertialBase::SaveOrRestoreHidenStates(bool is_save) {
    int index = 0;
    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        double* src, *dst;
        if (is_save) {
            src = it1->second->p_idepth;
            dst = idepth_short_save_pointer + index;
        } else {
            src = idepth_short_save_pointer + index;
            dst = it1->second->p_idepth;
        }
        memcpy(dst, src, sizeof(double));
        index++;
    }

    for (int order = 0; order < E_O_MARG_POS; ++order) {
        if (!order2pointers[order])continue;
        double* src, *dst;
        if (is_save) {
            src = order2pointers[order];
            dst = para_posebias_inside_save_pointer + order2p_global[order];
        } else {
            src = para_posebias_inside_save_pointer + order2p_global[order];
            dst = order2pointers[order];
        }
        memcpy(dst, src, globalSize(order2size[order]) * sizeof(double));
    }

    for (int i = 0; i < (int)short_phase_bias.size(); i++) {
        double* src, *dst;
        if (is_save) {
            src = short_phase_bias[i];
            dst = para_posebias_inside_save_pointer + order2p_global[O_PHASE_BIAS_SHORT] + short_phase_bias_idx_map[short_phase_bias[i]];
        } else {
            src = para_posebias_inside_save_pointer + order2p_global[O_PHASE_BIAS_SHORT] + short_phase_bias_idx_map[short_phase_bias[i]];
            dst = short_phase_bias[i];
        }
        memcpy(dst, src, 1 * sizeof(double));
    }

    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it)
        it->second->factor->SaveOrRestoreHidenStates(is_save);
}


void VisualInertialBase::SaveGradientOutside() {
    int shift1 = order2p_local[O_POSR];

    double* p0 = gradient_outside.data() + order2p_local[O_IDEPTHR] - shift1;
    for (auto& it : long_feature_factors)
        it->rhs_idepth_save[0] = *p0++;

    for (int order = O_POSR; order < O_IDEPTHR; ++order)
        memcpy(rhs_posbias_save_pointer + order2p_local[order], gradient_outside.data() + order2p_local[order] - shift1, order2size[order] * sizeof(double));
}


void VisualInertialBase::SaveGaussStepOutside() {
    int shift1 = order2p_local[O_POSR];
    Eigen::Map<flower_focus::VectorXd>(inc_outside_pointer, order2p_local[O_FULL] - order2p_local[O_POSR]) =
        Eigen::Map<flower_focus::VectorXd>(gauss_outside.data(), order2p_local[O_FULL] - order2p_local[O_POSR]);

    for (int order = O_POSR; order < O_IDEPTHR; ++order)
        memcpy(inc_posbias_pointer + order2p_local[order], inc_outside_pointer + order2p_local[order] - shift1, order2size[order] * sizeof(double));
}





void VisualInertialBase::EvaluateAlpha(double* alpha1) {
    alpha1[0] = 0;



    // LOG_OUT<<Eigen::Map<flower_focus::VectorXd>(gauss_posbias_pointer, order2p_local[O_IDEPTHR]).transpose()<<std::endl;;
    for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
        IMUGNSSPreFactor* imu_pre_factor = it->second;

        const auto& evaluate_idx = imu_pre_factor->evaluate_idx;
        const auto& orders = imu_pre_factor->orders;
        const auto& order_idx = imu_pre_factor->order_idx;
        const auto& parameter_size_sum = imu_pre_factor->parameter_size_sum;
        const auto& parameter_sizes = imu_pre_factor->parameter_sizes;

        Eigen::VectorXd imu_gradient(parameter_size_sum);
        Eigen::VectorXd imu_gauss(parameter_size_sum);
        imu_gradient.setZero();
        imu_gauss.setZero();

        for (int i = 0; i < (int)orders.size(); i++) {
            imu_gradient.segment(evaluate_idx[i], parameter_sizes[i]) = Eigen::Map<flower_focus::VectorXd>(rhs_posbias_save_pointer + order2p_local[orders[i]] + order_idx[i], parameter_sizes[i]);
            imu_gauss.segment(evaluate_idx[i], parameter_sizes[i]) = Eigen::Map<flower_focus::VectorXd>(gauss_posbias_pointer + order2p_local[orders[i]] + order_idx[i], parameter_sizes[i]);
        }



        imu_pre_factor->factor->UpdateInsideGaussStep(imu_gradient.data(), imu_gauss.data(), &gauss_newton_step_inside_square_norm, alpha1, &gradient_squared_norm_inside, &gradient_dot_gauss_newton_inside);

    }



    int v_index = 0;
    for (auto it = idepth_map_factors_short.begin(); it != idepth_map_factors_short.end(); ++it) {
        auto& pindexs = it->second->pindexs;
#if USE_STEREO
        auto& is_stereo =  it->second->is_stereo;
#endif
        for (int i = 1; i < (int)pindexs.size(); ++i) {
#if USE_ASSERT&&USE_STEREO
            if (pindexs[i] == 0)ASSERT(is_stereo[i]);
#endif
            double* p0 = visual_jacobians + v_index * VRNUM * 2;
            flower_focus::Vector2d model_residuals =
                Eigen::Map<flower_focus::Matrix2_6d>(p0 + 0) * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[pindex2order[pindexs[0]]]) +
                Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 2) * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[pindex2order[pindexs[i]]]) +
                Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4) * Eigen::Map<flower_focus::Vector1d>(it->second->rhs_idepth_save);
#if USE_STEREO
            if (is_stereo[i]) {
#if ESTIMATE_SCALE
                model_residuals += Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4 + 2) * Eigen::Map<flower_focus::Vector1d>(rhs_posbias_save_pointer + order2p_local[O_SCALE]);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                model_residuals += Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 4 + 2 + 2) * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[O_EXTRINSIC_RIGHT]);
#endif
            }
#endif
            alpha1[0] += model_residuals.squaredNorm();
            v_index++;
        }
    }
    ASSERT(v_index == visual_residual_short_num);

    for (auto& it : long_feature_factors) {
        auto& pindexs = it->pindexs;
#if USE_STEREO
        auto& is_stereo =  it->is_stereo;
#endif
        for (int i = 1; i < (int)pindexs.size(); ++i) {
#if USE_ASSERT&&USE_STEREO
            if (pindexs[i] == 0)ASSERT(is_stereo[i]);
#endif
            double* p0 = visual_jacobians + v_index * VRNUM * 2;
            flower_focus::Vector2d model_residuals =
                Eigen::Map<flower_focus::Matrix2_6d>(p0 + 0) * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[pindex2order[pindexs[0]]]) +
                Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 2) * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[pindex2order[pindexs[i]]]) +
                Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4) * Eigen::Map<flower_focus::Vector1d>(it->rhs_idepth_save);
#if USE_STEREO
            if (is_stereo[i]) {
#if ESTIMATE_SCALE
                model_residuals += Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4 + 2) * Eigen::Map<flower_focus::Vector1d>(rhs_posbias_save_pointer + order2p_local[O_SCALE]);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                model_residuals += Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 4 + 2 + 2) * Eigen::Map<flower_focus::Vector6d>(rhs_posbias_save_pointer + order2p_local[O_EXTRINSIC_RIGHT]);
#endif
            }
#endif
            alpha1[0] += model_residuals.squaredNorm();
            v_index++;
        }
    }
    ASSERT(v_index == visual_residual_short_num + visual_residual_long_num);

    if (last_marg_info)
        EvaluatePriorAlpha(alpha1);
    int index = 0;
    for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
        VisualFactor* visual_factor = it1->second;
        gradient_squared_norm_inside += visual_factor->rhs_idepth_save[0] * visual_factor->rhs_idepth_save[0];
        index++;
    }

    for (int orderi = 0; orderi < E_O_MARG_POS; ++orderi) {
        int sizei = order2size[orderi];
        double* pointer1 = rhs_posbias_save_pointer + order2p_local[orderi];
        for (int k = 0; k < sizei; k++)
            gradient_squared_norm_inside += pointer1[k] * pointer1[k];
    }

}


void VisualInertialBase::EvaluateCost(double* cost, double* model_cost_change) {

    if (cost) cost[0] = 0;
    if (model_cost_change) model_cost_change[0] = 0;

    setRPwc();
    if (cost) {
        for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
            IMUGNSSPreFactor* imu_pre_factori = it->second;
            imu_pre_factori->factor->EvaluateCost(imu_pre_factori->parameters.data(), cost);
        }

        int v_index = 0;
        for (auto it1 = idepth_map_factors_short.begin(); it1 != idepth_map_factors_short.end(); it1++) {
            VisualFactor* visual_factor = it1->second;
            int pindexi = visual_factor->pindexs[0];
            auto& pindexs = visual_factor->pindexs;
            auto& ptss = visual_factor->ptss;
#if USE_STEREO
            auto& is_stereo = visual_factor->is_stereo;
#endif
            Eigen::Vector3d ptsInWorld = Rwcs[pindexi] * (ptss[0] / visual_factor->p_idepth[0]) + Pwcs[pindexi];
            for (int i = 1; i < (int)pindexs.size(); ++i) {
                int pindexj = pindexs[i];
                Eigen::Vector3d pts_j = ptss[i];
                Eigen::Vector3d pts_camera_j = Rwcs[pindexj].transpose() * (ptsInWorld - Pwcs[pindexj]);
#if USE_STEREO
#if ESTIMATE_SCALE
                if (is_stereo[i])pts_camera_j = rC2C1 * pts_camera_j * scale_pointer[0] + tC2C1;
#else
                if (is_stereo[i])pts_camera_j = rC2C1 * pts_camera_j  + tC2C1;
#endif
#endif
                Eigen::Map<flower_focus::Vector2d>visual_residual_in(visual_residuals + v_index * 2);
                visual_residual_in = (pts_camera_j / pts_camera_j.z()).head<2>() - pts_j.head<2>();
                visual_residual_in *= SQRT_INFO_FEATURE;
#if USE_STEREO
                if (is_stereo[i])visual_residual_in *= STEREO_SCALE;
#endif
                {
                    double sq_norm, rho[3], sqrt_rho1_;
                    sq_norm = visual_residual_in.squaredNorm();
                    visual_costs[v_index] = sq_norm;
                    cauchy_loss_function->Evaluate(sq_norm, rho);
                    sqrt_rho1_ = sqrt(rho[1]);
                    if ((sq_norm == 0.0) || (rho[2] <= 0.0))
                        visual_residual_in *= sqrt_rho1_;
                    else {
                        const double D = 1.0 + 2.0 * sq_norm * rho[2] / rho[1];
                        const double alpha = 1.0 - sqrt(D);
                        visual_residual_in *= sqrt_rho1_ / (1 - alpha);
                    }
                    cost[0] +=  0.5 * cauchy_loss_function->Evaluate0(sq_norm);
                }

                v_index++;
            }
        }
        ASSERT(v_index == visual_residual_short_num);
        for (auto& it1 : long_feature_factors) {
            int pindexi = it1->pindexs[0];
            auto& pindexs = it1->pindexs;
            auto& ptss = it1->ptss;
#if USE_STEREO
            auto& is_stereo = it1->is_stereo;
#endif
            Eigen::Vector3d ptsInWorld = Rwcs[pindexi] * (ptss[0] / it1->p_idepth[0]) + Pwcs[pindexi];
            for (int i = 1; i < (int)pindexs.size(); ++i) {
                int pindexj = pindexs[i];
                Eigen::Vector3d pts_j = ptss[i];
                Eigen::Vector3d pts_camera_j = Rwcs[pindexj].transpose() * (ptsInWorld - Pwcs[pindexj]);
#if USE_STEREO
#if ESTIMATE_SCALE
                if (is_stereo[i])pts_camera_j = rC2C1 * pts_camera_j * scale_pointer[0] + tC2C1;
#else
                if (is_stereo[i])pts_camera_j = rC2C1 * pts_camera_j  + tC2C1;
#endif
#endif
                Eigen::Map<flower_focus::Vector2d>visual_residual_in(visual_residuals + v_index * 2);
                visual_residual_in = (pts_camera_j / pts_camera_j.z()).head<2>() - pts_j.head<2>();
                visual_residual_in *= SQRT_INFO_FEATURE;
#if USE_STEREO
                if (is_stereo[i])visual_residual_in *= STEREO_SCALE;
#endif
                {
                    double sq_norm, rho[3], sqrt_rho1_;
                    sq_norm = visual_residual_in.squaredNorm();
                    visual_costs[v_index] = sq_norm;
                    cauchy_loss_function->Evaluate(sq_norm, rho);
                    sqrt_rho1_ = sqrt(rho[1]);
                    if ((sq_norm == 0.0) || (rho[2] <= 0.0))
                        visual_residual_in *= sqrt_rho1_;
                    else {
                        const double D = 1.0 + 2.0 * sq_norm * rho[2] / rho[1];
                        const double alpha = 1.0 - sqrt(D);
                        visual_residual_in *= sqrt_rho1_ / (1 - alpha);
                    }
                    cost[0] +=  0.5 * cauchy_loss_function->Evaluate0(sq_norm);
                }
                v_index++;
            }
        }
        ASSERT(v_index == visual_residual_short_num + visual_residual_long_num);
        if (last_marg_info)
            EvaluatePriorCost(cost);
        last_cost = cost[0];
    }




    if (model_cost_change) {

        for (auto it = imu_map_factors.begin(); it != imu_map_factors.end(); ++it) {
            IMUGNSSPreFactor* imu_pre_factori = it->second;
            imu_pre_factori->factor->EvaluateModelCostChange(model_cost_change);

        }



        int index = 0;
        int v_index = 0;
        for (auto it = idepth_map_factors_short.begin(); it != idepth_map_factors_short.end(); ++it) {
            auto& pindexs = it->second->pindexs;
#if USE_STEREO
            auto& is_stereo =  it->second->is_stereo;
#endif
            for (int i = 1; i < (int)pindexs.size(); ++i) {
#if USE_ASSERT&&USE_STEREO
                if (pindexs[i] == 0)ASSERT(is_stereo[i]);
#endif
                double* p0 = visual_jacobians + v_index * VRNUM * 2;
                flower_focus::Vector2d model_residuals =
                    Eigen::Map<flower_focus::Matrix2_6d>(p0 + 0) * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[pindex2order[pindexs[0]]]) +
                    Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 2) * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[pindex2order[pindexs[i]]]) +
                    Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4) * Eigen::Map<flower_focus::Vector1d>(inc_idepth_short_pointer + index);
#if USE_STEREO
                if (is_stereo[i]) {
#if ESTIMATE_SCALE
                    model_residuals += Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4 + 2) * Eigen::Map<flower_focus::Vector1d>(inc_posbias_pointer + order2p_local[O_SCALE]);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                    model_residuals += Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 4 + 2 + 2) * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT]);
#endif
                }
#endif
                model_cost_change[0] += model_residuals.transpose() * (Eigen::Map<flower_focus::Vector2d>(visual_residuals_save + v_index * 2) - model_residuals / 2);
                Eigen::Map<flower_focus::Vector2d>(visual_model_residual_accum + v_index * 2) += model_residuals;
                v_index++;
            }
            index++;
        }
        ASSERT(v_index == visual_residual_short_num);


        index = order2p_local[O_IDEPTHR] - order2p_local[O_POSR];
        for (auto& it : long_feature_factors) {
            auto& pindexs = it->pindexs;
#if USE_STEREO
            auto& is_stereo =  it->is_stereo;
#endif
            for (int i = 1; i < (int)pindexs.size(); ++i) {
#if USE_ASSERT&&USE_STEREO
                if (pindexs[i] == 0)ASSERT(is_stereo[i]);
#endif
                double* p0 = visual_jacobians + v_index * VRNUM * 2;
                flower_focus::Vector2d model_residuals =
                    Eigen::Map<flower_focus::Matrix2_6d>(p0 + 0) * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[pindex2order[pindexs[0]]]) +
                    Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 2) * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[pindex2order[pindexs[i]]]) +
                    Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4) * Eigen::Map<flower_focus::Vector1d>(inc_outside_pointer  + index);
#if USE_STEREO
                if (is_stereo[i]) {
#if ESTIMATE_SCALE
                    model_residuals += Eigen::Map<flower_focus::Vector2d>(p0 + 6 * 4 + 2) * Eigen::Map<flower_focus::Vector1d>(inc_posbias_pointer + order2p_local[O_SCALE]);
#endif
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                    model_residuals += Eigen::Map<flower_focus::Matrix2_6d>(p0 + 6 * 4 + 2 + 2) * Eigen::Map<flower_focus::Vector6d>(inc_posbias_pointer + order2p_local[O_EXTRINSIC_RIGHT]);
#endif
                }
#endif
                model_cost_change[0] += model_residuals.transpose() * (Eigen::Map<flower_focus::Vector2d>(visual_residuals_save + v_index * 2) - model_residuals / 2);
                Eigen::Map<flower_focus::Vector2d>(visual_model_residual_accum + v_index * 2) += model_residuals;
                v_index++;
            }
            index++;
        }
        ASSERT(v_index == visual_residual_short_num + visual_residual_long_num);

        if (last_marg_info)
            EvaluatePriorModelCostChange(model_cost_change);

    }

}

bool VisualInertialBase::EvaluateLhsRhs(double const* const* parameters,
                                        double* rhs, double* lhs,
                                        double* gradient,
                                        double mu_) {

    if (nonlinear_quantity_deliver_accum > NLQ_THRESHOLD) {
        if (parameters == 0) {
            rhs[0] = 0;
            return false;
        };
        mu = mu_;
        ASSERT(mu != 0);
        history_flag = true;

        ResetMem();
        setRPwc();


        EvaluateIdepthsShort();
        IMUJacobianResidualUpdatelhsRhs();
        EvaluateIdepthsLong();

        MargeIdepthsShort();
        MargInsideBias();
        MargeShortPhaseBias();
        MargInsidePos();




        {
            if (ESTIMATE_SCALE && last_marg_info && fix_scale && !have_hist) {
                init_distance_constraint_residual_linerized = init_distance_constraint_residual;
                init_distance_constraint_model_residual_accum = 0;
            }

            if (last_marg_info && !have_hist) {
                yaw_constraint_residual_linerized = yaw_constraint_residual;
                yaw_constraint_model_residual_accum.setZero();
            }

            memcpy(visual_residuals_linerized, visual_residuals, sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);
            memset(visual_model_residual_accum, 0, sizeof(double) * (visual_residual_long_num + visual_residual_short_num) * 2);



        }


    } else {
        ASSERT(NLQ_THRESHOLD != 0);
        memset(rhs_posbias_pointer, 0, sizeof(double)*order2p_local[O_IDEPTHR]);
        memset(rhs_posbias_save_pointer, 0, sizeof(double)*order2p_local[O_IDEPTHR]);
        memset(inc_posbias_pointer, 0, sizeof(double)*order2p_local[O_IDEPTHR]);
        setRPwc();
        ComputeRhs();
        MargeRhs();

    }


    return true;

}

