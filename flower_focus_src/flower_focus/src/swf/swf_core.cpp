
#include "swf.h"
#include "../solver/solver.h"
#include "../solver/solver_residualblock.h"
#include "../factor/initial_factor.h"
#include "../gnss/gnss_imu_factor.h"




inline void SWFOptimization::AddFeatures(FeaturePerId& it_per_id, std::vector<VisualInertialBase*>& visual_inertial_bases, std::set<int>& factor_index_set) {
    int imu_j = it_per_id.start_frame - 1;
    int imu_i = it_per_id.start_frame ;

    //是否为跨块的特征点
    bool is_cross = (it_per_id.start_frame / SWF_SIZE_IN != (it_per_id.endFrame() - 1) / SWF_SIZE_IN)
                    || (idepths_all[0].find(it_per_id.feature_id) != idepths_all[0].end() && last_marg_info->keep_block_addr_set.find(&idepths_all[0][it_per_id.feature_id]) != last_marg_info->keep_block_addr_set.end());
    //如果是跨块的特征点，将特征点的参考帧设置为边界帧
    if ( is_cross && it_per_id.start_frame % SWF_SIZE_IN != 0) {
        imu_i = (imu_i / SWF_SIZE_IN + 1) * SWF_SIZE_IN;
        if (factor_index_set.find(imu_i / SWF_SIZE_IN - 1) != factor_index_set.end()) {
            //标识特征点不能被factor_index因子边缘化
            visual_inertial_bases[imu_i / SWF_SIZE_IN - 1]->deliver_idepth_pointer[&idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id]] = 0;
        }
    }
    int factor_index = imu_i / SWF_SIZE_IN;

    ASSERT(idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] != 0);
    //如果是相邻跨块的特征点，则将特征点视为非跨块的
    if (imu_i != 0 && (is_cross && (it_per_id.endFrame() - imu_i) <= LEAK_NUM)) {
        is_cross = false;
        factor_index -= 1;
        if (factor_index_set.find(imu_i / SWF_SIZE_IN - 1) != factor_index_set.end())
            visual_inertial_bases[imu_i / SWF_SIZE_IN - 1]->deliver_idepth_pointer.erase(&idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id]);
    }

    for (auto& it_per_frame : it_per_id.feature_per_frame) {
        imu_j++;
        if (imu_i == imu_j)continue;
        //跨块特征点改变参考帧
        if (imu_j - 1 == (imu_i / SWF_SIZE_IN + 1)*SWF_SIZE_IN) {
            ASSERT(is_cross);

            double* idepth_pointeri = &idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id];
            double* idepth_pointerj = &idepths_all[(imu_j - 1) / SWF_SIZE_IN][it_per_id.feature_id];

            if (idepth_pointerj[0] == 0) {
                Eigen::Vector3d ptsInW =
                    Rs[imu_i] * (
                        RIC[0] * (it_per_id.feature_per_frame[imu_i - it_per_id.start_frame].point /
                                  idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id]) + TIC[0]
                    ) + Ps[imu_i];
                Vector3d pts_cj = RIC[0].transpose() * ( Rs[imu_j - 1].transpose() * (ptsInW - Ps[imu_j - 1]) - TIC[0]);
                idepth_pointerj[0] = 1.0 / pts_cj.z();
            }
            if (factor_index_set.find(imu_i / SWF_SIZE_IN) != factor_index_set.end())
                visual_inertial_bases[imu_i / SWF_SIZE_IN]->deliver_idepth_pointer[idepth_pointeri] = idepth_pointerj;
            imu_i = imu_j - 1;
        }


        double* idepth_pointer = &idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id];

        ASSERT(imu_i <= it_per_id.endFrame());
        ASSERT(imu_j <= it_per_id.endFrame());
        Vector3d pts_i = it_per_id.feature_per_frame[imu_i - it_per_id.start_frame].point;

        if (idepth_pointer[0] == 0) {
            Eigen::Vector3d ptsInW =
                Rs[it_per_id.start_frame] * (RIC[0] * (it_per_id.feature_per_frame[0].point / idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id]) + TIC[0]) + Ps[it_per_id.start_frame];
            Vector3d pts_cj = RIC[0].transpose() * ( Rs[imu_i].transpose() * (ptsInW - Ps[imu_i]) - TIC[0]);
            idepth_pointer[0] = 1.0 / pts_cj.z();
            ASSERT(idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] != 0);
        }
        ASSERT(idepth_pointer[0] != 0);

        //非跨块特征点
        if (!is_cross) {
            if (factor_index_set.find(factor_index) == factor_index_set.end())continue;
            std::vector<double*>parameters{para_pose[imu_i], para_pose[imu_j], idepth_pointer};

            int pindexi, pindexj;
            std::unordered_map<double*, int> pos_pointer2indexs = visual_inertial_bases[factor_index]->pos_pointer2indexs;
            if (pos_pointer2indexs.find(parameters[0]) != pos_pointer2indexs.end()) pindexi = pos_pointer2indexs[parameters[0]];
            else pindexi = pos_pointer2indexs.size();

            if (pos_pointer2indexs.find(parameters[1]) != pos_pointer2indexs.end()) pindexj = pos_pointer2indexs[parameters[1]];
            else pindexj = pos_pointer2indexs.size();

            ASSERT(pindexi == imu_i - factor_index * SWF_SIZE_IN);
            ASSERT(pindexj == imu_j - factor_index * SWF_SIZE_IN);

            if (imu_i != imu_j)
                visual_inertial_bases[factor_index]->AddVisualFactorShort( it_per_id.feature_id,
                                                                           parameters[0], parameters[1], parameters[2],
                                                                           pindexi, pindexj,
                                                                           pts_i, it_per_frame.point
#if USE_STEREO
                                                                           , false
#endif
                                                                         );
#if USE_STEREO
            if (it_per_frame.is_stereo )//&& 1 / parameters[2][0] < 20 && imu_i == imu_j
                visual_inertial_bases[factor_index]->AddVisualFactorShort( it_per_id.feature_id,
                                                                           parameters[0], parameters[1], parameters[2],
                                                                           pindexi, pindexj,
                                                                           pts_i, it_per_frame.pointRight
#if USE_STEREO
                                                                           , true
#endif
                                                                         );
#endif
        } else {
            //跨块特征点
            int factor_index = (imu_j - 1) / SWF_SIZE_IN;
            if (factor_index_set.find(factor_index) == factor_index_set.end())continue;

            std::vector<double*>parameters{para_pose[imu_i], para_pose[imu_j], idepth_pointer};

            std::unordered_map<double*, int> pos_pointer2indexs = visual_inertial_bases[factor_index]->pos_pointer2indexs;
            int pindexi;
            ASSERT(pos_pointer2indexs.find(parameters[1]) != pos_pointer2indexs.end());
            if (pos_pointer2indexs.find(parameters[0]) != pos_pointer2indexs.end()) pindexi = pos_pointer2indexs[parameters[0]];
            else pindexi = pos_pointer2indexs.size();
            int pindexj = imu_j - factor_index * SWF_SIZE_IN;
            if (imu_i != imu_j)
                visual_inertial_bases[factor_index]->AddVisualFactorLong( it_per_id.feature_id,
                                                                          parameters[0], parameters[1], parameters[2],
                                                                          pindexi, pindexj,
                                                                          pts_i, it_per_frame.point
#if USE_STEREO
                                                                          , false
#endif
                                                                        );
#if USE_STEREO
            if (it_per_frame.is_stereo)// && 1 / parameters[2][0] < 20 && imu_i == imu_j
                visual_inertial_bases[factor_index]->AddVisualFactorLong( it_per_id.feature_id,
                                                                          parameters[0], parameters[1], parameters[2],
                                                                          pindexi, pindexj,
                                                                          pts_i, it_per_frame.pointRight
#if USE_STEREO
                                                                          , true
#endif
                                                                        );
#endif
        }
    }
}


void SWFOptimization::AddFactors(std::vector<VisualInertialBase*>& visual_inertial_bases, std::set<int>factor_index_set) {

    //添加IMU因子
    for (int j = 1; j < image_count; j++) {
        int factor_index = (j - 1) / SWF_SIZE_IN;
        if (factor_index_set.find(factor_index) == factor_index_set.end())continue;

        IMUGNSSFactor* factor = new IMUGNSSFactor(gnss_imu_infos[j]->imu_gnss_base);
        IMUGNSSPreFactor* imu_pre_factor = new IMUGNSSPreFactor(factor);
        std::vector<double*>parameters;

        if (ESTIMATE_SCALE)parameters.push_back(&scale_factor);
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)parameters.push_back(para_global_pos);
        if (ESTIMATE_EXTRINSIC) parameters.push_back(para_extrinsic);

        parameters.push_back(para_pose[j - 1]);
        parameters.push_back(para_speed_bias[j - 1]);
        parameters.push_back(para_pose[j]);
        parameters.push_back(para_speed_bias[j]);

        if (ESTIMATE_TD) parameters.push_back(&td);
        if (ESTIMATE_ACC_SCALE) parameters.push_back(acc_scale.data());
        if (ESTIMATE_PBG) parameters.push_back(Pbg.data());

        for (int i = 0; i < (int)gnss_imu_infos[j]->imu_gnss_base->gnss_phase_biases.size(); i++)
            parameters.push_back(gnss_imu_infos[j]->imu_gnss_base->gnss_phase_biases[i]);


        imu_pre_factor->pindexi = (j - 1) - factor_index * SWF_SIZE_IN;
        imu_pre_factor->pindexj = j - factor_index * SWF_SIZE_IN;
        imu_pre_factor->parameters = parameters;
        imu_pre_factor->phase_biases = gnss_imu_infos[j]->imu_gnss_base->gnss_phase_biases;
        visual_inertial_bases[factor_index]->AddIMUFactor(imu_pre_factor);
    }
    //添加视觉因子
    for (auto& it_per_id : f_manager.feature) {
        if (!it_per_id.valid)continue;
        AddFeatures(it_per_id, visual_inertial_bases, factor_index_set);
    }

}
//添加全局因子，全局因子会被反复利用
void SWFOptimization::UpdataGlobalFactors() {

    int max_factor_num = (image_count - 2) / SWF_SIZE_IN - 2;
    int old_max_factor_num = -1;
    for (int i = 0; i < max_factor_num; i++) {
        if (visual_inertial_bases_global[i] == 0) {
            old_max_factor_num = i;
            break;
        }
    }

    if (old_max_factor_num < 0)return;
    ASSERT(old_max_factor_num != max_factor_num);
    for (int i = old_max_factor_num; i < max_factor_num; i++)
        visual_inertial_bases_global[i] = new VisualInertialBase();

    std::set<int>factor_index_set;
    for (int i = old_max_factor_num; i < max_factor_num; i++)
        factor_index_set.insert(i);
    AddFactors(visual_inertial_bases_global, factor_index_set);
    if (old_max_factor_num == 0)visual_inertial_bases_global[0]->AddLastMargeInfo(last_marg_info);
    LOG_OUT << "add new\r\n";
}

void SWFOptimization::OptimizationOrMarginalization( int mode) {

    UpdataGlobalFactors();
    std::vector<VisualInertialBase*> visual_inertial_bases2;
    std::vector<VisualInertialBase*> visual_inertial_bases;
    visual_inertial_bases2.resize(image_count / SWF_SIZE_IN + 1, 0);
    idepth_pointers.clear();
    for (int i = 0; i < image_count / SWF_SIZE_IN + 1; i++)
        visual_inertial_bases2[i] = new VisualInertialBase();

    int old_leak_num = LEAK_NUM;
    std::set<int>factor_index_set;
    if (mode == MargeMode) {
        LEAK_NUM = 0;
        factor_index_set.insert(0);
        factor_index_set.insert(1);
        visual_inertial_bases.push_back(visual_inertial_bases2[0]);
        visual_inertial_bases.push_back(visual_inertial_bases2[1]);
    } else {
        for (int i = 0; i < image_count / SWF_SIZE_IN + 1; i++) {
            if (visual_inertial_bases_global[i] == 0) {
                factor_index_set.insert(i);
                visual_inertial_bases.push_back(visual_inertial_bases2[i]);
            } else
                visual_inertial_bases.push_back(visual_inertial_bases_global[i]);
        }
    }
    if (factor_index_set.size())AddFactors(visual_inertial_bases2, factor_index_set);

    LEAK_NUM = old_leak_num;
    for (int factor_index = 0; factor_index < (int)visual_inertial_bases.size(); factor_index++) {
        if (visual_inertial_bases[factor_index]->imu_map_factors.size() == 0) continue;
        for (auto it = visual_inertial_bases[factor_index]->idepth_map_factors_short.begin(); it != visual_inertial_bases[factor_index]->idepth_map_factors_short.end(); it++)
            idepth_pointers.insert(it->second->p_idepth);
        for (auto it = visual_inertial_bases[factor_index]->idepth_map_factors_long.begin(); it != visual_inertial_bases[factor_index]->idepth_map_factors_long.end(); it++)
            idepth_pointers.insert(it->second->p_idepth);

        for (auto it = visual_inertial_bases[factor_index]->deliver_idepth_pointer.begin(); it != visual_inertial_bases[factor_index]->deliver_idepth_pointer.end(); it++) {
            if (it->second)
                idepth_pointers.insert(it->second);
        }
    }

    for (auto it = idepth_pointers.begin(); it != idepth_pointers.end(); it++)
        (*it)[0] *= scale_factor;


    SolverInfo solver_info;
    VI_info_pointer_head = visual_inertial_bases[0];
    if (visual_inertial_bases[0]->last_marg_info == 0)visual_inertial_bases[0]->AddLastMargeInfo(last_marg_info);

    for (int factor_index = 0; factor_index <  (int)visual_inertial_bases.size(); factor_index++) {
        if (visual_inertial_bases[factor_index]->imu_map_factors.size() == 0) continue;
        if (factor_index != 0 && mode == MargeMode)continue;

        VI_info_pointer_tail = visual_inertial_bases[factor_index];

        if (factor_index != 0)
            visual_inertial_bases[factor_index]->prev_vibase = visual_inertial_bases[factor_index - 1];
        else visual_inertial_bases[factor_index]->prev_vibase = 0;

        if (factor_index != image_count / SWF_SIZE_IN && visual_inertial_bases[factor_index + 1] && visual_inertial_bases[factor_index + 1]->imu_map_factors.size())
            visual_inertial_bases[factor_index]->next_vibase = visual_inertial_bases[factor_index + 1];
        else visual_inertial_bases[factor_index]->next_vibase = 0;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
        visual_inertial_bases[factor_index]->extrinsic2_pointer = para_extrinsic2;
#endif
        visual_inertial_bases[factor_index]->FactorUpdates();

        solver_info.addResidualBlockInfo(new SolverResidualBlockInfo(
                                             new VisualInertialFactor(visual_inertial_bases[factor_index])
                                             , 0, visual_inertial_bases[factor_index]->outside_pointer_raw));

    }
    if (mode == MargeMode) {
        visual_inertial_bases[0]->next_vibase = 0;
        ASSERT(visual_inertial_bases[0]->next_vibase == 0);
        ASSERT(visual_inertial_bases[0]->prev_vibase == 0);
    }

    for (int factor_index = 0; factor_index <  (int)visual_inertial_bases.size(); factor_index++) {
        if (visual_inertial_bases_global[factor_index]) visual_inertial_bases_global[factor_index]->FactorUpdates();
    }



    TicToc t_marg;

    solver_info.init_solver();
    if (mode == NormalMode) {
        solver_info.solve(8);
        static double t_sum = 0;
        static int t_count = 0;
        t_sum += t_marg.toc();
        t_count++;
        LOG_OUT << t_marg.toc() << "," << t_sum / t_count << std::endl;
        if (needNreaolve && VI_info_pointer_tail->long_phase_bias.size()) {
            std::vector<double*>new_parameters;
            std::vector<int>new_sizes;
            flower_focus::MatrixXd lhs;
            Eigen::VectorXd rhs;
            ASSERT(VI_info_pointer_tail->phase_bias_alone_num == 0);
            VI_info_pointer_tail->GetMatrixNext(lhs,  rhs, new_parameters, new_sizes);
            std::vector<mea_t*>rovers;
            for (int i = 0; i < (int)gnss_imu_infos.size(); i++) {
                if (gnss_imu_infos[i] == 0)continue;
                for (int j = 0; j < (int)gnss_imu_infos[i]->Rovers.size(); j++)
                    rovers.push_back(gnss_imu_infos[i]->Rovers[j]);
            }
            ASSERT(rovers.size());
            gnss_imu_infos[image_count - 1]->LambdaSearch(rovers, new_parameters, new_sizes, lhs, rhs);
        }


    } else {
        solver_info.marginalization_process();
        std::vector<double*>new_parameters;
        std::vector<int>new_sizes;
        flower_focus::MatrixXd lhs;
        Eigen::VectorXd rhs;

        visual_inertial_bases[0]->GetMatrixNext(lhs,  rhs, new_parameters, new_sizes);

        MarginalizationInfo* marginalization_info = new MarginalizationInfo();
        if (USE_REPARAM) {
            flower_focus::MatrixXd lhs2 = flower_focus::MatrixXd::Zero(lhs.cols() + 6 * (SWF_SIZE_IN - 1), lhs.cols() + 6 * (SWF_SIZE_IN - 1));
            Eigen::VectorXd rhs2 = Eigen::VectorXd::Zero(lhs.cols() + 6 * (SWF_SIZE_IN - 1));
            lhs2.block(0, 0, lhs.cols(), lhs.cols()) = lhs;
            rhs2.segment(0, lhs.cols()) = rhs;
            for (int i = SWF_SIZE_IN + 1; i < SWF_SIZE_IN * 2; i++) {
                new_parameters.push_back(para_pose[i]);
                new_sizes.push_back(7);
            }
            lhs = lhs2;
            rhs = rhs2;
        }



        marginalization_info->setmarginalizeinfo(new_parameters, new_sizes, lhs, rhs, true);
        marginalization_info->getParameterBlocks(true);

        if (last_marg_info)
            delete last_marg_info;
        last_marg_info = marginalization_info;
        if (USE_REPARAM) {
            have_hist = true;
            set_last_marg_info_idx(SWF_SIZE_IN);
            last_marg_info->UpdateF();
        }


        for (int i = 0; i < (int)visual_inertial_bases_global.size(); i++) {
            if (visual_inertial_bases_global[i]) {
                delete visual_inertial_bases_global[i];
                visual_inertial_bases_global[i] = 0;
            }
        }
    }

    for (int i = 0; i < (int)visual_inertial_bases_global.size() - 1; i++) {
        if (visual_inertial_bases_global[i] && !visual_inertial_bases_global[i + 1])visual_inertial_bases_global[i]->next_vibase = 0;
    }

    for (auto it = idepth_pointers.begin(); it != idepth_pointers.end(); it++)
        (*it)[0] /= scale_factor;

    for (int i = 0; i < image_count / SWF_SIZE_IN + 1; i++) {
        if (visual_inertial_bases2[i]) delete visual_inertial_bases2[i];
    }

    LOG_OUT << "scale_factor:" << scale_factor << std::endl;

}



void SWFOptimization::set_last_marg_info_idx(int idx0) {


    std::set<double*>idepth_pointers;
    std::set<double*>phase_bias_pointers;


    for (int i = 0; i < (int)idepths_all.size(); i++) {
        for (auto it = idepths_all[i].begin(); it != idepths_all[i].end(); it++)
            idepth_pointers.insert(&(it->second));
    }
    for (int i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++)
            phase_bias_pointers.insert(&(it1->value));
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++)
            phase_bias_pointers.insert(&(it1->value));
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++)
            phase_bias_pointers.insert(&(it1->value));
        phase_bias_pointers.insert(&blackvalue);
    }


    if (last_marg_info) {
        int count2 = 0;
        int count1 = 0;
        last_marg_info->p_is.resize(SWF_SIZE_IN, -1);
        last_marg_info->v_i = -1;
        last_marg_info->extrinsic_i = -1;
        last_marg_info->lambda_is.clear();
        last_marg_info->R_WI_WC_i = -1;
        last_marg_info->scale_factor_i = -1;
        last_marg_info->gyr0 = gnss_imu_infos[idx0]->pre_integrations[gnss_imu_infos[idx0]->pre_integrations.size() - 1]->gyrj;
        for (int i = 0; i < (int)last_marg_info->keep_block_addr.size(); i++) {
            if (last_marg_info->keep_block_size[i] == 9) {
                count1++;
                last_marg_info->v_i = i;
                ASSERT(last_marg_info->keep_block_addr[i] == para_speed_bias[idx0]);
            }
            for (int j = idx0; j < SWF_SIZE_IN + idx0 + 1; j++) {
                if (last_marg_info->keep_block_addr[i] == para_pose[j]) {
                    count2++;
                    last_marg_info->p_is[j - idx0] = i;
                }
            }
            if (last_marg_info->keep_block_addr[i] != &scale_factor && last_marg_info->keep_block_addr[i] != &td && last_marg_info->keep_block_size[i] == 1) {
                if (idepth_pointers.find(last_marg_info->keep_block_addr[i]) != idepth_pointers.end())
                    last_marg_info->lambda_is.push_back(i);
                else
                    ASSERT(phase_bias_pointers.find(last_marg_info->keep_block_addr[i]) != phase_bias_pointers.end());
            }

            if (last_marg_info->keep_block_addr[i] == para_extrinsic)
                last_marg_info->extrinsic_i = i;
            if (last_marg_info->keep_block_addr[i] == para_global_pos)
                last_marg_info->R_WI_WC_i = i;
            if (last_marg_info->keep_block_addr[i] == &scale_factor)
                last_marg_info->scale_factor_i = i;



        }

        ASSERT(count1 == 1
               && count2 == (int)last_marg_info->p_is.size()
               && last_marg_info->v_i >= 0
               && (last_marg_info->extrinsic_i >= 0 || !ESTIMATE_EXTRINSIC)
               && (last_marg_info->scale_factor_i >= 0 || !fix_scale)
               && (last_marg_info->R_WI_WC_i >= 0 || !ESTIMATE_R_WI_WC || !ESTIMATE_P_WI_WC));;

        ASSERT(last_marg_info->lambda_is.size() + last_marg_info->p_is.size() + 1 <= last_marg_info->keep_block_addr.size());

    }

}
