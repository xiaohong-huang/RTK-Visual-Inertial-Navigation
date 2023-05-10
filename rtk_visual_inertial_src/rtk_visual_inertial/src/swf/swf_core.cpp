#include "swf.h"
#include "../factor/pose_local_parameterization.h"
#include "../factor/projection_factor.h"
#include "../factor/mag_factor.h"
#include "../factor/gnss_factor.h"
#include "../factor/initial_factor.h"
#include "../factor/mag_factor.h"
#include "../factor/pose0_factor.h"

#define ADDRESIDUAL(parameter_block,residual_define,arg2,sparseset){\
    bool keep=false;\
    if(mode==NormalMode||mode==NormalMode2){\
        residual_define\
        problem.AddResidualBlock(factor, loss_function, parameter_block); \
    }\
    else {\
        std::vector<double*>parameter_block_vector{parameter_block};\
        std::vector<int>dropset;std::vector<int>keepset;\
        for(int vi=0;vi<(int)parameter_block_vector.size();vi++){\
            if(mode==MargeExcludeMode)keep=true;\
            if(MargePoint.find(parameter_block_vector[vi])!=MargePoint.end()){\
                dropset.push_back(vi);\
                if(mode==MargeIncludeMode||mode==MargeIncludeMode2||mode==GNSSMargIncludeMode)keep=true;\
            }else{\
                keepset.push_back(vi);\
            }\
        }\
        if(keep||arg2){\
            if(mode==MargeIncludeMode2){\
                int ks;\
                for(int vi=0;vi<(int)keepset.size();vi++){\
                    for(ks=0;ks<(int)ceres::internal::parameter_head.size();ks++)\
                        if(ceres::internal::parameter_head[ks]==parameter_block_vector[keepset[vi]])break;\
                    if(ks==(int)ceres::internal::parameter_head.size()){\
                        ceres::internal::parameter_head.push_back(parameter_block_vector[keepset[vi]]);\
                    }\
                }\
                residual_define\
                problem.AddResidualBlock(factor, loss_function, parameter_block); \
            }else{\
                if(mode==MargeExcludeMode)dropset=keepset;\
                residual_define\
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(factor, loss_function,parameter_block_vector,dropset,sparseset);\
                marginalization_info->addResidualBlockInfo(residual_block_info);\
            }\
        }\
    }\
}


void SWFOptimization::AddParameter2Problem(ceres::Problem& problem, bool disable_constant) {
    for (int i = 0; i < image_count + rover_count; i++) {
        if (problem.HasParameterBlock(para_pose[i])) {
            ceres::LocalParameterization* local_parameterization = new PoseLocalParameterization();
            problem.AddParameterBlock(para_pose[i], SIZE_POSE, local_parameterization);
        }
    }

    for (int i = 0; i < NUM_OF_CAM; i++) {
        if (problem.HasParameterBlock(para_ex_Pose[i])) {
            ceres::LocalParameterization* local_parameterization = new PoseLocalParameterization();
            problem.AddParameterBlock(para_ex_Pose[i], SIZE_POSE, local_parameterization);
        }
    }

    if (disable_constant)return;

    for (int i = 0; i < NUM_OF_CAM; i++) {
        if (problem.HasParameterBlock(para_ex_Pose[i])) {
            if ((ESTIMATE_EXTRINSIC && Vs[g2f[rover_count - 1]].norm() > 0.2) || open_ex_estimation)open_ex_estimation = 1;
            else problem.SetParameterBlockConstant(para_ex_Pose[i]);
        }
    }


    if ((frame_types[rover_count + image_count - 1] == ImagFrame)) {
        if (problem.HasParameterBlock(para_bmg)) {
            problem.SetParameterBlockConstant(para_bmg);
        }
    }
}





void SWFOptimization::AddGnssResidual(int mode, std::set<double*>MargePoint, MarginalizationInfo* marginalization_info,
                                      ceres::Problem& problem,
                                      ceres::Solver::Options& options, int ir) {

    //gnss residual
    if (USE_GNSS) {
        bool have_base=false;
        ceres::LossFunction* loss_function = NULL;
        mea_t* rover = rovers[ir];
        double globalxyz[3];
        globalxyz[0] = para_pose[g2f[ir]][0] + rover->base_xyz[0];
        globalxyz[1] = para_pose[g2f[ir]][1] + rover->base_xyz[1];
        globalxyz[2] = para_pose[g2f[ir]][2] + rover->base_xyz[2];
        update_azel(globalxyz, rover);

#define A4s &blackvalue
#define B4s InitialBlackFactor*factor=new InitialBlackFactor(1);
        ADDRESIDUAL(A4s, B4s, false, std::vector<int> {})
        if (USE_RTK) {
            for (uint8_t i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                uint8_t sys = d->sys;
                for (uint8_t f = 0; f < NFREQ; f++) {
                    if (d->RTK_L[f] != 0 && rover_count > 1)assert(d->RTK_Npoint[f]);
                    if (!d->RTK_Npoint[f] )continue;
                    if (d->el < AZELMIN)continue;
                    assert(d->el >= AZELMIN);
                    have_base=true;
#define A44 para_pose[g2f[ir]],&(d->RTK_Npoint[f]->value),para_gnss_dt[0]+sys*2+f
#define B44 RTKCarrierPhaseFactor *factor = new RTKCarrierPhaseFactor(d->satellite_pos,d->RTK_L[f]*lams[d->sys][f],lams[d->sys][f],d->el,rover->br_time_diff,pow(d->RTK_Lstd[f]*lams[d->sys][f],2),rover->base_xyz,true,d->sys,f);
                    ADDRESIDUAL(A44, B44, false, std::vector<int> {})

                }
            }
        }
        if (USE_RTD) {
            for (uint8_t i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                uint8_t sys = d->sys;
                assert(sys == 0 || sys == 1 || sys == 2);
                for (uint8_t f = 0; f < NFREQ; f++) {
                    if (d->RTK_P[f] == 0.0 || d->SVH != 0  || d->RTK_Pstd[f] > 2)continue;
                    if (d->el < AZELMIN)continue;
                    assert(d->el >= AZELMIN);
                    have_base=true;
#define A54 para_pose[g2f[ir]],para_gnss_dt[0]+sys*2+f
#define B54 RTKPseudorangeFactor *factor = new RTKPseudorangeFactor(d->satellite_pos,d->RTK_P[f],  d->el,rover->br_time_diff,pow(d->RTK_Pstd[f],2), rover->base_xyz);
                    ADDRESIDUAL(A54, B54, false, std::vector<int> {})
                }
            }
        }


        for (uint8_t i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            if (d->SVH != 0)continue;
            if (d->el < AZELMIN)continue;
            if (d->SPP_P[0] != 0)assert(d->el >= AZELMIN);

            if (d->SPP_P[0] != 0.0 && d->SPP_Pstd[0] < 2 && !have_base) {
                assert(!USE_RTD);
                double sin_el = sin(d->el);
                //0.125,0.7,0.35,1 is described in Global Positioning System：Theory and Applications(Volume I),pp481
                double istd = sin_el * sin_el / sqrt(SQR(d->SPP_Pstd[0]) + (d->ion_var * 0.125 * 0.125 + d->trop_var * 0.7 * 0.7 + d->sat_var * 0.35 * 0.35 + 1));

                assert(istd > 0);
                //because we use the dogleg strategy, we set the higher istd of GNSS measurements in the start time, so that the states can be converge to good results faster.
                if (rover_count_accumulate - rover_count + ir < 100)istd *= 10;

#define A8 para_pose[g2f[ir]],para_gnss_dt[0]+6+d->sys*2+0
#define B8 SppPseudorangeFactor *factor = new SppPseudorangeFactor(d->satellite_pos,d->SPP_P[0],istd,rover->base_xyz);
                ADDRESIDUAL(A8, B8, false, std::vector<int> {})


            }
            if (USE_SPP_PHASE && d->SPP_L[0] != 0.0 && d->SPP_Npoint[0]) {
                double lam = lams[d->sys][0];
                double istd = 0;
                double sin_el = sin(d->el);
                istd = sin_el * sin_el / sqrt(SQR(d->SPP_Lstd[0] * lam) + (d->ion_var * 0.125 * 0.125 + d->trop_var * 0.7 * 0.7 + d->sat_var * 0.35 * 0.35));

                assert(istd > 0);
#define A8s para_pose[g2f[ir]],para_gnss_dt[0]+6+d->sys*2+0,&(d->SPP_Npoint[0]->value)
#define B8s SppCarrierPhaseFactor *factor = new SppCarrierPhaseFactor(d->satellite_pos,d->SPP_L[0]*lam ,istd,rover->base_xyz,lam);
                ADDRESIDUAL(A8s, B8s, false, std::vector<int> {})
            }

            if (USE_SPP_CORRECTION && d->SPP_P0[0] != 0.0 && d->SPP_Npoint_PCottections[0]) {
                double lam = lams[d->sys][0];
                double istd = 0;
                double sin_el = sin(d->el);

                assert(d->SPP_P[0] != 0);
                istd = sin_el * sin_el / sqrt(SQR(d->SPP_Pstd[0]) + (d->ion_var * 0.125 * 0.125 + d->trop_var * 0.7 * 0.7 + d->sat_var * 0.35 * 0.35));

                assert(istd > 0);
#define A8ss para_pose[g2f[ir]],para_gnss_dt[0]+6+d->sys*2+0,&(d->SPP_Npoint_PCottections[0]->value)
#define B8ss SppCarrierPhaseFactor *factor = new SppCarrierPhaseFactor(d->satellite_pos,d->SPP_P0[0] ,istd,rover->base_xyz,lam);
                ADDRESIDUAL(A8ss, B8ss, false, std::vector<int> {})
            }


        }
        if (USE_DOPPLER) {
            for (uint8_t i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                if (d->SPP_D[0] == 0.0 || d->SVH != 0)continue;
                if (d->SPP_Dstd[0] > 2)continue;
                if (d->el < AZELMIN)continue;
                double istd = sin(d->el) * sin(d->el) / (d->SPP_Dstd[0] * lams[d->sys][0]);
                assert(istd > 0);

#define A92 para_speed_bias[g2f[ir]],para_gnss_dt[0]+12,para_pose[g2f[ir]]
#define B92 SppDopplerFactor*factor=new SppDopplerFactor(d->satellite_vel,d->satellite_pos,para_pose[0],d->SPP_D[0]*lams[d->sys][0],istd,rover->base_xyz);
                ADDRESIDUAL(A92, B92, false, std::vector<int> {})
            }
        }
    }
}



void SWFOptimization::AddAllResidual(int mode, std::set<double*>MargePoint, MarginalizationInfo* marginalization_info,
                                     ceres::Problem& problem, ceres::Solver::Options& options,
                                     bool usefeature, bool useimu, bool usegnss) {



    ceres::LossFunction* loss_function = NULL;
    ceres::Problem problem2;

    if (mode == NormalMode || mode == NormalMode2 || mode == MargeIncludeMode2) {
        InitialBlackFactor* factor = new InitialBlackFactor(1);
        problem.AddResidualBlock(factor, 0, &blackvalue2);
    }

    if (mode == NormalMode || mode == NormalMode2 || mode == MargeIncludeMode2) {
        for (int i = 0; i < image_count + 1; i++) {
            IMUGNSSBase* IMUGNSSmeasurement = imu_gnss_factor[i];
            int index1, index2;
            if (image_count > 0) {
                if (i == 0) {
                    index1 = 0;
                    index2 = i2f[0];
                } else if (i < image_count) {
                    index1 = i2f[i - 1];
                    index2 = i2f[i];
                } else {
                    assert(i == image_count);
                    index1 = i2f[i - 1];
                    index2 = rover_count + image_count - 1;

                }
            } else {
                index1 = 0;
                index2 = rover_count - 1;
            }
            if (index2 - index1 == 1) {

                if (pre_integrations[index2]->sum_dt > IMUINTEGRATIONTIMEMAX)continue;
#define A10_ para_pose[index1], para_speed_bias[index1], para_pose[index2], para_speed_bias[index2]
#define B10_ IMUFactor* factor = new IMUFactor(pre_integrations[index2]);
                ADDRESIDUAL (A10_, B10_, false, std::vector<int> {})

            } else if (index2 - index1 > 1) {
                if (IMUGNSSmeasurement) {
                    assert(index2 - index1 > 1);
                    std::vector<double*>param{para_pose[index1], para_speed_bias[index1],
                                              para_pose[index2], para_speed_bias[index2]};
                    for (int j = 0; j < (int)(IMUGNSSmeasurement->gnss_phase_biases.size()); j++)
                        param.push_back(IMUGNSSmeasurement->gnss_phase_biases[j]);
                    assert(param.size() == IMUGNSSmeasurement->param.size());
                    for (int ki = 0; ki < (int)param.size(); ki++)assert(param[ki] == IMUGNSSmeasurement->param[ki]);
#define A20 param
#define B20 IMUGNSSFactor*factor=new IMUGNSSFactor(imu_gnss_factor[i]);
                    ADDRESIDUAL (A20, B20, false, std::vector<int> {})
                }
            }

        }
    }


    if (useimu && USE_IMU) {
        if (USE_MAG_CORRECT_YAW) {
            for (int i = 0; i < image_count; i++) {
                if (mags[i2f[i]].norm() > 1e-4 && frame_types[i2f[i]] == ImagFrame) {
#define A18 para_pose[i2f[i]],para_bmg
#define B18 MAGFactor*factor=new MAGFactor(mags[i2f[i]]);
                    ADDRESIDUAL(A18, B18, false, std::vector<int> {})
                }
            }
        }
        for (int j = 1; j < image_count + rover_count; j++) {
            if (pre_integrations[j]->sum_dt > IMUINTEGRATIONTIMEMAX)continue;
            if ((mode == NormalMode || mode == NormalMode2 || mode == MargeIncludeMode2))continue;
#define A10 para_pose[j-1], para_speed_bias[j-1], para_pose[j], para_speed_bias[j]
#define B10 IMUFactor* factor = new IMUFactor(pre_integrations[j]);
            ADDRESIDUAL (A10, B10, false, std::vector<int> {})
        }
    }
    if (usegnss) {
        for (int ir = 0; ir < rover_count; ir++) {
            if (rovers[ir]->marg_info_gnss) {
                if ((mode == NormalMode || mode == NormalMode2 || mode == MargeIncludeMode2) && (g2f[ir] != 0 && g2f[ir] != rover_count + image_count - 1))
                    continue;
#define A2_ rovers[ir]->marg_info_gnss->keep_block_addr
#define B2_ MarginalizationFactor *factor = new MarginalizationFactor(rovers[ir]->marg_info_gnss);
                ADDRESIDUAL(A2_, B2_, false, std::vector<int> {})
            }
        }
    }

    if (usefeature && USE_IMAGE) {

        for (auto& it_per_id : f_manager.feature) {
            it_per_id.used_num = it_per_id.feature_per_frame.size();
            if (it_per_id.used_num < FEATURE_CONTINUE)continue;
            int imu_j = it_per_id.start_frame ;
#if USE_INVERSE_DEPTH
            Vector3d pts_i = it_per_id.feature_per_frame[0].point;
#endif
            for (auto& it_per_frame : it_per_id.feature_per_frame) {

#if USE_INVERSE_DEPTH

                if (it_per_id.start_frame != imu_j) {
                    Vector3d pts_j = it_per_frame.point;
                    ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
#define A132 para_pose[i2f[it_per_id.start_frame]], para_pose[i2f[imu_j]], para_ex_Pose[0], &it_per_id.idepth_
#define B132 ProjectionTwoFrameOneCamFactor *factor = new ProjectionTwoFrameOneCamFactor(pts_i, pts_j);
                    ADDRESIDUAL (A132, B132, false, std::vector<int> {})
                }
                if (USE_STEREO && it_per_frame.is_stereo) {
                    Vector3d pts_j_right = it_per_frame.pointRight;
                    if (it_per_id.start_frame != imu_j) {
                        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
#define A133 para_pose[i2f[it_per_id.start_frame]], para_pose[i2f[imu_j]], para_ex_Pose[0], para_ex_Pose[1], &it_per_id.idepth_
#define B133 ProjectionTwoFrameTwoCamFactor *factor = new ProjectionTwoFrameTwoCamFactor(pts_i, pts_j_right);
                        ADDRESIDUAL (A133, B133, false, std::vector<int> {})
                    } else {
                        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
#define A134 para_ex_Pose[0], para_ex_Pose[1], &it_per_id.idepth_
#define B134 ProjectionOneFrameTwoCamFactor *factor= new ProjectionOneFrameTwoCamFactor(pts_i, pts_j_right);
                        ADDRESIDUAL (A134, B134, false, std::vector<int> {})
                    }
                }
#else
                ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
#define A13 para_pose[i2f[imu_j]], para_ex_Pose[0], it_per_id.ptsInWorld.data()
#define B13 projection_factor *factor = new projection_factor(it_per_frame.point);
                ADDRESIDUAL (A13, B13, false, std::vector<int> {})

                if (USE_STEREO && it_per_frame.is_stereo) {
                    ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
#define A14 para_pose[i2f[imu_j]],  para_ex_Pose[1], it_per_id.ptsInWorld.data()
#define B14 projection_factor *factor = new projection_factor(it_per_frame.pointRight);
                    ADDRESIDUAL (A14, B14, false, std::vector<int> {})

                }
#endif
                imu_j++;
            }
        }
    }



    if (last_marg_info && mode != GNSSMargIncludeMode) { //&&(frame_types[0]!=GnssFrame||D0616==0||!useimu)
#define A2 last_marg_info->keep_block_addr
#define B2 MarginalizationFactor *factor = new MarginalizationFactor(last_marg_info);
        std::vector<double*>a = std::vector<double*> {last_marg_info->keep_block_addr};
        std::cout << a.size() << "," << last_marg_info->keep_block_addr.size() << std::endl;
        ADDRESIDUAL(A2, B2, mode == MargeIncludeMode || mode == MargeIncludeMode2, std::vector<int> {})
    }


    if (mode == NormalMode || mode == NormalMode2 || mode == MargeIncludeMode2) {
        AddParameter2Problem(problem, mode == MargeIncludeMode2);
    }
    if (mode == MargeExcludeMode) {

        marginalization_info->marginalize_exclude();
        std::cout << "exclude:" << marginalization_info->n << "," << marginalization_info->m << std::endl;

    } else if (mode == MargeIncludeMode) {

        marginalization_info->marginalize(false, true);
        if (marginalization_info->n == 0) {
            if (last_marg_info) {
                delete last_marg_info;
                if (USE_GLOBAL_OPTIMIZATION)my_problem.RemoveResidualBlock(marg_residual_block_id);
            }
            last_marg_info = 0;
            printf("warning!!no last info\r\n");
            return;
        }
        if (!marginalization_info->n)printf("warning!!no last info\r\n");
        marginalization_info->getParameterBlocks();
        if (last_marg_info) {
            if (USE_GLOBAL_OPTIMIZATION)my_problem.RemoveResidualBlock(marg_residual_block_id);
            delete last_marg_info;
        }
        last_marg_info = marginalization_info;
        if (USE_GLOBAL_OPTIMIZATION) {
            MarginalizationFactor* factor = new MarginalizationFactor(last_marg_info);
            marg_residual_block_id = my_problem.AddResidualBlock(factor, 0, last_marg_info->keep_block_addr);
        }

    } else if (mode == NormalMode) {

        options.linear_solver_type = ceres::DENSE_SCHUR;
        options.trust_region_strategy_type = ceres::DOGLEG;
        options.max_num_iterations = 8;
        options.num_threads = 4;
        options.jacobi_scaling = 0;

        ceres::Solver::Summary summary;
        static double t, count;
        MyOrdering(problem, options);
        TicToc t_marg;
        ceres::Solve(options, &problem, &summary);
        double consumption_time = summary.minimizer_time_in_seconds * 1000 / summary.num_successful_steps * 8;
        t += consumption_time;
        count += 1;

        std::cout << "solving time:" << consumption_time << "," << t / count << std::endl;

    } else if (mode == NormalMode2) {

        options.linear_solver_type = ceres::DENSE_SCHUR;
        options.trust_region_strategy_type = ceres::DOGLEG;
        // options.initial_trust_region_radius=options.max_trust_region_radius=MAX_TRUST_REGION_RADIUS;
        options.max_num_iterations = 1;
        options.num_threads = 8;
        options.jacobi_scaling = false;
        ceres::Solver::Summary summary;
        TicToc t_marg;
        ceres::internal::is_optimize = false;
        MyOrdering(problem, options);
        ceres::Solve(options, &problem, &summary);
        ceres::internal::is_optimize = true;
        std::cout << "solving time:" << t_marg.toc() << std::endl;

    } else if (mode == MargeIncludeMode2) {

        if (ceres::internal::parameter_head.size() == 0) {
            if (last_marg_info) {
                delete last_marg_info;
                if (USE_GLOBAL_OPTIMIZATION) {
                    my_problem.RemoveResidualBlock(marg_residual_block_id);
                }
            }
            last_marg_info = 0;
            printf("warning!!no last info\r\n");
            return;
        } else {
            options.linear_solver_type = ceres::DENSE_SCHUR;
            options.trust_region_strategy_type = ceres::DOGLEG;
            // options.initial_trust_region_radius=options.max_trust_region_radius=MAX_TRUST_REGION_RADIUS;
            options.max_num_iterations = 1;
            options.num_threads = 8;
            options.jacobi_scaling = false;
            ceres::Solver::Summary summary;
            TicToc t_marg;
            ceres::internal::is_optimize = false;

            MyOrdering(problem, options);
            ceres::Solve(options, &problem, &summary);
            std::cout << "solving time:" << t_marg.toc() << std::endl;
            UpdateSchur(problem);
            ceres::internal::is_optimize = true;
            marginalization_info->setmarginalizeinfo(parameter_block_addr, parameter_block_global_size, A, b, true);
            marginalization_info->getParameterBlocks();
            if (last_marg_info) {
                delete last_marg_info;
                if (USE_GLOBAL_OPTIMIZATION) {
                    my_problem.RemoveResidualBlock(marg_residual_block_id);
                }
            }
            last_marg_info = marginalization_info;
            if (USE_GLOBAL_OPTIMIZATION) {
                MarginalizationFactor* factor = new MarginalizationFactor(last_marg_info);
                marg_residual_block_id = my_problem.AddResidualBlock(factor, 0, last_marg_info->keep_block_addr);
            }
        }

    }

}


void SWFOptimization::InitializeSqrtInfo() {


    Vector2Double();
    MarginalizationInfo* marginalization_info = new MarginalizationInfo();
    int mode = MargeExcludeMode;
    int i = 0;



#define A11 para_pose[i]
#define B11 InitialPoseFactor*factor=new InitialPoseFactor(Ps[i], Quaterniond(Rs[i]),sqrt_info_pose);
#define A12 para_speed_bias[i]
#define B12 InitialBiasFactor*factor=new InitialBiasFactor(Vs[i],Bas[i], Bgs[i],sqrt_info_bias);
#define A23 para_bmg
#define B23 InitialFactor33*factor=new InitialFactor33(Eigen::Vector3d(para_bmg),sqrt_info_pose4);

    std::set<double*>MargePoint;
    ceres::Problem problem;
    ceres::LossFunction* loss_function = NULL;
    Eigen::Matrix<double, 6, 6>sqrt_info_pose;
    sqrt_info_pose.setZero();
    sqrt_info_pose.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * ((rover_count > 0) ? 1e-3 : 2e2);
    sqrt_info_pose.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 2e2;
    MargePoint.insert(A11);
    ADDRESIDUAL (A11, B11, false, std::vector<int> {})

    {
        Eigen::Matrix<double, 9, 9>sqrt_info_bias;
        sqrt_info_bias.setZero();
        sqrt_info_bias.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * ((rover_count > 0) ? 1e-3 : 1e1);
        sqrt_info_bias.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 1;
        sqrt_info_bias.block<3, 3>(6, 6) = Eigen::Matrix<double, 3, 3>::Identity() * 1e1;
        MargePoint.insert(A12);
        ADDRESIDUAL (A12, B12, false, std::vector<int> {})
    }

    if (ESTIMATE_EXTRINSIC) {
        Eigen::Matrix<double, 6, 6>sqrt_info_pose2;
        sqrt_info_pose2.setZero();
        sqrt_info_pose2.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 1e3;
        sqrt_info_pose2.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 1e3;
        for (int i = 0; i < NUM_OF_CAM; i++) {
#define A16 para_ex_Pose[i]
#define B16 InitialPoseFactor*factor=new InitialPoseFactor(tic[i], Quaterniond(ric[i]),sqrt_info_pose2);
            MargePoint.insert(A16);
            ADDRESIDUAL (A16, B16, false, std::vector<int> {})
        }
    }


    if (USE_MAG_CORRECT_YAW) {
        Eigen::Matrix<double, 3, 3>sqrt_info_pose4;
        sqrt_info_pose4.setZero();
        sqrt_info_pose4.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 1;
        MargePoint.insert(A23);
        ADDRESIDUAL (A23, B23, false, std::vector<int> {})
    }



    marginalization_info->marginalize(true, true);
    marginalization_info->getParameterBlocks();
    if (last_marg_info) {
        delete last_marg_info;
        if (USE_GLOBAL_OPTIMIZATION) {
            my_problem.RemoveResidualBlock(marg_residual_block_id);
        }
    }
    last_marg_info = marginalization_info;

    if (USE_GLOBAL_OPTIMIZATION) {
        MarginalizationFactor* factor = new MarginalizationFactor(last_marg_info);
        marg_residual_block_id = my_problem.AddResidualBlock(factor, 0, last_marg_info->keep_block_addr);
        if (!my_problem.HasParameterBlock(&blackvalue2)) {
            InitialBlackFactor* factor = new InitialBlackFactor(1);
            my_problem.AddResidualBlock(factor, 0, &blackvalue2);
        }

        my_problem.AddParameterBlock(para_pose[0], SIZE_POSE, new PoseLocalParameterization());
        for (int i = 0; i < NUM_OF_CAM; i++) {
            my_problem.AddParameterBlock(para_ex_Pose[i], SIZE_POSE, new PoseLocalParameterization());
        }
    }
}





void SWFOptimization::MiddleMargGnssFrame() {
#define GNSS_MAX_NUM 30
    for (int IMUGNSSindex = 0; IMUGNSSindex < FEATURE_WINDOW_SIZE + 2; IMUGNSSindex++) {
        if (!imu_gnss_factor[IMUGNSSindex])continue;
        if (imu_gnss_factor[IMUGNSSindex]->gnss_speed_bias.size() <= GNSS_MAX_NUM)continue;
        Eigen::VectorXd V1, V2;

        TicToc tic;
        int gnss_num = imu_gnss_factor[IMUGNSSindex]->gnss_speed_bias.size();
        int index1, index2;
        if (image_count > 0) {
            if (IMUGNSSindex == 0) {
                index1 = 0;
                index2 = i2f[IMUGNSSindex];
            } else if (IMUGNSSindex < image_count) {
                index1 = i2f[IMUGNSSindex - 1];
                index2 = i2f[IMUGNSSindex];
            } else {
                assert(IMUGNSSindex == image_count);
                index1 = i2f[IMUGNSSindex - 1];
                index2 = rover_count + image_count - 1;
            }
        } else {
            index1 = 0;
            index2 = rover_count - 1;
        }
        index1 += 1;

        int marg_num = gnss_num - GNSS_MAX_NUM;
        // int marg_index0=GNSS_MAX_NUM/2;
        std::set<int>margeindex;
        for (int i = 0; i < marg_num; i++) {
            int index = GNSS_MAX_NUM / 2 + i + index1;
            assert(index < index2);
            margeindex.insert(index);
        }//debug

        MarginalizationInfo* gnss_middle_marginfo = MargGNSSFrames(margeindex, imu_gnss_factor[IMUGNSSindex]);

        if (imu_gnss_factor[IMUGNSSindex]->gnss_middle_marginfo) {
            MarginalizationFactor* factor = new MarginalizationFactor(imu_gnss_factor[IMUGNSSindex]->gnss_middle_marginfo);
            ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(factor, 0, \
                    imu_gnss_factor[IMUGNSSindex]->gnss_middle_marginfo->keep_block_addr, std::vector<int>({}), std::vector<int>({}));
            gnss_middle_marginfo->addResidualBlockInfo(residual_block_info);
        }
        PhaseBiasSaveAndReset();
        gnss_middle_marginfo->marginalize(true, true);
        gnss_middle_marginfo->getParameterBlocks();
        PhaseBiasRestore();
        int a = 0, b = 0, c = 0;
        for (int i = 0; i < (int)gnss_middle_marginfo->keep_block_size.size(); i++) {
            if (gnss_middle_marginfo->keep_block_size[i] == 9)a++;
            else if (gnss_middle_marginfo->keep_block_size[i] == 7)b++;
            else if (gnss_middle_marginfo->keep_block_size[i] == 3)c++;
        }
        assert(a == 2 && b == 2 && c == 0);
        pre_integrations[GNSS_MAX_NUM / 2 + index1]->sum_dt = IMUINTEGRATIONTIMEMAX + 100000;
        int index = GNSS_MAX_NUM / 2 + index1;
        for (int i = 0; i < marg_num; i++) {
            SlideWindowFrame(index, image_count + rover_count, false);
            free(rovers[f2g[index]]);//debug
            for (int ri = f2g[index]; ri < rover_count; ri++)rovers[ri] = rovers[ri + 1];
            rover_count--;
            UpdateVisualGnssIndex();
        }

        PhaseBiasSaveAndReset();
        rovers[f2g[GNSS_MAX_NUM / 2 + index1]]->marg_info_gnss->ResetLinearizationPoint(rovers[f2g[GNSS_MAX_NUM / 2 + index1]]->marg_info_gnss->keep_block_addr);
        rovers[f2g[GNSS_MAX_NUM / 2 + index1 - 1]]->marg_info_gnss->ResetLinearizationPoint(rovers[f2g[GNSS_MAX_NUM / 2 + index1 - 1]]->marg_info_gnss->keep_block_addr);
        PhaseBiasRestore();

        assert(pre_integrations[GNSS_MAX_NUM / 2 + index1]->sum_dt > IMUINTEGRATIONTIMEMAX + 100000);
        ResetImuGnssFactor(IMUGNSSindex, gnss_middle_marginfo);

        std::cout << "marge middle gnss spend time:" << tic.toc() << std::endl;



    }
}
