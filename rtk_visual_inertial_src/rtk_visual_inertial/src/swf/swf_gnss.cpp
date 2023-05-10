#include "swf.h"
#include "ceres/residual_block.h"
#include "../factor/gnss_factor.h"




void SWFOptimization::InputGnss(mea_t* rovernew2) {



    static int count = 0;
    if (count++ < 2)return;
    mea_t* rovernew = (mea_t*)malloc(sizeof(mea_t));
    memcpy(rovernew, rovernew2, sizeof(mea_t));

    rover_buf.push_back(rovernew);

}





void SWFOptimization::UpdateSchur(ceres::Problem& problem) {
    assert(!ceres::internal::is_optimize);
    TicToc tic;
    A = ceres::ConstMatrixRef(ceres::internal::lhs_out, ceres::internal::hs_row, ceres::internal::hs_row).selfadjointView<Eigen::Upper>();
    b = Eigen::Map<Eigen::VectorXd>(ceres::internal::rhs_out, ceres::internal::hs_row);
    int parameter_head_all_size = 0;
    parameter_block_addr.clear();
    parameter_block_global_size.clear();

    for (size_t i = 0; i < ceres::internal::parameter_head.size(); i++) {
        int size = problem.ParameterBlockSize(ceres::internal::parameter_head[i]);
        if (size == 7)parameter_head_all_size += 6;
        else parameter_head_all_size += size;
        parameter_block_addr.push_back(ceres::internal::parameter_head[i]);
        parameter_block_global_size.push_back(size);
    }

    int m = ceres::internal::hs_row - parameter_head_all_size;
    int n = parameter_head_all_size;

    if (m) {
        Eigen::MatrixXd Amm = A.block(0, 0, m, m);
        Eigen::VectorXd bmm = b.segment(0, m);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(Amm);
        Eigen::MatrixXd Amm_inv = saes2.eigenvectors() * Eigen::VectorXd((saes2.eigenvalues().array() > 1e-8).select(saes2.eigenvalues().array().inverse(), 0)).asDiagonal() * saes2.eigenvectors().transpose();
        Eigen::MatrixXd Anm = A.block(m, 0, n, m);
        Eigen::MatrixXd Amn = A.block(0, m, m, n);
        Eigen::MatrixXd Ann = A.block(m, m, n, n);
        Eigen::VectorXd bnn = b.segment(m, n);
        A = Ann - Anm * Amm_inv * Amn;
        b = bnn - Anm * Amm_inv * bmm;
    }

    ceres::internal::parameter_head.clear();
    std::cout << "update schur time:" << tic.toc() << std::endl;

}



void SWFOptimization::UpdateSchurHessianOnly(ceres::Problem& problem) {
    TicToc tic;

    int parameter_head_all_size = 0;
    parameter_block_addr.clear();
    parameter_block_global_size.clear();
    assert(ceres::internal::is_optimize);
    for (size_t i = 0; i < ceres::internal::parameter_head.size(); i++) {
        int size = problem.ParameterBlockSize(ceres::internal::parameter_head[i]);
        if (size == 7)parameter_head_all_size += 6;
        else parameter_head_all_size += size;
        parameter_block_addr.push_back(ceres::internal::parameter_head[i]);
        parameter_block_global_size.push_back(size);
    }

    int m = ceres::internal::hs_row - parameter_head_all_size;
    int n = parameter_head_all_size;

    if (m) {

        ceres::Matrix A2 = Eigen::Map<ceres::Matrix>(ceres::internal::lhs_out2, ceres::internal::hs_row, ceres::internal::hs_row);
        ceres::Matrix A3 = A2.block(m, m, n, n);
        A = A3 * A3.transpose();
        // std::cout<<A3*A3.transpose()<<std::endl;
    }

    ceres::internal::parameter_head.clear();
    std::cout << "UpdateSchurHessianOnly time:" << tic.toc() << std::endl;

}

void SWFOptimization::UpdateNParameterHead() {
    {
        if (!USE_RTK)return;
        ceres::internal::parameter_head.clear();
        for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
            for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
                it1->use = false;
            }
        }
        for (int ir = 0; ir < rover_count; ir++) {
            mea_t* rover = rovers[ir];
            for (uint8_t i = 0; i < rover->obs_count; i++) {
                ObsMea* d = rover->obs_data + i;
                for (uint8_t f = 0; f < NFREQ; f++) {
                    if (d->SVH == 0)
                        assert(d->el != 0);
                    if (d->RTK_Lstd[f] <= 0.004 * LTHRE && d->RTK_Npoint[f] &&
                            d->RTK_Npoint[f]->continue_count > CARRIER_PHASE_CONTINUE_THRESHOLD &&
                            ((d->half_flag[f]) & 8) && ((d->half_flag[f]) & 2) && d->el >= AZELMIN) { //
                        if (!d->RTK_Npoint[f]->use) {
                            ceres::internal::parameter_head.push_back(&(d->RTK_Npoint[f]->value));
                            d->RTK_Npoint[f]->use = true;
                            if (USE_GLOBAL_OPTIMIZATION)assert(my_problem.HasParameterBlock(&(d->RTK_Npoint[f]->value)));
                        }
                    }
                }
            }
        }

    }
}

void SWFOptimization::IntegerSolve() {
    TicToc tic;

    if (!USE_RTK || !need_Nresolve || !USE_N_RESOLVE)return;
    Vector2Double();
    if (!USE_DIRECT_N_RESOLVE) {

        UpdateNParameterHead();
        if (USE_GLOBAL_OPTIMIZATION) {
            my_options.max_num_iterations = 1;
            my_options.jacobi_scaling = false;
            ceres::internal::is_optimize = false;
            if (USE_MAG_CORRECT_YAW)my_problem.SetParameterBlockVariable(para_bmg);
            for (int i = 0; i < NUM_OF_CAM; i++) {
                my_problem.SetParameterBlockVariable(para_ex_Pose[i]);
            }
            ceres::Solver::Summary summary;
            TicToc t_marg;
            MyOrdering(my_problem, my_options);
            // my_options.initial_trust_region_radius=my_options.max_trust_region_radius=MAX_TRUST_REGION_RADIUS;
            ceres::Solve(my_options, &my_problem, &summary);
            std::cout << "solving time:" << t_marg.toc() << std::endl;

            my_options.max_num_iterations = MAX_NUM_ITERATIONS;
            my_options.jacobi_scaling = 0;
            // my_options.initial_trust_region_radius=my_options.max_trust_region_radius=MAX_TRUST_REGION_RADIUS;
            UpdateSchur(my_problem);
            ceres::internal::is_optimize = true;

        } else {
            ceres::Problem problem;
            ceres::Solver::Options options;
            AddAllResidual(NormalMode2, std::set<double*> {}, nullptr, problem, options, true, true, true);
            UpdateSchur(problem);
        }
    }

    LambdaSearch();

    Double2Vector();
    need_Nresolve = false;
    std::cout << "LAMBDASearch:" << tic.toc() << std::endl;

}

std::set<int>set_sat;

void SWFOptimization::GnssProcess(mea_t* rover) {

    {
        double globalxyz[3];
        globalxyz[0] = Ps[g2f[rover_count - 1]].x() + rover->base_xyz[0];
        globalxyz[1] = Ps[g2f[rover_count - 1]].y() + rover->base_xyz[1];
        globalxyz[2] = Ps[g2f[rover_count - 1]].z() + rover->base_xyz[2];
        update_azel(globalxyz, rover);
    }
    rover_count_accumulate++;
    rtk_fix = false;
    marg_flag = MargRoverOld;
    rovers[rover_count - 1] = rover;

    if (!init_gnss) {
        {
            double pos[3];
            Matrix3d Rwwg;
            std::cout << "ANCHOR_POINT:" << ANCHOR_POINT.transpose() << std::endl;
            ecef2pos(ANCHOR_POINT.data(), pos);
            xyz2enu(pos, Rwwg.data());//transform vector from ecef to enu
            InitRwgw = Rwwg.transpose();
        }
        pub_init = false;
        init_gnss = true;
        solver_flag = NonLinear;
        base_pos = Eigen::Map<Eigen::Vector3d> (rovers[0]->base_xyz);
        {
            Vector2Double();
            ceres::Problem problem;
            ceres::Solver::Options options;
            AddGnssResidual(NormalMode, std::set<double*> {}, 0, problem, options, rover_count - 1);
            options.linear_solver_type = ceres::DENSE_SCHUR;
            // options.trust_region_strategy_type = ceres::DOGLEG;
            options.initial_trust_region_radius = options.max_trust_region_radius = MAX_TRUST_REGION_RADIUS;
            options.max_num_iterations = 20;
            options.num_threads = 8;
            TicToc t_marg;
            ceres::Solver::Summary summary;
            ceres::Solve(options, &problem, &summary);
        }


        Vector3d        Psold[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
        Vector3d        Vsold[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
        Matrix3d        Rsold[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
        for (int i = 0; i < rover_count + image_count; i++) {
            Psold[i] = Ps[i];
            Vsold[i] = Vs[i];
            Rsold[i] = Rs[i];
            Ps[i][0] = para_pose[g2f[rover_count - 1]][0];
            Ps[i][1] = para_pose[g2f[rover_count - 1]][1];
            Ps[i][2] = para_pose[g2f[rover_count - 1]][2];
            Vs[i][0] = para_speed_bias[g2f[rover_count - 1]][0];
            Vs[i][1] = para_speed_bias[g2f[rover_count - 1]][1];
            Vs[i][2] = para_speed_bias[g2f[rover_count - 1]][2];
        }
        {
            double pos[3];
            Matrix3d Rwwg;
            Eigen::Vector3d xyz = base_pos + Ps[rover_count / 2];
            ecef2pos(xyz.data(), pos);
            xyz2enu(pos, Rwwg.data());//transform vector from ecef to enu
            Rwgw = Rwwg.transpose();

        }
        InitializePos(Rwgw);
        {
            Eigen::Vector3d newP = Ps[g2f[rover_count - 1]];
            Eigen::Matrix3d newR = Rs[g2f[rover_count - 1]];
            Eigen::Vector3d oldP = Psold[g2f[rover_count - 1]];
            Eigen::Matrix3d oldR = Rsold[g2f[rover_count - 1]];
            std::cout << newR* oldR.transpose() << std::endl;
            for (int i = 0; i < rover_count + image_count; i++) {
                Ps[i] = newR * oldR.transpose() * (Psold[i] - oldP) + newP;
                Vs[i] = newR * oldR.transpose() * Vsold[i];
                Rs[i] = newR * oldR.transpose() * Rsold[i];
            }

        }
        InitializeSqrtInfo();


    }

    GnssPreprocess(rovers[rover_count - 1]);

}


void SWFOptimization::GnssPreprocess(mea_t* data) {



    TicToc t_initialize;
    uint8_t i;
    ObsMea* d;
    if (USE_SPP_CORRECTION) {
        for (i = 0; i < data->obs_count; i++) {
            d = data->obs_data + i;
            if (d->SPP_P[0] != 0) {
                d->SPP_P0[0] = d->SPP_P[0];
                if (pseudorange_correction_variables[d->sat * 2 + 0].size()) {
                    auto it = pseudorange_correction_variables[d->sat * 2 + 0].end();
                    it--;
                    it->last_update_time = data->ros_time;
                    assert(pseudorange_correction_variables[d->sat * 2 + 0].size() == 1);
                    if (it->continue_count > EstimatePcorrectionPerio) {
                        d->SPP_P0[0] = 0;
                        d->SPP_P[0] += it->value * lams[d->sys][0];
                    } else {
                        std::cout << (int)(d->sat) << " is in MyOptimization for SPP corrections" << it->continue_count << "\r\n";
                    }
                }

            } else {
                d->SPP_P0[0] = 0;
            }
        }
    }


    for (i = 0; i < data->obs_count; i++) {
        d = data->obs_data + i;
        if (d->SVH)continue;
        for (uint f = 0; f < NFREQ; f++) {
            d->RTK_Npoint[f] = nullptr;
            d->SPP_Npoint[f] = nullptr;
            d->SPP_Npoint_PCottections[f] = nullptr;

            if (d->RTK_L[f] != 0) {
                if (rtk_phase_bias_variables[d->sat * 2 + f].size()) {
                    auto it = rtk_phase_bias_variables[d->sat * 2 + f].end();
                    it--;
                    if (data->ros_time - it->last_update_time < 10)
                        d->RTK_Npoint[f] = &(*it);
                }
            }
            if (d->SPP_L[f] != 0) {
                if (spp_phase_bias_variables[d->sat * 2 + f].size()) {
                    auto it = spp_phase_bias_variables[d->sat * 2 + f].end();
                    it--;
                    if (data->ros_time - it->last_update_time < 10)
                        d->SPP_Npoint[f] = &(*it);
                }
            }
            if (d-> SPP_P0[f] != 0) {
                if (pseudorange_correction_variables[d->sat * 2 + f].size()) {
                    auto it = pseudorange_correction_variables[d->sat * 2 + f].end();
                    it--;
                    if (data->ros_time - it->last_update_time < 10)
                        d->SPP_Npoint_PCottections[f] = &(*it);
                }
            }
        }
    }


    Vector2Double();
    int ir = rover_count - 1;
    mea_t* rover = data;

    std::vector<double>error1_rtk(data->obs_count * 2, 0);
    std::vector<double>error2_rtk[6];
    double median_error_rtk[6];

    std::vector<double>error1_spp(data->obs_count * 2, 0);
    std::vector<double>error2_spp[6];
    double median_error_spp[6];

    for (i = 0; i < data->obs_count; i++) {
        d = data->obs_data + i;
        if (d->SVH)continue;
        double* lam = lams[d->sys];
        uint8_t sys = d->sys;
        for (uint f = 0; f < NFREQ; f++) {
            if (d->el < AZELMIN) {
                d->RTK_L[f] = d->SPP_L[f] = d->SPP_P0[f] = 0;
            }
            if (d->RTK_Npoint[f]) {
                RTKCarrierPhaseFactor* factor = new RTKCarrierPhaseFactor(d->satellite_pos, d->RTK_L[f]*lam[f], lam[f], d->el, 0, 0, rover->base_xyz, false, d->sys, f);
                double residuals;
                std::vector<double*>parameter_blocks = std::vector<double*> {para_pose[g2f[ir]], &(d->RTK_Npoint[f]->value), para_gnss_dt[0] + sys * 2 + f};
                factor->Evaluate(parameter_blocks.data(), &residuals, 0);
                error1_rtk[i * 2 + f] = residuals;
                if (d->RTK_Npoint[f]->LLI == d->RTK_LLI[f]) {
                    error2_rtk[d->sys * 2 + f].push_back(residuals);
                }
                delete factor;
            }

            if (d->SPP_Npoint[f]) {
                RTKCarrierPhaseFactor* factor = new RTKCarrierPhaseFactor(d->satellite_pos, d->SPP_L[f]*lam[f], lam[f], d->el, 0, 0, rover->base_xyz, false, d->sys, f);
                double residuals;
                std::vector<double*>parameter_blocks = std::vector<double*> {
                    para_pose[g2f[ir]], &(d->SPP_Npoint[f]->value), para_gnss_dt[0] + 6 + d->sys * 2 + 0
                };
                factor->Evaluate(parameter_blocks.data(), &residuals, 0);
                error1_spp[i * 2 + f] = residuals;
                if (d->SPP_Npoint[f]->LLI == d->SPP_LLI[f]) {
                    error2_spp[d->sys * 2 + f].push_back(residuals);
                }
                delete factor;
            }
        }
    }
    for (int sys = 0; sys < 3; sys++) {
        for (int f = 0; f < NFREQ; f++) {
            if (error2_rtk[sys * 2 + f].size()) {
                std::sort(error2_rtk[sys * 2 + f].begin(), error2_rtk[sys * 2 + f].end());
                median_error_rtk[sys * 2 + f] = error2_rtk[sys * 2 + f][error2_rtk[sys * 2 + f].size() / 2];
                // if(error2_rtk[sys*2+f].size()<=4)median_error_rtk[sys*2+f]+=10;
            }
            if (error2_spp[sys * 2 + f].size()) {
                std::sort(error2_spp[sys * 2 + f].begin(), error2_spp[sys * 2 + f].end());
                median_error_spp[sys * 2 + f] = error2_spp[sys * 2 + f][error2_spp[sys * 2 + f].size() / 2];
                // if(error2_spp[sys*2+f].size()<=3)median_error_spp[sys*2+f]+=10;
            }
        }
    }
    if (not_fix_count > Phase_ALL_RESET_COUNT)std::cout << "                                                                               not_fix_count reach threshold Phase_ALL_RESET_COUNT\r\n";

    for (i = 0; i < data->obs_count; i++) {
        d = data->obs_data + i;
        if (d->SVH)continue;
        double* lam = lams[d->sys];
        uint8_t sys = d->sys;
        for (uint f = 0; f < NFREQ; f++) {
            if (f == 1)assert(d->RTK_L[f] == 0 && d->SPP_L[f] == 0);
            bool condition3 = false;
            bool condition4 = false;
            if (d->RTK_L[f] != 0) {
                if (USE_IMU && USE_RTK && solver_flag == NonLinear && rover_count > 1 && d->RTK_Npoint[f] && d->RTK_Npoint[f]->LLI == d->RTK_LLI[f]) {
                    double residuals = error1_rtk[i * 2 + f];
                    if (fabs(residuals - median_error_rtk[sys * 2 + f]) > lam[f] / 2) {
                        condition3 = true;
                        std::cout << "residual_rtk:" << residuals << "," << median_error_rtk[sys * 2 + f] << "," << d->RTK_Npoint[f]->value << "," << (int)(d->sat) << "," << (d->RTK_Lstd[f]) << "," << (int)(d->RTK_Npoint[f]->continue_count) << "," << (int)(d->el / PI * 180) << std::endl;
                    }
                }
            }

            if (d->SPP_L[f] != 0) {
                if (USE_IMU && USE_SPP_PHASE && solver_flag == NonLinear && rover_count > 1 && d->SPP_Npoint[f] && d->SPP_Npoint[f]->LLI == d->SPP_LLI[f]) {
                    double residuals = error1_spp[i * 2 + f];
                    if (abs((d->SPP_L[f] + d->SPP_Npoint[f]->value)*lam[f] - d->SPP_P[f])*sin(d->el)*sin(d->el) > 10) {
                        condition4 = true;
                        std::cout << "residual_spp2:" << (d->SPP_L[f] + d->SPP_Npoint[f]->value)*lam[f] << "," << d->SPP_P[f] << ","
                                  << (d->SPP_L[f] + d->SPP_Npoint[f]->value)*lam[f] - d->SPP_P[f] << "," << d->SPP_Pstd[0] << ","
                                  << (int)d->sat << "," << d->SPP_Npoint[f]->continue_count << "," << std::endl;
                    }
                    if (fabs(residuals - median_error_spp[sys * 2 + f]) > lam[f]) {
                        condition4 = true;
                        std::cout << "residual_spp:" << residuals << "," << median_error_spp[sys * 2 + f] << "," << d->SPP_Npoint[f]->value << "," << (int)(d->sat) << "," << (d->SPP_Lstd[f]*lams[d->sys][f]) << "," << (int)(d->SPP_Npoint[f]->continue_count) << "," << (int)(d->el / PI * 180) << std::endl;
                    }
                }
            }

            if (d->RTK_L[f] != 0) {
                if ((!d->RTK_Npoint[f]) || (d->RTK_Npoint[f]->LLI != d->RTK_LLI[f]) || condition3 || not_fix_count > Phase_ALL_RESET_COUNT) {
                    PBtype n;
                    n.sys = d->sys;
                    n.f = f;
                    n.value = 0;
                    n.continue_count = 0;
                    rtk_phase_bias_variables[d->sat * 2 + f].push_back(n);
                    auto it = rtk_phase_bias_variables[d->sat * 2 + f].end();
                    it--;
                    d->RTK_Npoint[f] = &(*it);
                    d->RTK_Npoint[f]->LLI = d->RTK_LLI[f];
                    d->RTK_Npoint[f]->half_flag = d->half_flag[f];

                }
                if (d->RTK_Npoint[f]) {
                    d->RTK_Npoint[f]->last_update_time = data->ros_time;
                    assert(d->half_flag[f] == d->RTK_Npoint[f]->half_flag);
                }
            }

            if (d->SPP_L[f] != 0) {
                if ((!d->SPP_Npoint[f]) || (d->SPP_Npoint[f]->LLI != d->SPP_LLI[f]) || condition3 || condition4) {
                    PBtype n;
                    n.sys = d->sys;
                    n.f = f;
                    n.value = 0;
                    n.continue_count = 0;
                    spp_phase_bias_variables[d->sat * 2 + f].push_back(n);
                    auto it2 = spp_phase_bias_variables[d->sat * 2 + f].end();
                    it2--;
                    d->SPP_Npoint[f] = &(*it2);
                    d->SPP_Npoint[f]->LLI = d->SPP_LLI[f];
                    d->SPP_Npoint[f]->half_flag = d->half_flag[f];

                }

                if (d->SPP_Npoint[f]) {
                    d->SPP_Npoint[f]->last_update_time = data->ros_time;
                }
            }

            if (d->SPP_P0[f] != 0) {
                if (!d->SPP_Npoint_PCottections[f]) {
                    PBtype n;
                    n.sys = d->sys;
                    n.f = f;
                    n.value = 0;
                    n.continue_count = 0;
                    pseudorange_correction_variables[d->sat * 2 + f].push_back(n);
                    auto it2 = pseudorange_correction_variables[d->sat * 2 + f].end();
                    it2--;
                    d->SPP_Npoint_PCottections[f] = &(*it2);

                }

                if (d->SPP_Npoint_PCottections[f]) {
                    d->SPP_Npoint_PCottections[f]->last_update_time = data->ros_time;
                }
            }

            if (d->RTK_Npoint[f])d->RTK_Npoint[f]->continue_count++;
            if (d->SPP_Npoint[f])d->SPP_Npoint[f]->continue_count++;
            if (d->SPP_Npoint_PCottections[f])d->SPP_Npoint_PCottections[f]->continue_count++;

        }

    }




    rover->marg_info_gnss = new MarginalizationInfo();

    std::set<double*>RemainPoint{para_pose[g2f[rover_count - 1]], para_speed_bias[g2f[rover_count - 1]], &blackvalue}; //




    for (uint8_t i = 0; i < rover->obs_count; i++) {
        ObsMea* d = rover->obs_data + i;
        for (uint8_t f = 0; f < NFREQ; f++) {
            if (d->RTK_Npoint[f]) {
                RemainPoint.insert(&(d->RTK_Npoint[f]->value));
            }
            if (d->SPP_Npoint[f]) {
                RemainPoint.insert(&(d->SPP_Npoint[f]->value));
            }
            if (d->SPP_Npoint_PCottections[f]) {
                RemainPoint.insert(&(d->SPP_Npoint_PCottections[f]->value));
            }
        }
    }
    ceres::Problem problem;
    ceres::Solver::Options options;

    PhaseBiasSaveAndReset();
    AddGnssResidual(MargeExcludeMode, RemainPoint, rover->marg_info_gnss, problem, options, rover_count - 1);
    rover->marg_info_gnss->marginalize(true, true);
    rover->marg_info_gnss->getParameterBlocks();
    PhaseBiasRestore();

    if (USE_SPP_PHASE || USE_RTK) {
        TicToc tic;
        //we use the dogleg strategy for global MyOptimization. For stability, we need to initialize the phase biases near the true value.
        //If the levengerg strategy is used for global MyOptimization, this block can be disabled.
        // MarginalizationFactor *factor = new MarginalizationFactor(rover->marg_info_gnss);
        // problem.AddResidualBlock(factor, 0, rover->marg_info_gnss->keep_block_addr);
        AddGnssResidual(NormalMode, std::set<double*> {}, 0, problem, options, rover_count - 1);
        // bool needoptimize=false;
        for (uint8_t i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            for (uint8_t f = 0; f < NFREQ; f++) {

                if (d->RTK_Npoint[f] && problem.HasParameterBlock(&(d->RTK_Npoint[f]->value))) {
                    if (d->RTK_Npoint[f]->continue_count > 10)problem.SetParameterBlockConstant(&(d->RTK_Npoint[f]->value));
                    // else needoptimize=true;
                }

                if (d->SPP_Npoint[f] && problem.HasParameterBlock(&(d->SPP_Npoint[f]->value))) {
                    if (d->SPP_Npoint[f]->continue_count > 10)problem.SetParameterBlockConstant(&(d->SPP_Npoint[f]->value));
                    // else needoptimize=true;
                }

                if (d->SPP_Npoint_PCottections[f] && problem.HasParameterBlock(&(d->SPP_Npoint_PCottections[f]->value))) {
                    if (d->SPP_Npoint_PCottections[f]->continue_count > 10)problem.SetParameterBlockConstant(&(d->SPP_Npoint_PCottections[f]->value));
                    // else needoptimize=true;
                }

            }
        }
        if (1) {
            problem.SetParameterBlockConstant(para_pose[rover_count + image_count - 1]);
            problem.SetParameterBlockConstant(para_speed_bias[rover_count + image_count - 1]);
            options.linear_solver_type = ceres::DENSE_SCHUR;
            options.initial_trust_region_radius = options.max_trust_region_radius = 1e15;
            options.max_num_iterations = 2;
            options.num_threads = 1;
            ceres::Solver::Summary summary;
            ceres::Solve(options, &problem, &summary);
            std::cout << summary.BriefReport() << "," << tic.toc() << std::endl;
        }

    }



    if (USE_GLOBAL_OPTIMIZATION) {
        MarginalizationFactor* factor = new MarginalizationFactor(rover->marg_info_gnss);
        rover->residualBlockId = my_problem.AddResidualBlock(factor, 0, rover->marg_info_gnss->keep_block_addr);
    }

    if (not_fix_count > Phase_ALL_RESET_COUNT)not_fix_count = 0;


}

void SWFOptimization::PhaseBiasSaveAndReset() {
    phase_bias_save.clear();
    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
            phase_bias_save.push_back(it1->value);
            it1->value = 0;
        }
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++) {
            phase_bias_save.push_back(it1->value);
            it1->value = 0;
        }
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++) {
            phase_bias_save.push_back(it1->value);
            it1->value = 0;
        }
    }
    phase_bias_save.push_back(blackvalue);
    blackvalue = 0;
}

void SWFOptimization::PhaseBiasRestore() {
    int index = 0;
    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
            it1->value = phase_bias_save[index++];
        }
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++) {
            it1->value = phase_bias_save[index++];
        }
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++) {
            it1->value = phase_bias_save[index++];
        }
    }
    blackvalue = phase_bias_save[index++];
    assert(index == (int)phase_bias_save.size());
    phase_bias_save.clear();
}



void SWFOptimization::MyOrdering(ceres::Problem& problem, ceres::Solver::Options& options) {
    if (!options.linear_solver_ordering) {
        options.linear_solver_ordering.reset(new ceres::ParameterBlockOrdering());
    }

    TicToc t;
    ceres::ParameterBlockOrdering* ordering = options.linear_solver_ordering.get();
    ordering->Clear();
    int ors = 1;
    double* pointer;
    std::set<double*>mark_add;
    mark_add.insert(last_marg_info->keep_block_addr.begin(), last_marg_info->keep_block_addr.end());
    mark_add.insert(ceres::internal::parameter_head.begin(), ceres::internal::parameter_head.end());

#define CONDITION problem.HasParameterBlock(pointer)&&!problem.IsParameterBlockConstant(pointer)&&mark_add.find(pointer)==mark_add.end()

    {
        assert(problem.HasParameterBlock(&blackvalue2));
        problem.SetParameterBlockVariable(&blackvalue2);
        ordering->AddElementToGroup(&blackvalue2, 0);
        mark_add.insert(&blackvalue2);
        std::vector<ceres::internal::ResidualBlock*>residual_blocks;
        problem.GetResidualBlocksForParameterBlock(&blackvalue2, &residual_blocks);
        assert(residual_blocks.size() == 1);
        residual_blocks[0]->is_use = true;
    }

    //feature points

    for (auto& it_per_id : f_manager.feature) {
#if USE_INVERSE_DEPTH
        pointer = &it_per_id.idepth_;
        if (CONDITION) {
            ordering->AddElementToGroup(pointer, 0);
            mark_add.insert(pointer);
        }
#else
        pointer = it_per_id.ptsInWorld.data();
        if (CONDITION) {
            ordering->AddElementToGroup(pointer, 0);
            mark_add.insert(pointer);
        }
#endif
    }

    //vel，acc bias，gyr bias
    int index = 0;
    for (int i = 0; i < rover_count + image_count; i++) {
        pointer = para_speed_bias[i];
        if (CONDITION) {
            if (index++ % 2 == 0) {
                ordering->AddElementToGroup(pointer, 0);
                mark_add.insert(pointer);
            }
        }
    }

    for (int i = 0; i < rover_count + image_count; i++) {
        pointer = para_speed_bias[i];
        if (CONDITION) {
            ordering->AddElementToGroup(pointer, ors++);
            mark_add.insert(pointer);
        }
    }

    //pos orientation
    for (int i = 0; i < rover_count + image_count; i++) {
        double* pointer = para_pose[i];
        if (CONDITION) {
            ordering->AddElementToGroup(pointer, ors++);
            mark_add.insert(pointer);
        }
    }

    //mag bias
    pointer = para_bmg;
    if (CONDITION) {
        ordering->AddElementToGroup(pointer, ors++);
        mark_add.insert(pointer);
    }


    //camera-IMU ex calibration
    for (int i = 0; i < NUM_OF_CAM; i++) {
        pointer = para_ex_Pose[i];
        if (CONDITION) {
            ordering->AddElementToGroup(pointer, ors++);
            mark_add.insert(pointer);
        }
    }

    pointer = &blackvalue;
    if (CONDITION) {
        ordering->AddElementToGroup(pointer, ors++);
        mark_add.insert(pointer);
    }

    //spp correction
    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++) {
            pointer = &it1->value;
            if (CONDITION) {
                ordering->AddElementToGroup(pointer, ors++);
                mark_add.insert(pointer);
            }
        }
    }
    //spp phase bias
    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++) {
            pointer = &it1->value;
            if (CONDITION) {
                ordering->AddElementToGroup(pointer, ors++);
                mark_add.insert(pointer);
            }
        }
    }

    //rb-sd phase bias
    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
            pointer = &it1->value;
            if (CONDITION) {
                ordering->AddElementToGroup(pointer, ors++);
                mark_add.insert(pointer);
            }
        }
    }



    mark_add.clear();
    mark_add.insert(ceres::internal::parameter_head.begin(), ceres::internal::parameter_head.end());
    //prior
    for (int i = 0; i < (int)last_marg_info->keep_block_addr.size(); i++) {
        pointer = last_marg_info->keep_block_addr[i];
        if (CONDITION) {
            ordering->AddElementToGroup(pointer, ors++);
        }
    }


    //remained states
    for (int i = 0; i < (int)ceres::internal::parameter_head.size(); i++) {
        double* pointer = ceres::internal::parameter_head[i];
        if (problem.HasParameterBlock(pointer) && !problem.IsParameterBlockConstant(pointer)) {
            ordering->AddElementToGroup(pointer, ors++);

        }
        assert(problem.HasParameterBlock(pointer));
    }



}
