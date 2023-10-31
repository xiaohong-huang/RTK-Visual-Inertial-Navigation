

#include "swf.h"
#include "../utility/visualization.h"
#include <thread>
#include <queue>
#include "../factor/pose0_factor.h"
#include "../factor/pose_local_parameterization.h"
#include "../factor/projection_factor.h"



SWFOptimization::SWFOptimization(): f_manager{Rs} {
    printf("init begins");
    ClearState();
    prev_time = -1;
    prev_time2 = -1;

    cur_time = 0;
    open_ex_estimation = 0;
    rover_count = 0;
    fix = false;
    imu_initialize = false;
    //initializing the global solver.
    my_options.linear_solver_type = ceres::DENSE_SCHUR;
    my_options.max_num_iterations = MAX_NUM_ITERATIONS;
    my_options.jacobi_scaling = 0;
    my_options.trust_region_strategy_type = ceres::DOGLEG;
    my_options.num_threads = 4;
    my_options.linear_solver_ordering.reset(new ceres::ParameterBlockOrdering());
    mag_mean.setZero();


}


void SWFOptimization::SetParameter() {
    for (int i = 0; i < NUM_OF_CAM; i++) {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
        cout << " exitrinsic cam " << i << endl  << ric[i] << endl << tic[i].transpose() << endl;
    }
    f_manager.setRic(ric);

    cout << "set G " << G.transpose() << endl;
    feature_tracker.readIntrinsicParameter(CAM_NAMES);
    projection_factor::sqrt_info = FOCAL_LENGTH / FEATUREWEIGHTINVERSE * Matrix2d::Identity();
    ProjectionTwoFrameOneCamFactor::sqrt_info = FOCAL_LENGTH / FEATUREWEIGHTINVERSE * Matrix2d::Identity();
    ProjectionTwoFrameTwoCamFactor::sqrt_info = FOCAL_LENGTH / FEATUREWEIGHTINVERSE * Matrix2d::Identity();
    ProjectionOneFrameTwoCamFactor::sqrt_info = FOCAL_LENGTH / FEATUREWEIGHTINVERSE * Matrix2d::Identity();

}

//need to fix for reseting the system.
void SWFOptimization::ClearState() {
    for (int i = 0; i < FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1; i++) {

        if (para_pose[i])delete para_pose[i];
        if (para_speed_bias[i])delete para_speed_bias[i];
        if (para_gnss_dt[i])delete para_gnss_dt[i];
        para_gnss_dt[i] = new double[13];
        para_pose[i] = new double[SIZE_POSE];
        para_speed_bias[i] = new double[SIZE_SPEEDBIAS];
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();

        for (int s = 0; s < 13; s++)
            para_gnss_dt[i][s] = 0;


        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr) {
            delete pre_integrations[i];
        }
        pre_integrations[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++) {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }

    first_imu = true,
    solver_flag = Initial;

    if (last_marg_info != nullptr)
        delete last_marg_info;

    last_marg_info = nullptr;

    f_manager.ClearState();

    rover_count = 0;
    image_count = 0;
    fix = false;

    acc_mean.setZero();
    gyr_mean.setZero();
    mag_mean.setZero();
    gyr_count = 0;
    Rwgw.setIdentity();


    para_bmg[0] = para_bmg[1] = para_bmg[2] = 0;
    base_pos.setZero();
    init_gnss = false;
    gnss_fix_solution_count = 0;
    pub_init = false;
    rtk_fix = false;
    InitRwgw.setIdentity();
    mag_yaw = 0;
}

//getting the gnss index or visual index according to the frame index, or reverse.
int SWFOptimization::ImageRoverId2FrameId(int image_rover_index, int mode) {
    int image_rover_index_tmp = image_rover_index;
    image_rover_index++;
    for (int i = 0; i < image_count + rover_count; i++) {
        if (frame_types[i] == mode) {
            image_rover_index--;
            if (image_rover_index == 0) {
                if (rover_count == 0 || image_count == 0) {
                    assert(image_rover_index_tmp == i);
                }
                return i;
            }
        }
    }
    assert(0);
    return -1;
}

//getting the pointer of the states.
void SWFOptimization::Vector2Double() {
    for (int i = 0; i < rover_count + image_count; i++) {
        para_pose[i][0] = Ps[i].x();
        para_pose[i][1] = Ps[i].y();
        para_pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_pose[i][3] = q.x();
        para_pose[i][4] = q.y();
        para_pose[i][5] = q.z();
        para_pose[i][6] = q.w();


        para_speed_bias[i][0] = Vs[i].x();
        para_speed_bias[i][1] = Vs[i].y();
        para_speed_bias[i][2] = Vs[i].z();

        para_speed_bias[i][3] = Bas[i].x();
        para_speed_bias[i][4] = Bas[i].y();
        para_speed_bias[i][5] = Bas[i].z();

        para_speed_bias[i][6] = Bgs[i].x();
        para_speed_bias[i][7] = Bgs[i].y();
        para_speed_bias[i][8] = Bgs[i].z();

    }

    for (int i = 0; i < NUM_OF_CAM; i++) {
        para_ex_Pose[i][0] = tic[i].x();
        para_ex_Pose[i][1] = tic[i].y();
        para_ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_ex_Pose[i][3] = q.x();
        para_ex_Pose[i][4] = q.y();
        para_ex_Pose[i][5] = q.z();
        para_ex_Pose[i][6] = q.w();
    }

    for (auto& it_per_id : f_manager.feature) {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < FEATURE_CONTINUE)continue;
    }



}


//saving the states from pointer.
void SWFOptimization::Double2Vector() {

    for (int i = 0; i < rover_count + image_count; i++) {

        Rs[i] = Quaterniond(para_pose[i][6], para_pose[i][3], para_pose[i][4], para_pose[i][5]).normalized().toRotationMatrix();

        Ps[i] = Vector3d(para_pose[i][0], para_pose[i][1], para_pose[i][2] ) ;

        Vs[i] = Vector3d(para_speed_bias[i][0], para_speed_bias[i][1], para_speed_bias[i][2]);

        Bas[i] = Vector3d(para_speed_bias[i][3], para_speed_bias[i][4], para_speed_bias[i][5]);

        Bgs[i] = Vector3d(para_speed_bias[i][6], para_speed_bias[i][7], para_speed_bias[i][8]);

    }



    for (int i = 0; i < NUM_OF_CAM; i++) {
        tic[i] = Vector3d(para_ex_Pose[i][0], para_ex_Pose[i][1], para_ex_Pose[i][2]);
        ric[i] = Quaterniond(para_ex_Pose[i][6],
                             para_ex_Pose[i][3],
                             para_ex_Pose[i][4],
                             para_ex_Pose[i][5]).normalized().toRotationMatrix();
    }

    for (auto& it_per_id : f_manager.feature) {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < FEATURE_CONTINUE)continue;
#if USE_INVERSE_DEPTH
        it_per_id.ptsInWorld = Rs[i2f[it_per_id.start_frame]] * (
                                   ric[0] * (it_per_id.feature_per_frame[0].point / it_per_id.idepth_) + tic[0] - Pbg
                               ) + Ps[i2f[it_per_id.start_frame]];
        if (it_per_id.idepth_ < 0) it_per_id.solve_flag = 2;
        else  it_per_id.solve_flag = 1;
#else
        Vector3d pts_cj = ric[0].transpose() * ( Rs[i2f[it_per_id.start_frame]].transpose() * (it_per_id.ptsInWorld - Ps[i2f[it_per_id.start_frame]]) + Pbg - tic[0]);
        it_per_id.idepth_ = 1.0 / pts_cj.z();
        if (it_per_id.idepth_ < 0) it_per_id.solve_flag = 2;
        else  it_per_id.solve_flag = 1;
#endif
    }



}




//
void SWFOptimization::SlideWindowFrame(int frameindex, int windowsize, bool updateIMU) {

    if (frameindex != 0) {
        for (unsigned int i = 0; i < dt_buf[frameindex + 1].size(); i++) { 
            pre_integrations[frameindex]->push_back(dt_buf[frameindex + 1][i], linear_acceleration_buf[frameindex + 1][i], angular_velocity_buf[frameindex + 1][i]);
            dt_buf[frameindex].push_back(dt_buf[frameindex + 1][i]);
            linear_acceleration_buf[frameindex].push_back(linear_acceleration_buf[frameindex + 1][i]);
            angular_velocity_buf[frameindex].push_back(angular_velocity_buf[frameindex + 1][i]);
        }
        std::swap(pre_integrations[frameindex], pre_integrations[frameindex + 1]);
        dt_buf[frameindex].swap(dt_buf[frameindex + 1]);
        linear_acceleration_buf[frameindex].swap(linear_acceleration_buf[frameindex + 1]);
        angular_velocity_buf[frameindex].swap(angular_velocity_buf[frameindex + 1]);
        if (frameindex >= 1 && updateIMU) {
            if (USE_GLOBAL_OPTIMIZATION) {
                IMUFactor* factor = new IMUFactor(pre_integrations[frameindex + 1]);
                my_problem.AddResidualBlock(factor, 0, para_pose[frameindex - 1],
                                            para_speed_bias[frameindex - 1],
                                            para_pose[frameindex + 1],
                                            para_speed_bias[frameindex + 1]);
            }

        }
    }



    for (int i = frameindex; i < windowsize - 1; i++) {
        headers[i] = headers[i + 1];
        frame_types[i] = frame_types[i + 1];
        Rs[i] = Rs[i + 1];
        mags[i] = mags[i + 1];
        Ps[i] = Ps[i + 1];
        Vs[i] = Vs[i + 1];
        Bas[i] = Bas[i + 1];
        Bgs[i] = Bgs[i + 1];

        std::swap(para_pose[i], para_pose[i + 1]);
        std::swap(para_speed_bias[i], para_speed_bias[i + 1]);
        std::swap(pre_integrations[i], pre_integrations[i + 1]);

        dt_buf[i].swap(dt_buf[i + 1]);
        linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
        angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);
    }
    delete para_pose[windowsize - 1];
    delete para_speed_bias[windowsize - 1];
    delete pre_integrations[windowsize - 1];

    para_pose[windowsize - 1] = new double[SIZE_POSE];
    para_speed_bias[windowsize - 1] = new double[SIZE_SPEEDBIAS];
    pre_integrations[windowsize - 1] = 0;

    dt_buf[windowsize - 1].clear();
    linear_acceleration_buf[windowsize - 1].clear();
    angular_velocity_buf[windowsize - 1].clear();
}


//marginalizing the select frames.
//param margeindex is the set of frame indexes that are selected to be marginalized.
void SWFOptimization::MargFrames(std::set<int> margeindex) {
    Vector2Double();

    std::set <double*>MargePoints = FindMargSet( margeindex);
    std::set <double*>MargePoints2 = MargePoints;


    if (USE_GLOBAL_OPTIMIZATION) {
        if (marg_flag == MargImagOld) {
            for (auto& it_per_id : f_manager.feature) {
                if (margeindex.find(i2f[it_per_id.start_frame]) != margeindex.end()) {

#if USE_INVERSE_DEPTH
                    if (my_problem.HasParameterBlock(&it_per_id.idepth_)) {
                        MargePoints2.insert(&it_per_id.idepth_);
                        assert(it_per_id.used_num >= FEATURE_CONTINUE);
                        assert(it_per_id.start_frame == 0);
                    }
#else
                    if (my_problem.HasParameterBlock(it_per_id.ptsInWorld.data())) {
                        MargePoints2.insert(it_per_id.ptsInWorld.data());
                        assert(it_per_id.used_num >= FEATURE_CONTINUE);
                        assert(it_per_id.start_frame == 0);
                    }

#endif
                }
            }
            GlobalMarge(MargePoints2);
        } else if (marg_flag == MargImagSecondNew) {
            MarginalizationInfo* marginalization_info = new MarginalizationInfo();
            ceres::Problem problem;
            ceres::Solver::Options options;
            assert(margeindex.size() == 1);
            AddAllResidual(MargeIncludeMode, MargePoints2, marginalization_info, problem, options, false, false, false);
        } else if (marg_flag == MargRoverOld) {
            MarginalizationInfo* marginalization_info = new MarginalizationInfo();
            ceres::Problem problem;
            ceres::Solver::Options options;
            assert(margeindex.size() == 1);
            AddAllResidual(MargeIncludeMode, MargePoints2, marginalization_info, problem, options, false, true, true);
        }
    } else {
        if (marg_flag == MargImagOld) {
            for (auto& it_per_id : f_manager.feature) {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                if (it_per_id.used_num < FEATURE_CONTINUE || margeindex.find(i2f[it_per_id.start_frame]) == margeindex.end())
                    continue;
                assert(it_per_id.start_frame == 0);
#if USE_INVERSE_DEPTH
                MargePoints2.insert(&it_per_id.idepth_);
#else
                MargePoints2.insert(it_per_id.ptsInWorld.data());
#endif
                if (marg_flag == MargImagSecondNew)
                    assert(0);
            }
        }

        MarginalizationInfo* marginalization_info = new MarginalizationInfo();
        ceres::Problem problem;
        ceres::Solver::Options options;
        if (marg_flag == MargImagOld)
            AddAllResidual(MargeIncludeMode2, MargePoints2, marginalization_info, problem, options, true, true, true);
        else if (marg_flag == MargImagSecondNew)
            AddAllResidual(MargeIncludeMode, MargePoints2, marginalization_info, problem, options, false, false, false);
        else if (marg_flag == MargRoverOld)
            AddAllResidual(MargeIncludeMode, MargePoints2, marginalization_info, problem, options, false, true, true);
        else assert(0);
    }


    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it = rtk_phase_bias_variables[i].begin(), it_next = rtk_phase_bias_variables[i].begin(); it != rtk_phase_bias_variables[i].end(); it = it_next) {
            it_next++;
            if (!it->use) {
                rtk_phase_bias_variables[i].erase(it);
            }
        }
        for (auto it = spp_phase_bias_variables[i].begin(), it_next = spp_phase_bias_variables[i].begin(); it != spp_phase_bias_variables[i].end(); it = it_next) {
            it_next++;
            if (!it->use) {
                spp_phase_bias_variables[i].erase(it);
            }
        }

    }

    if (USE_GLOBAL_OPTIMIZATION) {
        for (auto it = MargePoints.begin(); it != MargePoints.end(); it++) {
            if (my_problem.HasParameterBlock(*it))my_problem.RemoveParameterBlock(*it);
        }
    }
#if USE_INVERSE_DEPTH
    for (auto& it_per_id : f_manager.feature) {
        it_per_id.used_num = it_per_id.feature_per_frame.size();

        if ( it_per_id.used_num >= FEATURE_CONTINUE  && (
                    (margeindex.find(i2f[it_per_id.start_frame]) != margeindex.end()) ||
                    (margeindex.find(i2f[it_per_id.endFrame()]) != margeindex.end() && it_per_id.used_num == 2))  ) {
            if (my_problem.HasParameterBlock(&it_per_id.idepth_)) {
                my_problem.RemoveParameterBlock(&it_per_id.idepth_);
            }
        }
    }
#endif

    for (auto& it_per_id : f_manager.feature) {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < FEATURE_CONTINUE)continue;
#if USE_INVERSE_DEPTH
        if (my_problem.HasParameterBlock(&it_per_id.idepth_)) {
            std::vector<ceres::internal::ResidualBlock*>residual_blocks;
            my_problem.GetResidualBlocksForParameterBlock(&it_per_id.idepth_, &residual_blocks);
            if (residual_blocks.size() == 0) {
                my_problem.RemoveParameterBlock(&it_per_id.idepth_);
                assert(0);
            }
        }
#else
        if (my_problem.HasParameterBlock(it_per_id.ptsInWorld.data())) {
            std::vector<ceres::internal::ResidualBlock*>residual_blocks;
            my_problem.GetResidualBlocksForParameterBlock(it_per_id.ptsInWorld.data(), &residual_blocks);
            if (residual_blocks.size() == 0) {
                my_problem.RemoveParameterBlock(it_per_id.ptsInWorld.data());
                assert(0);
            }
        }
#endif
    }

}
//getting the states (pointer) that will being marginzlized.
//param margeindex is the set of frame indexes that are selected to be marginalized.
std::set <double*> SWFOptimization::FindMargSet(std::set<int> margeindex) {
    std::set<double*>MargePoints;
    for (auto it = margeindex.begin(); it != margeindex.end(); it++) {
        MargePoints.insert(para_pose[*it]);
        MargePoints.insert(para_speed_bias[*it]);
    }


    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it = rtk_phase_bias_variables[i].begin(); it != rtk_phase_bias_variables[i].end(); it++) {
            it->use = false;
        }
        for (auto it = spp_phase_bias_variables[i].begin(); it != spp_phase_bias_variables[i].end(); it++) {
            it->use = false;
        }
        for (auto it = pseudorange_correction_variables[i].begin(); it != pseudorange_correction_variables[i].end(); it++) {
            it->use = false;
        }
    }

    for (int ir = 0; ir < rover_count; ir++) {
        if (margeindex.find(g2f[ir]) != margeindex.end())continue;
        mea_t* rover = rovers[ir];
        for (uint8_t i = 0; i < rover->obs_count; i++) {
            ObsMea* d = rover->obs_data + i;
            for (uint8_t f = 0; f < NFREQ; f++) {
                if (d->RTK_Npoint[f] )d->RTK_Npoint[f]->use = true;
                if (d->SPP_Npoint[f] )d->SPP_Npoint[f]->use = true;
                if (d->SPP_Npoint_PCottections[f] )d->SPP_Npoint_PCottections[f]->use = true;
            }
        }
    }

    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it1 = rtk_phase_bias_variables[i].begin(); it1 != rtk_phase_bias_variables[i].end(); it1++) {
            if (!it1->use) {
                MargePoints.insert(&it1->value);
                // std::cout<<"margepoint:"<<i<<std::endl;
            }
        }
        for (auto it1 = spp_phase_bias_variables[i].begin(); it1 != spp_phase_bias_variables[i].end(); it1++) {
            if (!it1->use) {
                MargePoints.insert(&it1->value);
                // std::cout<<"margepoint:"<<i<<std::endl;
            }
        }
        for (auto it1 = pseudorange_correction_variables[i].begin(); it1 != pseudorange_correction_variables[i].end(); it1++) {
            if (!it1->use) {
                MargePoints.insert(&it1->value);
                // std::cout<<"margepoint:"<<i<<std::endl;
            }
        }
    }
    return MargePoints;
}

//marginalizing the selected gnss frames.
MarginalizationInfo* SWFOptimization::MargGNSSFrames(std::set<int> margeindex, IMUGNSSBase* IMUGNSSmeasurement) {
    Vector2Double();
    std::set <double*>MargePoints = FindMargSet( margeindex);


    MarginalizationInfo* marginalization_info = new MarginalizationInfo();
    ceres::Problem problem;
    ceres::Solver::Options options;

    AddAllResidual(GNSSMargIncludeMode, MargePoints, marginalization_info, problem, options, false, true, true);


    for (uint8_t i = 0; i < MAXSATNUM * 2; i++) {
        for (auto it = rtk_phase_bias_variables[i].begin(), it_next = rtk_phase_bias_variables[i].begin(); it != rtk_phase_bias_variables[i].end(); it = it_next) {
            it_next++;
            if (!it->use) {
                rtk_phase_bias_variables[i].erase(it);
            }
        }
        for (auto it = spp_phase_bias_variables[i].begin(), it_next = spp_phase_bias_variables[i].begin(); it != spp_phase_bias_variables[i].end(); it = it_next) {
            it_next++;
            if (!it->use) {
                spp_phase_bias_variables[i].erase(it);
            }
        }
    }

    if (USE_GLOBAL_OPTIMIZATION) {
        for (auto it = MargePoints.begin(); it != MargePoints.end(); it++) {
            if (my_problem.HasParameterBlock(*it)) {
                assert(my_problem.ParameterBlockSize(*it) == 1);
                my_problem.RemoveParameterBlock(*it);
                IMUGNSSmeasurement->residualBlockId = 0;
            }
        }
    }

    return marginalization_info;


}



void SWFOptimization::SlideWindow() {
    TicToc t_marg;

    if (image_count > FEATURE_WINDOW_SIZE) {
        if (imag_marg_index != 0) {
            marg_flag = MargImagSecondNew;
        } else {
            marg_flag = MargImagOld;
        }
    } else if (frame_types[0] == GnssFrame && rover_count > 10) {
        return;
        marg_flag = MargRoverOld;
    } else {
        return;
    }

    std::set<int>margeindex;

    if (marg_flag == MargImagOld) {
        Eigen::Vector3d P0 = Ps[i2f[0]];
        Eigen::Matrix3d R0 = Rs[i2f[0]];
        Eigen::Vector3d P1 = Ps[i2f[1]];
        Eigen::Matrix3d R1 = Rs[i2f[1]];

        for (int i = 0; i < i2f[1]; i++) {
            margeindex.insert(i);
        }
        std::cout << "MargImagOld:" << image_count << "," << rover_count << std::endl;
        MargFrames(margeindex);
        int frame_counts = image_count + rover_count;
        for (int i = 0; i < i2f[1]; i++) {
            if (frame_types[0] == GnssFrame) {
                free(rovers[0]);
                for (int ri = 0; ri < rover_count; ri++) {
                    rovers[ri] = rovers[ri + 1];
                }
                rover_count--;
            }
            SlideWindowFrame( 0, frame_counts, USE_IMU);
            frame_counts--;
        }

        SlideWindowOld( P0, R0, P1, R1, tic[0], ric[0], Pbg);
        if (imu_gnss_factor[0]) {
            delete imu_gnss_factor[0];
            imu_gnss_factor[0] = 0;
        }
        for (int i = 2; i < image_count + 1; i++)std::swap(imu_gnss_factor[i - 1], imu_gnss_factor[i]);
        if (imu_gnss_factor[image_count]) {
            delete imu_gnss_factor[image_count];
            imu_gnss_factor[image_count] = 0;
        }
        image_count--;
        UpdateVisualGnssIndex();
        if (USE_GLOBAL_OPTIMIZATION && USE_INVERSE_DEPTH) {
            AddFeature2Problem();
        }
    } else if (marg_flag == MargImagSecondNew) {
        margeindex.insert(i2f[image_count - 2]);
        std::cout << "MargImagSecondNew:" << image_count << "," << rover_count << std::endl;
        MargFrames(margeindex);
        SlideWindowFrame(i2f[image_count - 2], image_count + rover_count, USE_IMU);
        SlideWindowNew();

        if (imu_gnss_factor[image_count - 1] || imu_gnss_factor[image_count - 2]) {
            if (!imu_gnss_factor[image_count - 2])imu_gnss_factor[image_count - 2] = new IMUGNSSBase(
                    para_pose[i2f[image_count - 2] - 1], para_speed_bias[i2f[image_count - 2] - 1], &my_problem);
            IMUGNSSBase* IMUGNSSmeasurement = imu_gnss_factor[image_count - 2];
            int index1 = i2f[image_count - 2] - 1, index2 = i2f[image_count - 1] - 1; //21,10
            for (int j = index1 + 1; j < index2; j++) { //11
                IMUGNSSmeasurement->AddMargInfo(rovers[f2g[j + 1]]->marg_info_gnss, pre_integrations[j],
                                                para_pose[j], para_speed_bias[j]);
            }
            IMUGNSSmeasurement->SetLastImuFactor(pre_integrations[index2], para_pose[index2], para_speed_bias[index2]);
            std::swap(imu_gnss_factor[image_count - 1], imu_gnss_factor[image_count - 2]);
        }

        for (int i = image_count - 1; i < image_count + 1; i++)std::swap(imu_gnss_factor[i - 1], imu_gnss_factor[i]);
        if (imu_gnss_factor[image_count]) {
            delete imu_gnss_factor[image_count];
            imu_gnss_factor[image_count] = 0;
        }

        image_count--;
        UpdateVisualGnssIndex();
        if (USE_GLOBAL_OPTIMIZATION && USE_INVERSE_DEPTH) {
            AddFeature2Problem();
        }

    } else if (marg_flag == MargRoverOld) {
        margeindex.insert(g2f[0]);
        std::cout << "MargRoverOld:" << image_count << "," << rover_count << std::endl;
        MargFrames(margeindex);
        SlideWindowFrame(g2f[0], image_count + rover_count, USE_IMU);
        free(rovers[0]);
        for (int ri = 0; ri < rover_count; ri++)rovers[ri] = rovers[ri + 1];
        rover_count--;
        UpdateVisualGnssIndex();
        if (USE_GLOBAL_OPTIMIZATION) {
            if (frame_types[0] == GnssFrame) {
                MarginalizationFactor* factor = new MarginalizationFactor(rovers[0]->marg_info_gnss);
                rovers[0]->residualBlockId = my_problem.AddResidualBlock(factor, 0, rovers[0]->marg_info_gnss->keep_block_addr);
            }
            if (USE_IMU) {
                assert(pre_integrations[1]);
                IMUFactor* factor = new IMUFactor(pre_integrations[1]);
                my_problem.AddResidualBlock(factor, 0, para_pose[0], para_speed_bias[0], para_pose[1], para_speed_bias[1]);
            }
        }
        ResetImuGnssFactor(0, 0);

    }
    std::cout << "marge time:" << t_marg.toc() << std::endl;



}



//updating the frame index, gnss frame index, and visual frame index.
void SWFOptimization::UpdateVisualGnssIndex() {
    for (int k = 0; k < image_count; k++) {
        i2f[k] = ImageRoverId2FrameId(k, ImagFrame);
        f2i[i2f[k]] = k;
    }
    for (int k = 0; k < rover_count; k++) {
        g2f[k] = ImageRoverId2FrameId(k, GnssFrame);
        f2g[g2f[k]] = k;
    }
}


//regenerating the gnss-imu factor.
void SWFOptimization::ResetImuGnssFactor(int IMUGNSSindex, MarginalizationInfo* gnss_middle_marginfo) {

    if (imu_gnss_factor[IMUGNSSindex]) {
        if (USE_GLOBAL_OPTIMIZATION)
            if (imu_gnss_factor[IMUGNSSindex]->residualBlockId)
                my_problem.RemoveResidualBlock(imu_gnss_factor[IMUGNSSindex]->residualBlockId);
        delete imu_gnss_factor[IMUGNSSindex];
    }
    imu_gnss_factor[IMUGNSSindex] = 0;
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

    if (index2 - index1 > 1) {
        imu_gnss_factor[IMUGNSSindex] = new IMUGNSSBase(para_pose[index1], para_speed_bias[index1], &my_problem);
        IMUGNSSBase* IMUGNSSmeasurement = imu_gnss_factor[IMUGNSSindex];
        for (int j = index1 + 1; j < index2; j++) {
            IMUGNSSmeasurement->AddMargInfo(rovers[f2g[j]]->marg_info_gnss, pre_integrations[j],
                                            para_pose[j], para_speed_bias[j]);
        }
        if (gnss_middle_marginfo)IMUGNSSmeasurement->AddMidMargInfo( gnss_middle_marginfo);
        IMUGNSSmeasurement->SetLastImuFactor(pre_integrations[index2], para_pose[index2], para_speed_bias[index2]);

    }
}




//adding the new gnss, imu information to the corresponding gnss-imu factor.
void SWFOptimization::UpdateImuGnssFactor() {



    if (image_count + rover_count >= 3) {
        int j = image_count + rover_count - 2;
        if (frame_types[j] == GnssFrame) {
            int index = frame_types[j + 1] == ImagFrame ? (image_count - 1) : image_count;
            if (!imu_gnss_factor[index])imu_gnss_factor[index] = new IMUGNSSBase(
                    para_pose[image_count + rover_count - 3], para_speed_bias[image_count + rover_count - 3], &my_problem);
            IMUGNSSBase* IMUGNSSmeasurement = imu_gnss_factor[index];
            IMUGNSSmeasurement->AddMargInfo(rovers[f2g[j]]->marg_info_gnss, pre_integrations[j],
                                            para_pose[j], para_speed_bias[j]);
            IMUGNSSmeasurement->SetLastImuFactor(pre_integrations[j + 1], para_pose[j + 1], para_speed_bias[j + 1]);
        }
    }

}
//finding the avaliable gnss observations.
uint8_t getVariableUseSingleNum(mea_t* obs_data) {
    uint8_t i, j = 0;
    ObsMea* datai = obs_data->obs_data;
    for (i = 0; i < obs_data->obs_count; i++) {
        ObsMea* d = datai + i;
        if (d->SVH == 0 && d->SPP_P[0] != 0 && d->SPP_Pstd[0] < 2) {
            j += 1;
        }
    }
    return j;
}

//main process
void SWFOptimization::MeasurementProcess() {
#define condition1 (USE_IMAGE&&USE_GNSS&&USE_IMU&&!feature_buf.empty()&&!rover_buf.empty())
#define condition2 (!USE_IMAGE&&USE_GNSS&&!rover_buf.empty())
#define condition3 (USE_IMAGE&&!USE_GNSS&&!feature_buf.empty())


    while (condition1 || condition2 || condition3) {
        std::cout << rover_buf.size() << "," << feature_buf.size() << std::endl;
        TicToc t_process;
        pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1> > > > > feature;
        if (condition1 || condition3)
            feature = feature_buf.front();

        mea_t* rover = rover_buf.front();;
        FrameType modeflag = ErroFrame;

        if (condition1) {
            if (!rover_buf.empty() && rover && rover->ros_time < feature.first) {
                cur_time = rover->ros_time;
                modeflag = GnssFrame;
            } else {
                cur_time = feature.first;
                modeflag = ImagFrame;
            }
            if (!rover_buf.empty() && rover && modeflag == GnssFrame && getVariableUseSingleNum(rover) < 8) {
                rover_buf.pop_front();
                std::cout << "                                                                                                         pop GNSS\r\n";
                continue;
            }
            if (!rover_buf.empty() && rover && modeflag == ImagFrame && rover->ros_time - cur_time <= 0.005 && getVariableUseSingleNum(rover) >= 8) {
                feature_buf.pop();
                std::cout << "                                                                                                         pop IMAGE\r\n";
                continue;
            }
        }
        if (condition2) {
            solver_flag = NonLinear;
            cur_time = rover->ros_time;
            modeflag = GnssFrame;
            if (!rover_buf.empty() && rover && modeflag == GnssFrame && getVariableUseSingleNum(rover) < 8) {
                rover_buf.pop_front();
                std::cout << "                                                                                                         pop GNSS\r\n";
                continue;
            }
        }
        if (condition3) {
            cur_time = feature.first;
            modeflag = ImagFrame;
        }

        if (cur_time <= prev_time2 || cur_time - prev_time2 <= 0.005||(USE_IMU&&acc_buf.front().first>cur_time)) { //
            if (modeflag == GnssFrame) rover_buf.pop_front();
            else feature_buf.pop();
            continue;
        }

        assert(modeflag != ErroFrame);
        if (!ImuAvailable(cur_time) && USE_IMU)return;



        if (modeflag == GnssFrame) {
            rover_buf.pop_front();
            rover_count++;
            need_Nresolve = true;
            frame_types[rover_count + image_count - 1] = GnssFrame;

        } else {
            feature_buf.pop();
            image_count++;
            frame_types[rover_count + image_count - 1] = modeflag;
        }

        MagProcess(cur_time);
        ImuIntegrate();


        if (USE_GLOBAL_OPTIMIZATION)
            my_problem.AddParameterBlock(para_pose[image_count + rover_count - 1], SIZE_POSE, new PoseLocalParameterization());

        prev_time2 = cur_time;

        UpdateVisualGnssIndex();

        // if(modeflag==GnssFrame)
        assert(fabs(headers[image_count + rover_count - 1] - cur_time) < 1e-4);
        headers[image_count + rover_count - 1] = cur_time;
        // else headers[image_count+rover_count-1] = feature.first;

        assert(rover_count + image_count < FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE);

        if (init_gnss) {
            double pos[3];
            Matrix3d Rwwg;
            Eigen::Vector3d xyz = base_pos + Ps[rover_count / 2];
            ecef2pos(xyz.data(), pos);
            xyz2enu(pos, Rwwg.data());//transform vector from ecef to enu
            Rwgw = Rwwg.transpose();

        }
        if (modeflag == ImagFrame)ImagePreprocess(feature.second);
        else GnssProcess(rover);
        UpdateImuGnssFactor();

        if (solver_flag != Initial && (modeflag == ImagFrame || frame_types[rover_count + image_count - 2] == GnssFrame)) {
            MyOptimization();
            IntegerSolve();
        }

        std::cout << std::endl;
        for (int i = 0; i < image_count + 1; i++) {
            if (imu_gnss_factor[i])
                std::cout << imu_gnss_factor[i]->update_round * 1.0 / imu_gnss_factor[i]->evaluate_round << ",";
        }
        if (modeflag == ImagFrame)ImagePostprocess();


        SlideWindow();
        MiddleMargGnssFrame();
        if (modeflag == ImagFrame)
            f_manager.removeOut2(image_count, my_problem);

        Ps[rover_count + image_count] = Ps[rover_count + image_count - 1];
        Rs[rover_count + image_count] = Rs[rover_count + image_count - 1];
        Vs[rover_count + image_count] = Vs[rover_count + image_count - 1];
        Bas[rover_count + image_count] = Bas[rover_count + image_count - 1];
        Bgs[rover_count + image_count] = Bgs[rover_count + image_count - 1];
        headers[rover_count + image_count] = headers[rover_count + image_count - 1];


        if (solver_flag != Initial) {
            PubData();
        }

        if (modeflag == ImagFrame) {
            static double t_process2 = 0;
            static int t_count = 0;
            double ts = t_process.toc();
            t_process2 += ts;
            t_count += 1;
            printf("process measurement time: %f   ,%f   ,%d\n", ts, t_process2 / t_count, (int)modeflag);
        }

    }



}


//publicating and saving results.
void SWFOptimization::PubData() {


    for (int k = 0; k < image_count; k++)i2f[k] = ImageRoverId2FrameId(k, ImagFrame);
    for (int k = 0; k < rover_count; k++)g2f[k] = ImageRoverId2FrameId(k, GnssFrame);

    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(headers[image_count + rover_count - 1]);
    if (!pub_init) {
        resetpot(*this, header);
        pub_init = true;
    }
    printStatistics(*this, 0);
    pubOdometry(*this, header);
    pubCameraPose(*this, header);
    pubPointCloud(*this, header);
    pubKeyframe(*this);



}



