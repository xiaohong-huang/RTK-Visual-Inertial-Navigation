

#include "swf.h"
#include "../utility/visualization.h"
#include <thread>
#include <queue>



SWFOptimization::SWFOptimization() {
    printf("init begins");
    ClearState();
    prev_time = -1;
    cur_time = 0;
    R_WI_WC.setIdentity();
    InitRwgw.setIdentity();

}


void SWFOptimization::SetParameter() {

    LOG_OUT << "set G " << G.transpose() << endl;
    feature_tracker.readIntrinsicParameter(CAM_NAMES);

}



//need to fix for reseting the system.
void SWFOptimization::ClearState() {
    mag_mean.setZero();
    visual_inertial_bases_global.resize(SWF_SIZE_OUT + 1, 0);
    idepths_all.resize(SWF_SIZE_OUT + 1);
    Ps.resize(SWF_WINDOW_SIZE + 1);
    Vs.resize(SWF_WINDOW_SIZE + 1);
    Rs.resize(SWF_WINDOW_SIZE + 1);
    Bas.resize(SWF_WINDOW_SIZE + 1);
    Bgs.resize(SWF_WINDOW_SIZE + 1);
    headers.resize(SWF_WINDOW_SIZE + 1);
    para_pose.resize(SWF_WINDOW_SIZE + 1, 0);
    para_speed_bias.resize(SWF_WINDOW_SIZE + 1, 0);
    gnss_imu_infos.resize(SWF_WINDOW_SIZE + 1, 0);
    old_time_shift = 0;
    P_WI_WC.setZero();
    old_P.setZero();

    for (int i = 0; i < SWF_WINDOW_SIZE + 1; i++) {

        if (para_pose[i]) delete para_pose[i];
        if (para_speed_bias[i]) delete para_speed_bias[i];
        if (gnss_imu_infos[i] != nullptr) delete gnss_imu_infos[i];

        para_pose[i] = new double[SIZE_POSE];
        para_speed_bias[i] = new double[SIZE_SPEEDBIAS];

        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        gnss_imu_infos[i] = nullptr;
    }

    for (int i = 0; i < SWF_SIZE_OUT + 1; i++) {
        if (visual_inertial_bases_global[i] != nullptr)
            delete visual_inertial_bases_global[i];
        visual_inertial_bases_global[i] = nullptr;
    }


    solver_flag = Initial;
    last_marg_info = nullptr;
    image_count = 0;
    acc_count = 0;

    f_manager.ClearState();
    acc_mean.setZero();

    acc_scale = Eigen::Vector3d({1, 1, 1});
    gyr_scale = Eigen::Vector3d({1, 1, 1});

    if (last_marg_info != nullptr)delete last_marg_info;
}


//getting the pointer of the states.
void SWFOptimization::Vector2Double() {



    Quaterniond q {R_WI_WC};
    para_global_pos[0] = P_WI_WC.x();
    para_global_pos[1] = P_WI_WC.y();
    para_global_pos[2] = P_WI_WC.z();

    para_global_pos[3] = q.x();
    para_global_pos[4] = q.y();
    para_global_pos[5] = q.z();
    para_global_pos[6] = q.w();



    for (int i = 0; i < image_count; i++) {
        gnss_imu_infos[i]->Vector2Double();
        Eigen::Vector3d Pc = Ps[i] + Rs[i] * TIC[0] - P_WI_WC;
        Eigen::Matrix3d Rc = Rs[i] * RIC[0];
        Pc = R_WI_WC.transpose() * Pc / scale_factor;
        Rc = R_WI_WC.transpose() * Rc;





        para_pose[i][0] = Pc.x();
        para_pose[i][1] = Pc.y();
        para_pose[i][2] = Pc.z();
        Quaterniond q{Rc};
        para_pose[i][3] = q.x();
        para_pose[i][4] = q.y();
        para_pose[i][5] = q.z();
        para_pose[i][6] = q.w();
        Eigen::Vector3d Vg = Vs[i] + Rs[i] * (Utility::skewSymmetric(gnss_imu_infos[i]->pre_integrations[gnss_imu_infos[i]->pre_integrations.size() - 1]->gyrj - Bgs[i]) * Pbg);
        para_speed_bias[i][0] = Vg.x();
        para_speed_bias[i][1] = Vg.y();
        para_speed_bias[i][2] = Vg.z();

        para_speed_bias[i][3] = Bas[i].x();
        para_speed_bias[i][4] = Bas[i].y();
        para_speed_bias[i][5] = Bas[i].z();

        para_speed_bias[i][6] = Bgs[i].x();
        para_speed_bias[i][7] = Bgs[i].y();
        para_speed_bias[i][8] = Bgs[i].z();

    }

    if (ESTIMATE_EXTRINSIC) {
        para_extrinsic[0] = TIC[0].x();
        para_extrinsic[1] = TIC[0].y();
        para_extrinsic[2] = TIC[0].z();
        para_extrinsic[3] = QIC[0].x();
        para_extrinsic[4] = QIC[0].y();
        para_extrinsic[5] = QIC[0].z();
        para_extrinsic[6] = QIC[0].w();
    }
#if USE_STEREO
    {
        para_extrinsic2[0] = TC2C1.x();
        para_extrinsic2[1] = TC2C1.y();
        para_extrinsic2[2] = TC2C1.z();
        Eigen::Quaterniond q(RC2C1);
        para_extrinsic2[3] = q.x();
        para_extrinsic2[4] = q.y();
        para_extrinsic2[5] = q.z();
        para_extrinsic2[6] = q.w();
    }

#endif
}


//saving the states from pointer.
void SWFOptimization::Double2Vector() {

    R_WI_WC = Quaterniond(para_global_pos[6], para_global_pos[3], para_global_pos[4], para_global_pos[5]).normalized().toRotationMatrix();
    P_WI_WC = Vector3d(para_global_pos[0], para_global_pos[1], para_global_pos[2] ) ;

    if (ESTIMATE_EXTRINSIC) {
        TIC[0] = Vector3d(para_extrinsic[0], para_extrinsic[1], para_extrinsic[2] ) ;
        QIC[0] = Quaterniond(para_extrinsic[6], para_extrinsic[3], para_extrinsic[4], para_extrinsic[5]).normalized();
        RIC[0] = QIC[0].toRotationMatrix();
    }
#if USE_STEREO
    {
        TC2C1 = Vector3d(para_extrinsic2[0], para_extrinsic2[1], para_extrinsic2[2] ) ;
        RC2C1 = Quaterniond(para_extrinsic2[6], para_extrinsic2[3], para_extrinsic2[4], para_extrinsic2[5]).normalized().toRotationMatrix();

        RIC[1] = RIC[0] * RC2C1.transpose();
        TIC[1] = TIC[0] - RIC[1] * TC2C1;
    }
#endif
    for (int i = 0; i < image_count; i++) {
        gnss_imu_infos[i]->Double2Vector();

        Eigen::Vector3d Pc;
        Eigen::Matrix3d Rc;

        Rc = Quaterniond(para_pose[i][6], para_pose[i][3], para_pose[i][4], para_pose[i][5]).normalized().toRotationMatrix();
        Pc = Vector3d(para_pose[i][0], para_pose[i][1], para_pose[i][2] ) ;

        Pc = R_WI_WC * Pc * scale_factor + P_WI_WC;
        Rc = R_WI_WC * Rc;

        Rs[i] = Rc * RIC[0].transpose();
        Ps[i] = Pc - Rs[i] * TIC[0];
        Eigen::Vector3d Vg = Vector3d(para_speed_bias[i][0], para_speed_bias[i][1], para_speed_bias[i][2]);
        Vs[i] = Vg - Rs[i] * (Utility::skewSymmetric(gnss_imu_infos[i]->pre_integrations[gnss_imu_infos[i]->pre_integrations.size() - 1]->gyrj - Bgs[i]) * Pbg);
        Bas[i] = Vector3d(para_speed_bias[i][3], para_speed_bias[i][4], para_speed_bias[i][5]);
        Bgs[i] = Vector3d(para_speed_bias[i][6], para_speed_bias[i][7], para_speed_bias[i][8]);

    }

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) if (NUM_THREADS > 1)
    for (int j = 1; j < image_count; j++)
        gnss_imu_infos[j]->reprogate(Bas[j - 1], Bgs[j - 1], acc_scale, gyr_scale);

    for (auto& it_per_id : f_manager.feature) {
        if (!it_per_id.valid)continue;
        it_per_id.solve_flag = 1;
        for (int i = 0; i < SWF_SIZE_OUT; i++) {
            if (idepths_all[i].find(it_per_id.feature_id) != idepths_all[i].end() && idepths_all[i][it_per_id.feature_id] < 1 / 500.0)
                it_per_id.solve_flag = 2;
        }
    }

}


//
void SWFOptimization::SlideWindowFrame(int frameindex, int windowsize, bool updateIMU) {

    if (frameindex != 0) {
        gnss_imu_infos[frameindex]->merge_gnss_imu_factor(gnss_imu_infos[frameindex + 1]);

        std::swap(gnss_imu_infos[frameindex], gnss_imu_infos[frameindex + 1]);
    }

    for (int i = frameindex; i < windowsize - 1; i++) {
        headers[i] = headers[i + 1];
        Rs[i] = Rs[i + 1];
        Ps[i] = Ps[i + 1];
        Vs[i] = Vs[i + 1];
        Bas[i] = Bas[i + 1];
        Bgs[i] = Bgs[i + 1];


        std::swap(para_pose[i], para_pose[i + 1]);
        std::swap(para_speed_bias[i], para_speed_bias[i + 1]);
        std::swap(gnss_imu_infos[i], gnss_imu_infos[i + 1]);


    }
    delete para_pose[windowsize - 1];
    delete para_speed_bias[windowsize - 1];
    delete gnss_imu_infos[windowsize - 1];

    para_pose[windowsize - 1] = new double[SIZE_POSE];
    para_speed_bias[windowsize - 1] = new double[SIZE_SPEEDBIAS];
    gnss_imu_infos[windowsize - 1] = 0;

}


//marginalizing the select frames.
//param margeindex is the set of frame indexes that are selected to be marginalized.
void SWFOptimization::MargFrames() {

    if (marg_flag == MargImagOld) {
        double oldkminmu = kMinMu;
        if (boolsmallmuwhenmarge)
            kMinMu = 1e-6;

        OptimizationOrMarginalization( MargeMode);
        have_hist = 1;

        if (boolsmallmuwhenmarge)
            kMinMu = oldkminmu;
    }
}




void SWFOptimization::SlideWindow() {
    TicToc t_marg;
    if (solver_flag != Initial) {
        if (imag_marg_index != 0) marg_flag = MargImagSecondNew;
        else marg_flag = MargImagOld;
    } else
        return;
    assert(image_count >= 2);
    if (gnss_imu_infos[image_count - 2]->Rovers.size() > 5)marg_flag = MargImagOld;


    if (marg_flag == MargImagOld && image_count <= SWF_WINDOW_SIZE)
        return;

    if (marg_flag == MargImagOld) {
        Eigen::Vector3d P0 = Ps[0];
        Eigen::Matrix3d R0 = Rs[0];
        Eigen::Vector3d P1 = Ps[1];
        Eigen::Matrix3d R1 = Rs[1];

        MargFrames();

        if (visual_inertial_bases_global[0])delete visual_inertial_bases_global[0];
        visual_inertial_bases_global[0] = 0;
        for (int i = 0; i < (int)visual_inertial_bases_global.size() - 1; i++)
            visual_inertial_bases_global[i] = visual_inertial_bases_global[i + 1];

        for (int i = 0; i < SWF_SIZE_IN; i++)
            SaveKefPos(i);

        int frame_counts = image_count;
        for (int i = 0; i < SWF_SIZE_IN; i++) {
            SlideWindowFrame( 0, frame_counts, 1);
            frame_counts--;
        }
        SlideWindowOld( P0, R0, P1, R1, TIC[0], RIC[0]);
        image_count -= SWF_SIZE_IN;

#if USE_ASSERT
        for (int i = 0; i < (int)last_marg_info->keep_block_size.size(); i++) {
            if (last_marg_info->keep_block_size[i] == 1)
                ASSERT(fabs(last_marg_info->keep_block_addr[i][0] - last_marg_info->keep_block_data[i][0]) < 1e4);
        }
#endif

    } else if (marg_flag == MargImagSecondNew) {

        SaveLocalPos(image_count - 2, image_count - 3);

        Eigen::Vector3d P0 = Ps[(image_count - 2)];
        Eigen::Matrix3d R0 = Rs[(image_count - 2)];
        Eigen::Vector3d P1 = Ps[(image_count - 1)];
        Eigen::Matrix3d R1 = Rs[(image_count - 1)];
        SlideWindowFrame((image_count - 2), image_count, 1);
        SlideWindowNew(P0, R0, P1, R1, TIC[0], RIC[0]);
        image_count--;
    }
    LOG_OUT << "marge time:" << t_marg.toc() << std::endl;

#if USE_ASSERT
    if (!USE_GNSS) {
        for (int i = 0; i < (int)last_marg_info->keep_block_addr.size(); i++) {
            double* pointer = last_marg_info->keep_block_addr[i];
            if (last_marg_info->keep_block_size[i] == 1) {
                bool found = false;
                for (auto& it_per_id : f_manager.feature) {
                    if (!it_per_id.valid)continue;
                    bool condition = (idepths_all[0].find(it_per_id.feature_id) != idepths_all[0].end() && pointer == &idepths_all[0][it_per_id.feature_id]) || pointer == &scale_factor || pointer == &blackvalue;
                    if (ESTIMATE_TD)
                        condition = condition || (pointer == &td);
                    if (condition)
                        found = true;
                }
                ASSERT(found);
            }
        }
        for (int i = 0; i < (int)idepths_all.size(); i++) {
            for (auto it = idepths_all[i].begin(); it != idepths_all[i].end(); it++)
                ASSERT(it->second != 0);
        }
    }

#endif

}




void SWFOptimization::SaveLocalPos(int m_index, int r_index) {
    if (!enable_output)return;
    PosInfo pos_info;

    Eigen::Vector3d Pc_m = Ps[m_index] + Rs[m_index] * TIC[0] - P_WI_WC;
    Eigen::Matrix3d Rc_m = Rs[m_index] * RIC[0];
    Pc_m = R_WI_WC.transpose() * (Pc_m / scale_factor);
    Rc_m = R_WI_WC.transpose() * Rc_m;

    Eigen::Vector3d Pc_r = Ps[r_index] + Rs[r_index] * TIC[0] - P_WI_WC;
    Eigen::Matrix3d Rc_r = Rs[r_index] * RIC[0];
    Pc_r = R_WI_WC.transpose() * (Pc_r / scale_factor);
    Rc_r = R_WI_WC.transpose() * Rc_r;

    pos_info.t = Rc_r.transpose() * (Pc_m - Pc_r);
    pos_info.R = Rc_r.transpose() * Rc_m;
    // pos_info.scale_factor = scale_factor;
    pos_info.time_stamp = headers[m_index];
    pos_save[para_pose[r_index]].push_back(pos_info);
}

void SWFOptimization::SaveKefPos(int r_index) {
    if (!enable_output)return;
    PosInfo pos_info;
    Eigen::Vector3d Pc_r = Ps[r_index] + Rs[r_index] * TIC[0] - P_WI_WC;
    Eigen::Matrix3d Rc_r = Rs[r_index] * RIC[0];
    Pc_r = R_WI_WC.transpose() * (Pc_r / scale_factor);
    Rc_r = R_WI_WC.transpose() * Rc_r;
    pos_info.R = Rc_r;
    pos_info.t = Pc_r;
    pos_info.time_stamp = headers[r_index];
    if (pos_save_all.size())
        ASSERT(pos_save_all[pos_save_all.size() - 1].time_stamp < pos_info.time_stamp);
    pos_save_all.push_back(pos_info);
    if (pos_save.find(para_pose[r_index]) != pos_save.end()) {
        auto& tmp = pos_save[para_pose[r_index]];
        for (int j = 0; j < (int)tmp.size(); j++) {

            PosInfo pos_info;
            pos_info.time_stamp = tmp[j].time_stamp;
            pos_info.R = Rc_r * tmp[j].R;
            pos_info.t = Rc_r * tmp[j].t + Pc_r;
            if (pos_save_all.size())
                ASSERT(pos_save_all[pos_save_all.size() - 1].time_stamp < pos_info.time_stamp);
            pos_save_all.push_back(pos_info);
        }
        pos_save.erase(para_pose[r_index]);
    }
}
std::vector<PosInfo> SWFOptimization::RetriveAllPose() {

    for (int i = 0; i < image_count; i++)
        SaveKefPos(i);
    return pos_save_all;
}


//main process
void SWFOptimization::MeasurementProcess() {


    while ((USE_GNSS ? (!feature_buf.empty() && !acc_buf.empty() && !rover_buf.empty()) : (!feature_buf.empty() && !acc_buf.empty()))) {

        TicToc t_process;
        pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>> feature;

        feature = feature_buf.front();
        if (USE_GNSS) {
            mea_t* rover = rover_buf.front();
            if (image_count==0&& getVariableUseSingleNum(rover) < 8) {
                rover_buf.pop_front();
                while (rover->ros_time > feature_buf.front().first)
                    feature_buf.pop();
                return;
            }


            if (!rover_buf.empty() && rover && rover->ros_time < feature.first) {
                cur_time = rover->ros_time;
                // std::cout << rover->ros_time << std::endl;
                if (!ImuAvailable(cur_time)) return;
                rover_buf.pop_front();
                if (getVariableUseSingleNum(rover) < 8)continue;
                if (image_count == 0)continue;
                ImuIntegrate();
                gnss_imu_infos[image_count]->AddGNSSInfo(rover, Ps[image_count], Vs[image_count], Rs[image_count], Bas[image_count], Bgs[image_count]);
                if (USE_RTK)needNreaolve = true;
                double pos[3];
                Eigen::Matrix3d Rwwg;
                Eigen::Vector3d xyz = base_pos + gnss_imu_infos[image_count]->Pgs[0];
                ecef2pos(xyz.data(), pos);
                xyz2enu(pos, Rwwg.data());//transform vector from ecef to enu
                Rwgw = Rwwg.transpose();


                if (!init_gnss) {

                    ASSERT(gnss_imu_infos[image_count]->Vgs.size() == 1);
                    Eigen::Vector3d Vg = gnss_imu_infos[image_count]->Vgs[0];
                    Eigen::Vector3d Pg = gnss_imu_infos[image_count]->Pgs[0];

                    InitializePose();
                    for (int i = 0; i < image_count + 1; i++) {
                        Rs[i] = Rwgw * Rs[i];
                        Ps[i] = Pg - Rs[i] * Pbg;
                        Vs[i] = Vg - Rs[i] * (Utility::skewSymmetric(gnss_imu_infos[i]->pre_integrations[0]->gyrj - Bgs[i]) * Pbg);
                    }
                    ASSERT(scale_factor == 1);

                    P_WI_WC = Pg - Rs[0] * Pbg + Rs[0] * TIC[0];
                    // P_WI_WC.setZero();
                    gnss_imu_infos[image_count]->Rs[0] = Rs[image_count];
                    gnss_imu_infos[image_count]->Vector2Double();
                    old_P = Ps[image_count];



                    InitializeSqrtInfo();
                    init_gnss = true;
                }

                continue;
            } else {
                if (rover->ros_time - cur_time <= 0.01) {
#if REAL_TIME
                    mutex_feature.lock();
#endif
                    feature_buf.pop();
#if REAL_TIME
                    mutex_feature.unlock();
#endif
                    continue;
                }
            }

        }
#if REAL_TIME
        if (feature_buf.size() > 20) {
            mutex_feature.lock();
            feature_buf.pop();
            std::cout << "pop\r\n" << std::endl;
            feature = feature_buf.front();
            mutex_feature.unlock();

        }
#endif
        cur_time = feature.first;

        if (!ImuAvailable(cur_time))
            return;
#if REAL_TIME
        mutex_feature.lock();
#endif
        feature_buf.pop();
#if REAL_TIME
        mutex_feature.unlock();
#endif
        ImuIntegrate();
        image_count++;


        headers[image_count - 1] = cur_time;
        ASSERT(fabs(headers[image_count - 1] - cur_time) < 1e-4);
        ASSERT(image_count <= SWF_WINDOW_SIZE + 1);
        ImagePreprocess(feature.second);

        if (solver_flag != Initial )
            MyOptimization();

        ImagePostprocess();
        if (enable_output)
            save_result(*this, headers[image_count - 1]);

        SlideWindow();

#if USE_ASSERT
        if (last_marg_info && !USE_GNSS)
            for (int i = 0; i < (int)last_marg_info->keep_block_addr.size(); i++) {
                double* pointer = last_marg_info->keep_block_addr[i];
                if (last_marg_info->keep_block_size[i] == 1) {
                    bool found = false;
                    for (auto& it_per_id : f_manager.feature) {
                        if (!it_per_id.valid)continue;
                        bool condition = (idepths_all[0].find(it_per_id.feature_id) != idepths_all[0].end() && pointer == &idepths_all[0][it_per_id.feature_id]) || pointer == &scale_factor || pointer == &blackvalue;
                        if (ESTIMATE_TD)
                            condition = condition || (pointer == &td);
                        if (condition)
                            found = true;
                    }
                    ASSERT(found);
                }
            }

#endif

        Ps[image_count] = Ps[image_count - 1];
        Rs[image_count] = Rs[image_count - 1];
        Vs[image_count] = Vs[image_count - 1];
        Bas[image_count] = Bas[image_count - 1];
        Bgs[image_count] = Bgs[image_count - 1];
        headers[image_count] = headers[image_count - 1];

        if (solver_flag != Initial) {
            static double travel_distance = 0;
            travel_distance += (Ps[image_count - 1] - old_P).norm();
            if (!has_excitation && travel_distance > 1 && Vs[image_count - 1].norm() > 0.1 && image_count > 11) {
                has_excitation = true;
                enable_output = true;
                if (!have_hist) {
                    if (ESTIMATE_SCALE)fix_scale = true;
                    if (visual_inertial_bases_global[0]) visual_inertial_bases_global[0]->ResetInit();
                }
            }
        }
        if (have_hist)enable_output = true;
        old_P = Ps[image_count - 1];


        if (solver_flag != Initial)
            PubData();

        {
            int rover_count = 0;
            for (int i = 0; i < (int)gnss_imu_infos.size(); i++) {
                if (gnss_imu_infos[i])
                    rover_count += gnss_imu_infos[i]->Rovers.size();
            }
            static double t_process2 = 0;
            static int t_count = 0;
            double ts = t_process.toc();
            t_process2 += ts;
            t_count += 1;
            printf("process measurement time: %f   ,%f   ,%f  ,\t %d\t %d \n", headers[image_count], ts, t_process2 / t_count, (int)feature_buf.size(), rover_count);
            LOG_OUT << "process measurement time: " << ts << "," << t_process2 / t_count << std::endl << std::endl;
            LOG_OUT.flush();

            // static bool f;
            // if (!f) {
            //     f = true;
            //     std::ofstream fout(RESULT_PATH + "_time", std::ios::out);
            //     fout << "time1,time2";
            //     fout << std::endl;
            //     fout.close();
            // }


            // ofstream foutC(RESULT_PATH + "_time", ios::app);
            // foutC << ts << "," << t_process2 / t_count << std::endl;

        }

        if (last_marg_info) {
            int count = 0;
            for (int i = 0; i < (int)last_marg_info->keep_block_addr.size(); i++) {
                ASSERT(!(last_marg_info->keep_block_size[i] == 9 && last_marg_info->keep_block_addr[i] != para_speed_bias[0]));
                if (last_marg_info->keep_block_size[i] == 1)
                    count++;
            }
            LOG_OUT << "prior:" << count << std::endl;
            for (int i = 0; i < (int)idepths_all.size(); i++) {
                for (auto it = idepths_all[i].begin(); it != idepths_all[i].end(); it++)
                    ASSERT(it->second != 0);
            }
        }
    }

}




//publicating and saving results.
void SWFOptimization::PubData() {

#if DEBUG
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(headers[image_count - 1]);
    if (!pub_init && solver_flag == NonLinear) {
        resetpot(*this, header);
        pub_init = true;
    }
    printStatistics(*this, 0);

    if (enable_output) {
        pubOdometry(*this, header);
        pubCameraPose(*this, header);
        pubPointCloud(*this, header);
    }
#endif
}

