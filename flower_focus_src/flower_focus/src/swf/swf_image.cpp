#include "swf.h"
#include <thread>
#include<fstream>
void SWFOptimization::InputImage(double t, const cv::Mat& _img, const cv::Mat& _img1) {
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;

    TicToc t_process;

    if (_img1.empty())
        featureFrame = feature_tracker.trackImage(t, _img);
    else
        featureFrame = feature_tracker.trackImage(t, _img, _img1);

    // std::cout<<t_process.toc()<<std::endl;
    static int t_count;
    static double t_sum;
    t_sum += t_process.toc();
    t_count += 1;
#if REAL_TIME
    mutex_feature.lock();
    mutex_imu_feature.lock();
#endif

    if (first_observe_time != 0 && t >= first_observe_time + SKIP_TIME + AVERAGE_TIME + 0.1)
        feature_buf.push(make_pair(t, featureFrame));

#if REAL_TIME
    mutex_imu_feature.unlock();
    mutex_feature.unlock();
#endif
#if !REAL_TIME
    MeasurementProcess();
#endif

}




void SWFOptimization::ImagePreprocess(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image) {

    f_manager.addFeature(image_count - 1, image);
    f_manager.removeOut(image_count);
    f_manager.triangulate( Ps, Rs);

    if (solver_flag == Initial) {
        if (image_count == FEATURE_CONTINUE + 1) {
            if (USE_GNSS && !init_gnss)return;
            if (!init_gnss) {
                ASSERT(!USE_GNSS);
                InitializeSqrtInfo();
                printf("Initialization finish!");
            }

        }
    }
}



void SWFOptimization::ImagePostprocess() {
    if (solver_flag != Initial) {
        OutliersRejection( image_count);
        std::set<double*>feature_marge_pointer;
        std::set<double*>prior_pointer;
        prior_pointer.insert(last_marg_info->keep_block_addr.begin(), last_marg_info->keep_block_addr.end());
        f_manager.GetFailurePointers(feature_marge_pointer, prior_pointer);


        if (!feature_marge_pointer.empty()) {
            last_marg_info->marginalize_pointers(feature_marge_pointer, true);
            if (USE_REPARAM) {
                set_last_marg_info_idx(0);
                last_marg_info->UpdateF();
            }

            if (visual_inertial_bases_global[0])
                visual_inertial_bases_global[0]->ResetInit();
        }

        f_manager.removeFailures(visual_inertial_bases_global);
    }

    if (f_manager.CheckParallax(image_count - 1)) imag_marg_index = 0;
    else imag_marg_index = image_count - 2;


}

void SWFOptimization::MyOptimization( ) {

    TicToc t_marg;

    Vector2Double();
    OptimizationOrMarginalization(NormalMode);
    Double2Vector();
}



double SWFOptimization::ReprojectionError(Matrix3d& Rj, Vector3d& Pj, Matrix3d& ricj, Vector3d& ticj, Vector3d pts_w, Vector3d& uvj, double weight) {
    Vector3d pts_cj = ricj.transpose() * (Rj.transpose() * (pts_w - Pj) - ticj);
    Vector2d residual = (pts_cj / pts_cj.z()).head<2>() - uvj.head<2>();
    double rx = residual.x();
    double ry = residual.y();
    return sqrt(rx * rx + ry * ry) * weight;
}



void SWFOptimization::OutliersRejection(int image_count) {

    std::vector<FeaturePerId>new_features;
    std::set<int>outlier_id;
    int evaluate_frame_begin = (image_count - FRAME_NUM_FOR_REJECTION) / SWF_SIZE_IN * SWF_SIZE_IN;
    int count_left = 0, count_right = 0;
    for (auto& it_per_id : f_manager.feature) {
        for (auto& it_per_frame : it_per_id.feature_per_frame) {
            if (it_per_frame.is_stereo)count_right++;
            count_left++;
        }

    }
    LOG_OUT << "count_left:" << count_left << "\t" << "count_right:" << count_right << std::endl;


    for (auto& it_per_id : f_manager.feature) {
        if (!it_per_id.valid)continue;

        count_left++;

        double err = 0;
        int errCnt = 0;
        int imu_j = it_per_id.start_frame - 1;
        Eigen::Vector3d ptsInW;


        bool is_cross = (it_per_id.start_frame / SWF_SIZE_IN != (it_per_id.endFrame() - 1) / SWF_SIZE_IN);
        int imu_i;
        if (!is_cross)
            imu_i = it_per_id.start_frame;
        else {
            if (it_per_id.endFrame() >= evaluate_frame_begin) {
                if (it_per_id.start_frame >= evaluate_frame_begin) {
                    if ( it_per_id.start_frame % SWF_SIZE_IN != 0)imu_i = (it_per_id.start_frame / SWF_SIZE_IN + 1) * SWF_SIZE_IN;
                    else imu_i = it_per_id.start_frame;
                } else
                    imu_i = evaluate_frame_begin;
            } else
                imu_i = it_per_id.start_frame;
        }

        if (it_per_id.endFrame() > evaluate_frame_begin)ASSERT(imu_i >= evaluate_frame_begin);
        ASSERT(imu_i >= it_per_id.start_frame && imu_i <= it_per_id.endFrame());

        if (imu_i < evaluate_frame_begin || it_per_id.endFrame() <= evaluate_frame_begin)continue;

        if (is_cross)ASSERT(imu_i % SWF_SIZE_IN == 0);

        double idepth = idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id];
        ptsInW = Rs[imu_i] * (RIC[0] * (it_per_id.feature_per_frame[imu_i - it_per_id.start_frame].point / idepth) + TIC[0]) + Ps[imu_i];

        ASSERT(idepth != 0);

        for (auto& it_per_frame : it_per_id.feature_per_frame) {
            imu_j++;
            if (imu_j < evaluate_frame_begin)continue;

            double tmp_error = ReprojectionError(Rs[imu_j],  Ps[imu_j], RIC[0], TIC[0],
                                                 ptsInW, it_per_frame.point,
                                                 SQRT_INFO_FEATURE
                                                );
            if (tmp_error  > OUTLIER_SHRESHOLD1 * 3) {
                err += 10000;
                errCnt++;
                break;
            }
            err += tmp_error;
            errCnt++;
#if USE_STEREO
            if (it_per_frame.is_stereo && STEREO_SCALE) {
                double tmp_error = ReprojectionError(Rs[imu_j],  Ps[imu_j], RIC[1], TIC[1],
                                                     ptsInW, it_per_frame.pointRight,
                                                     SQRT_INFO_FEATURE
                                                    );
                if (tmp_error > OUTLIER_SHRESHOLD1 * 3) {
                    it_per_frame.is_stereo = false;
                    continue;
                }
            }
#endif
        }
        ASSERT(errCnt > 0);

        double ave_err = err / errCnt;

        if (ave_err > OUTLIER_SHRESHOLD1) {
            outlier_id.insert(it_per_id.feature_id);
            it_per_id.solve_flag = 2;
            if (it_per_id.start_frame < evaluate_frame_begin && evaluate_frame_begin > 0) {
                int remain_size = evaluate_frame_begin - it_per_id.start_frame + 1;
                if (remain_size >= FEATURE_CONTINUE)RegenerateFeatures(it_per_id, remain_size, new_features);
            }
        }
    }
    for (int i = 0; i < (int)new_features.size(); i++) {
        f_manager.feature.push_back(new_features[i]);
        f_manager.feature_map[f_manager.feature.back().feature_id] = &f_manager.feature.back();
    }
}

void SWFOptimization::RegenerateFeatures(FeaturePerId& it_per_id, int remain_size, std::vector<FeaturePerId>& new_features) {
    static int new_id = -1;

    FeaturePerId tmp(new_id--, it_per_id.start_frame);
    tmp.valid = 1;
    tmp.solve_flag = 1;
    ASSERT(remain_size < (int)it_per_id.feature_per_frame.size());
    for (int ij = 0; ij < remain_size; ij++)
        tmp.feature_per_frame.push_back(it_per_id.feature_per_frame[ij]);
    int found = 0;
    for (int ij = 0; ij < tmp.endFrame() / SWF_SIZE_IN + 1; ij++) {
        if (idepths_all[ij].find(it_per_id.feature_id) != idepths_all[ij].end()) {
            idepths_all[ij][tmp.feature_id] = idepths_all[ij][it_per_id.feature_id];
            if (tmp.feature_per_frame.size() >= FEATURE_CONTINUE) {
                double* old_pointer = &idepths_all[ij][it_per_id.feature_id];
                double* new_pointer = &idepths_all[ij][tmp.feature_id];
                if (last_marg_info->keep_block_addr_set.find(old_pointer) != last_marg_info->keep_block_addr_set.end()) {
                    found++;
                    int idx = last_marg_info->keep_block_addr_set[old_pointer];
                    last_marg_info->keep_block_addr[idx] = new_pointer;

                    last_marg_info->parameter_block_size[new_pointer] =
                        last_marg_info->parameter_block_size[old_pointer];
                    last_marg_info->parameter_block_idx[new_pointer] =
                        last_marg_info->parameter_block_idx[old_pointer];
                    last_marg_info->parameter_block_data[new_pointer] =
                        last_marg_info->parameter_block_data[old_pointer];

                    last_marg_info->parameter_block_size.erase(old_pointer);
                    last_marg_info->parameter_block_idx.erase(old_pointer);
                    last_marg_info->parameter_block_data.erase(old_pointer);

                    ASSERT(ij == 0);
                    last_marg_info->keep_block_addr_set.clear();
                    for (int i = 0; i < (int)last_marg_info->keep_block_addr.size(); i++)
                        last_marg_info->keep_block_addr_set[last_marg_info->keep_block_addr[i]] = i;
                }
            }
        }

    }
    ASSERT(found <= 1);

    std::set<int>factor_index_set;
    for (int i = 0; i < (int)visual_inertial_bases_global.size(); i++) {
        if (visual_inertial_bases_global[i])factor_index_set.insert(i);
    }
    AddFeatures(tmp, visual_inertial_bases_global, factor_index_set);


    LOG_OUT << "new feature:" << it_per_id.feature_id << "," << tmp.feature_id << "," << it_per_id.feature_per_frame.size() << "," << tmp.feature_per_frame.size() << "," << tmp.start_frame << std::endl;
#if USE_ASSERT
    bool is_cross = (tmp.start_frame / SWF_SIZE_IN != (tmp.endFrame() - 1) / SWF_SIZE_IN) || tmp.start_frame == 0;
    int imu_i = tmp.start_frame;
    if (is_cross && tmp.start_frame % SWF_SIZE_IN != 0)
        imu_i = (tmp.start_frame / SWF_SIZE_IN + 1) * SWF_SIZE_IN;
    ASSERT(idepths_all[imu_i / SWF_SIZE_IN][tmp.feature_id] != 0);
#endif

    new_features.push_back(tmp);
}

void SWFOptimization::SlideWindowOld(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0) {
    f_manager.removeBack(P0, R0, P1, R1, tic0, ric0);
}


void SWFOptimization::SlideWindowNew(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0) {
    f_manager.removeFront(image_count - 1, P0, R0, P1, R1, tic0, ric0);
}





