#include "swf.h"
#include "ceres/residual_block.h"
#include "../factor/projection_factor.h"
#include "../factor/mag_factor.h"
#include "../factor/pose0_factor.h"
void SWFOptimization::InputImage(double t, const cv::Mat& _img, const cv::Mat& _img1) {
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    static int inputImageCnt;

    if (inputImageCnt < AVERAGE_IMAGE) {
        inputImageCnt++;
        return;
    }

    if (USE_IMAGE) {
        TicToc t_process;
        if (_img1.empty())
            featureFrame = feature_tracker.trackImage(t, _img);
        else
            featureFrame = feature_tracker.trackImage(t, _img, _img1);

        feature_buf.push(make_pair(t, featureFrame));
    }
    MeasurementProcess();

}
void SWFOptimization::ImagePreprocess(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image) {
    assert(frame_types[image_count + rover_count - 1] == ImagFrame);
    if (f_manager.addFeatureCheckParallax(image_count - 1, image))imag_marg_index = 0;
    else  imag_marg_index = image_count - 2;

    f_manager.removeOut(image_count, my_problem);
    if ((solver_flag == Initial || !USE_IMU) && frame_types[image_count + rover_count - 1] == ImagFrame) {
        transformdataimage(Frame2IR, ImagFrame);
        f_manager.initFramePoseByPnP(image_count - 1, Ps_image, Rs_image, tic, ric, Pbg);
        transformdataimage(IR2Frame, ImagFrame);
    }
    if (frame_types[image_count + rover_count - 1] == ImagFrame) {
        transformdataimage(Frame2IR, ImagFrame);
        f_manager.triangulate( Ps_image, Rs_image, tic, ric, Pbg);
    }

    if (solver_flag == Initial) {
        if (frame_types[image_count + rover_count - 1] == ImagFrame) {
            if (image_count == 1 + FEATURE_CONTINUE) {
                if (!init_gnss)
                    InitializeSqrtInfo();
                solver_flag = NonLinear;
                Vector2Double();
                if (USE_GLOBAL_OPTIMIZATION) {
                    AddFeature2Problem();
                } else {
                    ceres::Problem problem;
                    ceres::Solver::Options options;
                    AddAllResidual(NormalMode, std::set<double*> {}, nullptr, problem, options, true, true, true);
                }
                Double2Vector();
                printf("Initialization finish!");
            }
        }
    }
}


void SWFOptimization::AddFeature2Problem() {
    for (auto& it_per_id : f_manager.feature) {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < FEATURE_CONTINUE)continue;
        int imu_j = it_per_id.start_frame ;
#if USE_INVERSE_DEPTH
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;
#endif
        for (auto& it_per_frame : it_per_id.feature_per_frame) {


            if (!it_per_frame.is_in_optimize) {
                it_per_frame.is_in_optimize = true;
#if USE_INVERSE_DEPTH
                if (it_per_id.start_frame != imu_j) {
                    ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
                    Vector3d pts_j = it_per_frame.point;
                    ProjectionTwoFrameOneCamFactor* factor = new ProjectionTwoFrameOneCamFactor(pts_i, pts_j);
                    my_problem.AddResidualBlock(factor, loss_function,  para_pose[i2f[it_per_id.start_frame]], para_pose[i2f[imu_j]], para_ex_Pose[0], &it_per_id.idepth_);
                }
                if (USE_STEREO && it_per_frame.is_stereo) {
                    Vector3d pts_j_right = it_per_frame.pointRight;
                    if (it_per_id.start_frame != imu_j) {
                        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
                        ProjectionTwoFrameTwoCamFactor* factor = new ProjectionTwoFrameTwoCamFactor(pts_i, pts_j_right);
                        my_problem.AddResidualBlock(factor, loss_function,  para_pose[i2f[it_per_id.start_frame]], para_pose[i2f[imu_j]], para_ex_Pose[0], para_ex_Pose[1], &it_per_id.idepth_);
                    } else {
                        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
                        ProjectionOneFrameTwoCamFactor* factor = new ProjectionOneFrameTwoCamFactor(pts_i, pts_j_right);
                        my_problem.AddResidualBlock(factor, loss_function,  para_ex_Pose[0], para_ex_Pose[1], &it_per_id.idepth_);
                    }
                }
#else
                ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
                projection_factor* factor = new projection_factor(it_per_frame.point);
                my_problem.AddResidualBlock(factor, loss_function, para_pose[i2f[imu_j]], para_ex_Pose[0], it_per_id.ptsInWorld.data());
                if (USE_STEREO && it_per_frame.is_stereo) {
                    ceres::LossFunction* loss_function = new ceres::CauchyLoss(1.0);
                    projection_factor* factor = new projection_factor(it_per_frame.pointRight);
                    my_problem.AddResidualBlock(factor, loss_function, para_pose[i2f[imu_j]],  para_ex_Pose[1], it_per_id.ptsInWorld.data());
                }
#endif



            }
            imu_j++;
        }
    }
}
void SWFOptimization::ImagePostprocess() {
    if (solver_flag != Initial) {
        set<int> removeIndex;
        OutliersRejection(removeIndex, image_count);
        f_manager.removeFailures(my_problem);

    }

    if (imag_marg_index == 0 || f_manager.CheckParallax(image_count - 1))imag_marg_index = 0;
    else  imag_marg_index = image_count - 2;


}

int SWFOptimization::CheckProblem() {

    int constnum = 0;
    for (int i = 0; i < image_count + 1; i++) {
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
                if (index2 - index1 > 1)
                    std::cout << "GNSS multiple frame\r\n";
            }
        } else {
            index1 = 0;
            index2 = rover_count - 1;
        }
        if (index2 - index1 > 1) {
            assert(imu_gnss_factor[i]);
            assert((int)(imu_gnss_factor[i]->gnss_poses.size()) == index2 - index1 - 1);
            for (int j = index1 + 1; j < index2; j++) {
                assert(!my_problem.HasParameterBlock(para_speed_bias[j]));
                assert(!my_problem.HasParameterBlock(para_pose[j]));

            }
        }
        assert(my_problem.HasParameterBlock(para_speed_bias[index1]));
        assert(my_problem.HasParameterBlock(para_pose[index1]));
        assert(my_problem.HasParameterBlock(para_speed_bias[index2]));
        assert(my_problem.HasParameterBlock(para_pose[index2]));

    }

    for (int i = 0; i < NUM_OF_CAM; i++) {
        if ((ESTIMATE_EXTRINSIC && Vs[g2f[rover_count - 1]].norm() > 0.2) || open_ex_estimation) {
            open_ex_estimation = 1;
            if (my_problem.HasParameterBlock(para_ex_Pose[i]))
                my_problem.SetParameterBlockVariable(para_ex_Pose[i]);
        } else {
            if (my_problem.HasParameterBlock(para_ex_Pose[i])) {
                my_problem.SetParameterBlockConstant(para_ex_Pose[i]);
                constnum++;
            }

        }
    }
    if (USE_MAG_CORRECT_YAW) {
        if (frame_types[rover_count + image_count - 1] == ImagFrame) {
            if (my_problem.HasParameterBlock(para_bmg)) {
                my_problem.SetParameterBlockConstant(para_bmg);
                constnum++;
            }
        } else {
            if (my_problem.HasParameterBlock(para_bmg))
                my_problem.SetParameterBlockVariable(para_bmg);
        }
    }


    return constnum;

}

void SWFOptimization::MyOptimization( ) {
    if (USE_DIRECT_N_RESOLVE && need_Nresolve)
        UpdateNParameterHead();
    if (USE_GLOBAL_OPTIMIZATION) {
        Vector2Double();
        ceres::Solver::Summary summary;
        ceres::internal::is_optimize = true;


        int constnum = CheckProblem();
        AddFeature2Problem();



        MyOrdering(my_problem, my_options);
        my_options.max_num_iterations = MAX_NUM_ITERATIONS;



        TicToc t_marg;
        static double t, count;
        ceres::Solve(my_options, &my_problem, &summary);
        if (summary.final_cost > 1e10) {
            ceres::Solve(my_options, &my_problem, &summary);
        }
        assert(summary.final_cost < 1e10);


        t += t_marg.toc();
        count += 1;
        std::cout << "solving time:" << t_marg.toc() << "," << t / count << std::endl;

        std::cout << summary.BriefReport() << std::endl;
        assert(my_options.linear_solver_ordering);
        if (my_options.linear_solver_ordering) {
            std::cout << my_options.linear_solver_ordering->NumElements() << "," << constnum << "," << my_problem.NumParameterBlocks() << std::endl;
            assert(my_options.linear_solver_ordering->NumElements() + constnum == my_problem.NumParameterBlocks());
        }


        if (USE_DIRECT_N_RESOLVE && need_Nresolve)
            UpdateSchurHessianOnly(my_problem);
        Double2Vector();
    } else {
        Vector2Double();
        ceres::Problem problem;
        ceres::Solver::Options options;
        AddAllResidual(NormalMode, std::set<double*> {}, nullptr, problem, options, true, true, true);

        if (USE_DIRECT_N_RESOLVE && need_Nresolve)
            UpdateSchurHessianOnly(problem);
        Double2Vector();
    }
}



double SWFOptimization::ReprojectionError(Matrix3d& Rj, Vector3d& Pj, Matrix3d& ricj, Vector3d& ticj, Vector3d pts_w, Vector3d& uvj) {
    Vector3d pts_cj = ricj.transpose() * (Rj.transpose() * (pts_w - Pj) - ticj);
    Vector2d residual = (pts_cj / pts_cj.z()).head<2>() - uvj.head<2>();
    double rx = residual.x();
    double ry = residual.y();
    return sqrt(rx * rx + ry * ry);
}



void SWFOptimization::OutliersRejection(set<int>& removeIndex, int image_count) {

    for (auto& it_per_id : f_manager.feature) {
        double err = 0;
        int errCnt = 0;
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < FEATURE_CONTINUE)
            continue;
        int imu_j = it_per_id.start_frame - 1;
#if USE_INVERSE_DEPTH
        Eigen::Vector3d ptsInW = Rs[i2f[it_per_id.start_frame]] * (
                                     ric[0] * (it_per_id.feature_per_frame[0].point / it_per_id.idepth_) + tic[0] - Pbg
                                 ) + Ps[i2f[it_per_id.start_frame]];
#else
        Eigen::Vector3d ptsInW = it_per_id.ptsInWorld;
#endif

        for (auto& it_per_frame : it_per_id.feature_per_frame) {
            imu_j++;
            Eigen::Vector3d P = Ps[i2f[imu_j]] - Rs[i2f[imu_j]] * Pbg;



            double tmp_error = ReprojectionError(Rs[i2f[imu_j]], P, ric[0], tic[0],
                                                 ptsInW, it_per_frame.point);
            err += tmp_error;
            errCnt++;
            // need to rewrite projecton factor.........
            if (it_per_frame.is_stereo && USE_STEREO) {
                double tmp_error = ReprojectionError(Rs[i2f[imu_j]], P, ric[1], tic[1],
                                                     ptsInW, it_per_frame.pointRight);
                err += tmp_error;
                errCnt++;
            }
        }
        double ave_err = err / errCnt;
        if (ave_err * FOCAL_LENGTH / FEATUREWEIGHTINVERSE > 2) {
            it_per_id.solve_flag = 2;
            removeIndex.insert(it_per_id.feature_id);
        }

    }

}




void SWFOptimization::SlideWindowOld(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0, Eigen::Vector3d Pbg) {
    f_manager.removeBack(P0, R0, P1, R1, tic0, ric0, Pbg, my_problem);
}


void SWFOptimization::SlideWindowNew() {
    f_manager.removeFront(image_count - 1, my_problem);
}


void SWFOptimization::transformdataimage(Direction direction, int mode) {
    assert(mode == ImagFrame);
    if (direction == Frame2IR) {
        for (int i = 0; i < image_count; i++) {
            int j = ImageRoverId2FrameId(i, mode);
            Ps_image[i] = Ps[j];
            Vs_image[i] = Vs[j];
            Rs_image[i] = Rs[j];
        }

    } else {
        for (int i = 0; i < image_count; i++) {
            int j = ImageRoverId2FrameId(i, mode);
            Ps[j] = Ps_image[i];
            Vs[j] = Vs_image[i];
            Rs[j] = Rs_image[i];
        }
    }
}

void SWFOptimization::GlobalMarge(std::set<double*>MargePoints) {

    TicToc t_marg2;
    MarginalizationInfo* marginalization_info = new MarginalizationInfo();
    std::vector<ceres::internal::ResidualBlock*>residual_blocks_all;
    std::set<double*>parameter_head;
    std::vector<double*>parameter_blocks_all;
    my_problem.GetParameterBlocks(&parameter_blocks_all);
    for (int i = 0; i < (int)parameter_blocks_all.size(); i++)my_problem.SetParameterBlockConstant(parameter_blocks_all[i]);
    my_problem.GetResidualBlocks(&residual_blocks_all);
    for (int i = 0; i < (int)residual_blocks_all.size(); i++) {
        if (residual_blocks_all[i] != marg_residual_block_id)
            residual_blocks_all[i]->is_use = false;
        else
            residual_blocks_all[i]->is_use = true;
    }

    for (auto it = MargePoints.begin(); it != MargePoints.end(); it++) {
        std::vector<ceres::internal::ResidualBlock*>residual_blocks;
        if (my_problem.HasParameterBlock(*it))my_problem.GetResidualBlocksForParameterBlock(*it, &residual_blocks);
        for (int i = 0; i < (int)residual_blocks.size(); i++) {
            if (residual_blocks[i]->is_use)continue;
            residual_blocks[i]->is_use = true;
            std::vector<double*> parameter_blocks;
            my_problem.GetParameterBlocksForResidualBlock(residual_blocks[i], &parameter_blocks);
            for (int j = 0; j < (int)parameter_blocks.size(); j++) {
                if (MargePoints.find(parameter_blocks[j]) != MargePoints.end())continue;
                parameter_head.insert(parameter_blocks[j]);
            }
        }
    }

    if (last_marg_info) {
        for (int j = 0; j < (int)last_marg_info->keep_block_addr.size(); j++) {
            if (MargePoints.find(last_marg_info->keep_block_addr[j]) != MargePoints.end())continue;
            parameter_head.insert(last_marg_info->keep_block_addr[j]);
        }
    }

    for (auto it = parameter_head.begin(); it != parameter_head.end(); it++) {
        my_problem.SetParameterBlockVariable(*it);
        ceres::internal::parameter_head.push_back(*it);
    }
    for (auto it = MargePoints.begin(); it != MargePoints.end(); it++) {
        if (my_problem.HasParameterBlock(*it))my_problem.SetParameterBlockVariable(*it);
    }


    if (ceres::internal::parameter_head.size() == 0) {
        if (last_marg_info) {
            delete last_marg_info;
            if (USE_GLOBAL_OPTIMIZATION) {
                my_problem.RemoveResidualBlock(marg_residual_block_id);
            }
        }
        last_marg_info = 0;
        printf("warning!!no last info\r\n");
    } else {
        my_options.max_num_iterations = 1;
        my_options.jacobi_scaling = false;
        // my_options.initial_trust_region_radius=my_options.max_trust_region_radius=MAX_TRUST_REGION_RADIUS;
        ceres::internal::is_optimize = false;
        ceres::Solver::Summary summary;
        TicToc t_marg;
        MyOrdering(my_problem, my_options);
        ceres::Solve(my_options, &my_problem, &summary);
        std::cout << "solving time:" << t_marg.toc() << std::endl;

        my_options.max_num_iterations = MAX_NUM_ITERATIONS;
        my_options.jacobi_scaling = 0;
        // my_options.initial_trust_region_radius=my_options.max_trust_region_radius=MAX_TRUST_REGION_RADIUS;

        UpdateSchur(my_problem);
        ceres::internal::is_optimize = true;

        marginalization_info->setmarginalizeinfo(parameter_block_addr, parameter_block_global_size, A, b, true);
        marginalization_info->getParameterBlocks();
    }

    for (int i = 0; i < (int)residual_blocks_all.size(); i++)residual_blocks_all[i]->is_use = true;
    for (int i = 0; i < (int)parameter_blocks_all.size(); i++)my_problem.SetParameterBlockVariable(parameter_blocks_all[i]);
    if (last_marg_info) {
        delete last_marg_info;
        if (USE_GLOBAL_OPTIMIZATION) {
            my_problem.RemoveResidualBlock(marg_residual_block_id);
        }
    }
    last_marg_info = marginalization_info;
    MarginalizationFactor* factor = new MarginalizationFactor(last_marg_info);
    marg_residual_block_id = my_problem.AddResidualBlock(factor, 0, last_marg_info->keep_block_addr);
}




