#include "swf.h"
#include "../utility/visualization.h"
#include<fstream>
#include "../factor/initial_factor.h"

void SWFOptimization::InitializePose() {
    for (int i = 0; i < SWF_WINDOW_SIZE + 1; i++)
        Bgs[i].setZero();

    LOG_OUT << "initial bgs:" << Bgs[0];
    printf("averge acc %f %f %f\n", acc_mean.x(), acc_mean.y(), acc_mean.z());

    if (mag_mean.norm() < 1) {
        mag_mean(0) = 0;
        mag_mean(1) = 1;
        mag_mean(2) = 0;
    }


    Matrix3d Rwb0;
    Eigen::Vector3d z0 = acc_mean.normalized();
    Eigen::Vector3d x0 = (Utility::skewSymmetric(mag_mean) * z0).normalized();
    Eigen::Vector3d y0 = (Utility::skewSymmetric(z0) * x0).normalized();
    Rwb0.block(0, 0, 1, 3) = x0.transpose();
    Rwb0.block(1, 0, 1, 3) = y0.transpose();
    Rwb0.block(2, 0, 1, 3) = z0.transpose();


    for (int i = 0; i < SWF_WINDOW_SIZE + 1; i++)
        Rs[i] = Rwb0;
    LOG_OUT << "init R0: " << endl
            << Utility::R2ypr(Rs[0]).transpose() << ","
            << Utility::R2ypr(Rwb0).transpose() << ","
            << endl;
}



void SWFOptimization::InputIMU(double t, const Vector3d& linearAcceleration, const Vector3d& angularVelocity) {
    if (first_observe_time == 0) {
#if REAL_TIME
        mutex_imu_feature.lock();
#endif
        first_observe_time = t;
#if REAL_TIME
        mutex_imu_feature.unlock();
#endif
    }

    if (t < first_observe_time + SKIP_TIME)
        return;
    if (t < first_observe_time + SKIP_TIME + AVERAGE_TIME) {
        acc_mean += linearAcceleration;
        acc_count++;
        return;
    }
    static bool imu_initialize = false;
    if (!imu_initialize) {
        imu_initialize = true;
        acc_mean /= acc_count;
        InitializePose();
    }
#if REAL_TIME
    mutex_imu.lock();
#endif
    acc_buf.push(make_pair(t, linearAcceleration));
    gyr_buf.push(make_pair(t, angularVelocity));
#if REAL_TIME
    mutex_imu.unlock();
#endif
}


bool SWFOptimization::GetImuInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>>& accVector,
                                     vector<pair<double, Eigen::Vector3d>>& gyrVector) {
    if (acc_buf.empty()) {
        printf("not receive imu\n");
        return false;
    }
    if (t1 <= acc_buf.back().first) {
#if REAL_TIME
        mutex_imu.lock();
#endif
        while (acc_buf.front().first <= t0) {
            acc_buf.pop();
            gyr_buf.pop();
        }
        while (acc_buf.front().first < t1) {
            accVector.push_back(acc_buf.front());
            acc_buf.pop();
            gyrVector.push_back(gyr_buf.front());
            gyr_buf.pop();
        }
        accVector.push_back(acc_buf.front());
        gyrVector.push_back(gyr_buf.front());
#if REAL_TIME
        mutex_imu.unlock();
#endif
    } else {
        printf("wait for imu\n");
        return false;
    }
    return true;
}


bool SWFOptimization::ImuAvailable(double t) {

    double tend = acc_buf.back().first;
    if (!acc_buf.empty() && t <= tend)
        return true;
    else
        return false;
}


void SWFOptimization::ImuIntegrate() {
    if (AMP) {
        static double ACC_N0 = ACC_N, INT_N0 = INT_N, ACC_W0 = ACC_W, GYR_N0 = GYR_N, GYR_W0 = GYR_W;
        if (image_count + 1 < SWF_WINDOW_SIZE / AMP_NUM) {
            ACC_N = ACC_N0 * AMP;
            ACC_W = ACC_W0 * AMP;
            GYR_N = GYR_N0 * AMP;
            GYR_W = GYR_W0 * AMP;
        } else {
            static bool is_first = true;
            if (is_first) {
                is_first = false;
                ACC_N = ACC_N0; INT_N = INT_N0; ACC_W = ACC_W0; GYR_N = GYR_N0; GYR_W = GYR_W0;
                for (int i = 0; i < (int)gnss_imu_infos.size(); i++) {
                    if (gnss_imu_infos[i])
                        gnss_imu_infos[i]->resetNoise(Bas[i], Bgs[i], acc_scale, gyr_scale);

                }
                if (visual_inertial_bases_global[0]) visual_inertial_bases_global[0]->ResetInit();
            }
        }
    }
    vector<pair<double, Eigen::Vector3d>> accVector, gyrVector;

    GetImuInterval(prev_time, cur_time, accVector, gyrVector);

    ASSERT(accVector.size() >= 1);
    ASSERT(cur_time >= prev_time);

    // ASSERT(accVector.size() > 2);
    double last_time_stamp = accVector[accVector.size() - 1].first;
    double second_last_time_stamp = accVector[accVector.size() - 2].first;
    int acc_size = accVector.size();
    if (abs(last_time_stamp - cur_time) > abs(second_last_time_stamp - cur_time)) {
        if (acc_size > 2)
            acc_size -= 1;
    }

    // ASSERT(acc_size >= 2);
    Ps[image_count] += Vs[image_count] * old_time_shift;
    Rs[image_count] *= Utility::deltaQ((gyr_0 - Bgs[image_count]) * old_time_shift).toRotationMatrix();
    static double old_imu_time = -1;
    IntegrationBase* pre_integration = 0;
    for (int i = 0; i < acc_size; i++) {
        double t = accVector[i].first;
        if (old_imu_time < 0) old_imu_time = t;
        if (t - old_imu_time == 0 && acc_size == 1 && image_count > 5)t = old_imu_time + 1e-3;
        double dt = t - old_imu_time;
        if (dt <= 0)continue;
        ASSERT(dt >= 0);
        old_imu_time = t;
        if (dt > 0) {
            IMUProcess(dt, accVector[i].second, gyrVector[i].second, pre_integration);
             
        }
    }
    ASSERT(pre_integration);
    double dt_j = old_imu_time - cur_time;
    double dt_i = old_time_shift;
    if (dt_i == 0)dt_i = 1e-10;
    if (dt_j == 0)dt_j = 1e-10;

    gnss_imu_infos[image_count]->Setdt(dt_i, dt_j);

    old_time_shift = dt_j;
    Ps[image_count] -= Vs[image_count] * old_time_shift;
    Rs[image_count] *= Utility::deltaQ(-(gyr_0 - Bgs[image_count]) * old_time_shift).toRotationMatrix();

    prev_time = cur_time;

}

void SWFOptimization::IMUProcess( double dt, const Vector3d& linear_acceleration, const Vector3d& angular_velocity, IntegrationBase*& pre_integration) {
    static bool first_imu = true;
    if (first_imu) {
        first_imu = false;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }
    if (!pre_integration) {
        pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[image_count], Bgs[image_count], acc_scale, gyr_scale};
        if (!gnss_imu_infos[image_count])
            gnss_imu_infos[image_count] = new GNSSIMUInfo(pre_integration);
        else
            gnss_imu_infos[image_count]->AddIMUinfo(pre_integration);

    }

    pre_integration->push_back(dt, linear_acceleration, angular_velocity);
    if (image_count != 0) {
        int j = image_count;
        Vector3d un_acc_0 = Rs[j] * ((acc_0.array() * acc_scale.array()).matrix() - Bas[j]) - Rwgw * G;
        Vector3d un_gyr = 0.5 * ((gyr_0 + angular_velocity).array() * gyr_scale.array()).matrix() - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * ((linear_acceleration.array() * acc_scale.array()).matrix() - Bas[j]) - Rwgw * G;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;

    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;

}






void SWFOptimization::InitializeSqrtInfo() {


    solver_flag = NonLinear;
    Vector2Double();
    MarginalizationInfo* marginalization_info = new MarginalizationInfo();

    {
        Eigen::Matrix<double, 6, 6>sqrt_info_pose;
        sqrt_info_pose.setZero();
        if (USE_GNSS && !ESTIMATE_P_WI_WC)
            sqrt_info_pose.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 0;
        else
            sqrt_info_pose.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 1e3;
        if (!ESTIMATE_R_WI_WC)
            sqrt_info_pose.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 0;
        else
            sqrt_info_pose.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 1e3;

        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(
            new InitialPoseFactor(Eigen::Vector3d(para_pose[0]), Quaterniond(Rs[0]*RIC[0]), sqrt_info_pose),
            0, std::vector<double*> {para_pose[0]}, std::vector<int> {}, std::vector<int> {});
        marginalization_info->addResidualBlockInfo(residual_block_info);
        if (USE_GNSS)ASSERT(Eigen::Vector3d(para_pose[0]).norm() < 1e-2);

    }

    {
        Eigen::Matrix<double, 9, 9>sqrt_info_bias;
        sqrt_info_bias.setZero();
        sqrt_info_bias.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 1e1;
        sqrt_info_bias.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 1;
        sqrt_info_bias.block<3, 3>(6, 6) = Eigen::Matrix<double, 3, 3>::Identity() * 1e1;
        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(
            new InitialBiasFactor(Eigen::Vector3d(para_speed_bias[0]), Bas[0], Bgs[0], sqrt_info_bias),
            0, std::vector<double*> {para_speed_bias[0]}, std::vector<int> {}, std::vector<int> {});
        marginalization_info->addResidualBlockInfo(residual_block_info);
    }
    if (ESTIMATE_EXTRINSIC) {
        Eigen::Matrix<double, 6, 6>sqrt_info_pose;
        sqrt_info_pose.setZero();
        sqrt_info_pose.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 100;
        sqrt_info_pose.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 100;
        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(
            new InitialPoseFactor(TIC[0], QIC[0], sqrt_info_pose),
            0, std::vector<double*> {para_extrinsic}, std::vector<int> {}, std::vector<int> {});
        marginalization_info->addResidualBlockInfo(residual_block_info);
    }
// #if ESTIMATE_EXTRINSIC2&&USE_STEREO
//     {
//         Eigen::Matrix<double, 6, 6>sqrt_info_pose;
//         sqrt_info_pose.setZero();
//         sqrt_info_pose.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 1;
//         sqrt_info_pose.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * 1;
//         ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(
//             new InitialPoseFactor(TC2C1, Eigen::Quaterniond(RC2C1), sqrt_info_pose),
//             0, std::vector<double*> {para_extrinsic2}, std::vector<int> {}, std::vector<int> {});
//         marginalization_info->addResidualBlockInfo(residual_block_info);
//     }
// #endif
    if (ESTIMATE_ACC_SCALE) {
        Eigen::Matrix<double, 3, 3>sqrt_info_imu_scale;
        sqrt_info_imu_scale.setZero();
        sqrt_info_imu_scale.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 200;
        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(
            new InitialFactor33(Eigen::Vector3d({1, 1, 1}), sqrt_info_imu_scale),
            0, std::vector<double*> {acc_scale.data()}, std::vector<int> {}, std::vector<int> {});
        marginalization_info->addResidualBlockInfo(residual_block_info);
    }

    if (ESTIMATE_PBG) {
        Eigen::Matrix<double, 3, 3>sqrt_info_pbg;
        sqrt_info_pbg.setZero();
        sqrt_info_pbg.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * 200;
        ResidualBlockInfo* residual_block_info = new ResidualBlockInfo(
            new InitialFactor33(Pbg, sqrt_info_pbg),
            0, std::vector<double*> {Pbg.data()}, std::vector<int> {}, std::vector<int> {});
        marginalization_info->addResidualBlockInfo(residual_block_info);
    }




    marginalization_info->marginalize(true);
    marginalization_info->getParameterBlocks(true);

    if (USE_REPARAM)marginalization_info->UpdateF();


    if (last_marg_info) delete last_marg_info;
    last_marg_info = marginalization_info;

}
