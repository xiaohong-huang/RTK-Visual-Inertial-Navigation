
#pragma once
#include "../parameter/parameters.h"
#include "integration_base.h"
#include "../gnss/common_function.h"
#include "ceres/residual_block.h"
#include "ceres/ceres.h"
#include "../utility/utility.h"
#include "../gnss/gnss_imu_factor.h"

extern std::list<PBtype> rtk_phase_bias_variables[MAXSATNUM * 2];
//storing the rover-only phase biases
extern std::list<PBtype> spp_phase_bias_variables[MAXSATNUM * 2];
//storing the pseudorange correction
extern std::list<PBtype> pseudorange_correction_variables[MAXSATNUM * 2];

class GNSSIMUInfo  {
  public:
    GNSSIMUInfo() = delete;
    GNSSIMUInfo(IntegrationBase* pre_integration_) {
        pre_integrations.push_back(pre_integration_);
        ASSERT(imu_gnss_base == 0);
        imu_gnss_base = new IMUGNSSBase();
        imu_gnss_base->AddImuFactor(pre_integrations[pre_integrations.size() - 1]);
    }

    ~GNSSIMUInfo() {
        for (int i = 0; i < (int)pre_integrations.size(); i++)
            if (pre_integrations[i])delete pre_integrations[i];
        for (int i = 0; i < (int)Rovers.size(); i++) {
            if (Rovers[i])delete Rovers[i];
        }


        for (int i = 0; i < (int)para_pose.size(); i++)
            if (para_pose[i])delete para_pose[i];
        for (int i = 0; i < (int)para_speed_bias.size(); i++)
            if (para_speed_bias[i])delete para_speed_bias[i];

    }


    void AddIMUinfo(IntegrationBase* pre_integration_) {

        pre_integrations.push_back(pre_integration_);
        ASSERT(pre_integrations.size() > 1);

        imu_gnss_base->AddImuFactor(pre_integrations[pre_integrations.size() - 1]);
    }

    void RegenerateGNSSIMUbase() {
        ASSERT(imu_gnss_base);
        delete imu_gnss_base;
        imu_gnss_base = new IMUGNSSBase();
        imu_gnss_base->AddImuFactor(pre_integrations[0]);
        ASSERT(pre_integrations.size() == Rovers.size() + 1);
        for (int i = 0; i < (int)pre_integrations.size() - 1; i++) {
            imu_gnss_base->AddMargInfo(Rovers[i]->marg_info_gnss, para_pose[i], para_speed_bias[i]);
            imu_gnss_base->AddImuFactor(pre_integrations[i + 1]);
        }
    }

    void AddGNSSInfo(mea_t* rover, Eigen::Vector3d Pb, Eigen::Vector3d Vb, Eigen::Matrix3d Rb, Eigen::Vector3d Ba0, Eigen::Vector3d Bg0) {

        Eigen::Vector3d Pg = Pb + Rb * Pbg;
        Eigen::Vector3d Vg = Vb + Rb * (Utility::skewSymmetric(pre_integrations[pre_integrations.size() - 1]->gyrj - Bg0) * Pbg);

        Rovers.push_back(rover);
        ASSERT(pre_integrations.size() == Rovers.size());

        Pgs.push_back(Pg);
        Vgs.push_back(Vg);
        Rs.push_back(Rb);
        Bas.push_back(Ba0);
        Bgs.push_back(Bg0);

        para_pose.push_back(new double [7]);
        para_speed_bias.push_back(new double [9]);
        memset(para_pose[para_pose.size() - 1], 0, sizeof(double) * 7);
        memset(para_speed_bias[para_speed_bias.size() - 1], 0, sizeof(double) * 9);

        GnssProcess(rover);

        if (!imu_gnss_base)imu_gnss_base = new IMUGNSSBase();
        imu_gnss_base->AddMargInfo(Rovers[Rovers.size() - 1]->marg_info_gnss,
                                   para_pose[Rovers.size() - 1], para_speed_bias[Rovers.size() - 1]);
    }


    void resetNoise(const Eigen::Vector3d& Ba0, const Eigen::Vector3d& Bg0, const Eigen::Vector3d& acc_scale, const Eigen::Vector3d& gyr_scale) {
        for (int i = 0; i < (int)pre_integrations.size(); i++) {
            if (pre_integrations[i]) {
                pre_integrations[i]->noise = Eigen::Matrix<double, 18, 18>::Zero();
                pre_integrations[i]->noise.block<3, 3>(0, 0) =  (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
                pre_integrations[i]->noise.block<3, 3>(3, 3) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
                pre_integrations[i]->noise.block<3, 3>(6, 6) =  (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
                pre_integrations[i]->noise.block<3, 3>(9, 9) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
                pre_integrations[i]->noise.block<3, 3>(12, 12) =  (ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
                pre_integrations[i]->noise.block<3, 3>(15, 15) =  (GYR_W * GYR_W) * Eigen::Matrix3d::Identity();
                if (i == 0)pre_integrations[i]->repropagate(Ba0, Bg0, acc_scale, gyr_scale);
                else pre_integrations[i]->repropagate(Bas[i - 1], Bgs[i - 1], acc_scale, gyr_scale);
            }
        }
    }


    void Setdt(double dt_i, double dt_j) {
        pre_integrations[pre_integrations.size() - 1]->dt_j = dt_j;
        pre_integrations[pre_integrations.size() - 1]->dt_i = dt_i;
        if (pre_integrations.size() > 1)ASSERT(pre_integrations[pre_integrations.size() - 1]->dt_i == pre_integrations[pre_integrations.size() - 2]->dt_j);
    }


    void reprogate(const Eigen::Vector3d& Ba0, const Eigen::Vector3d& Bg0, const Eigen::Vector3d& acc_scale, const Eigen::Vector3d& gyr_scale) {
        for (int j = 0; j < (int)pre_integrations.size(); j++) {
            Eigen::Vector3d Ba, Bg;
            if (j == 0) {
                Ba = Ba0;
                Bg = Bg0;
            } else {
                Ba = Bas[j - 1];
                Bg = Bgs[j - 1];
            }
            if (pre_integrations[j]->dt_buf.size() > 5000)continue;
            if ((Ba - pre_integrations[j]->linearized_ba).norm() > 0.01 || (Bg - pre_integrations[j]->linearized_bg).norm() > 0.001 || (acc_scale - pre_integrations[j]->acc_scale).norm() > 0.001 || (gyr_scale - pre_integrations[j]->gyr_scale).norm() > 0.001)
                pre_integrations[j]->repropagate(Ba, Bg, acc_scale, gyr_scale
                                                );
        }
    }


    void merge_gnss_imu_factor(GNSSIMUInfo* nextfactor) {
        IntegrationBase* pre_integration_tail = pre_integrations[pre_integrations.size() - 1];
        auto& next_dt_buf = nextfactor->pre_integrations[0]->dt_buf;
        auto& next_acc_buf = nextfactor->pre_integrations[0]->acc_buf;
        auto& next_gyr_buf = nextfactor->pre_integrations[0]->gyr_buf;

        for (int i = 0; i < (int)next_dt_buf.size(); i++)
            pre_integration_tail->push_back(next_dt_buf[i], next_acc_buf[i], next_gyr_buf[i]);
        pre_integration_tail->dt_j = nextfactor->pre_integrations[0]->dt_j;
        delete nextfactor->pre_integrations[0];
        ASSERT(nextfactor->pre_integrations.size() == nextfactor->Rovers.size() + 1);
        imu_gnss_base->history_flag = false;
        for (int i = 1; i < (int)nextfactor->pre_integrations.size(); i++) {
            ASSERT(nextfactor->pre_integrations[i]);
            Rovers.push_back(nextfactor->Rovers[i - 1]);
            ASSERT(pre_integrations.size() == Rovers.size());

            Pgs.push_back(nextfactor->Pgs[i - 1]);
            Vgs.push_back(nextfactor->Vgs[i - 1]);
            Rs.push_back(nextfactor->Rs[i - 1]);
            Bas.push_back(nextfactor->Bas[i - 1]);
            Bgs.push_back(nextfactor->Bgs[i - 1]);
            para_pose.push_back(nextfactor->para_pose[i - 1]);
            para_speed_bias.push_back(nextfactor->para_speed_bias[i - 1]);
            nextfactor->para_pose[i - 1] = 0;
            nextfactor->para_speed_bias[i - 1] = 0;
            nextfactor->Rovers[i - 1] = 0;

            if (!imu_gnss_base)imu_gnss_base = new IMUGNSSBase();
            imu_gnss_base->AddMargInfo(Rovers[Rovers.size() - 1]->marg_info_gnss,
                                       para_pose[Rovers.size() - 1], para_speed_bias[Rovers.size() - 1]);

            AddIMUinfo(nextfactor->pre_integrations[i]);
            nextfactor->pre_integrations[0] = 0;
        }
        imu_gnss_base->history_flag = false;

        nextfactor->pre_integrations.clear();

    }


    void GnssProcess(mea_t* rover);
    void AddCeresGnssResidual(std::set<double*>MargePoint, MarginalizationInfo* marginalization_info,
                          ceres::Problem& problem,
                          ceres::Solver::Options& options, int ir);
    void Vector2Double();
    void Double2Vector();
    void GnssPreprocess(mea_t* data);
    void PhaseBiasSaveAndReset();
    void PhaseBiasRestore();
    void AddGnssResidual(int mode, std::set<double*>MargePoint, MarginalizationInfo* marginalization_info,int ir);
    void CycleSlipDetection(mea_t* data);
    void GenerateGNSSPriorFactor(mea_t* rove);

    std::vector<Eigen::Vector3d>        Bas;
    std::vector<Eigen::Vector3d>        Bgs;
    std::vector<Eigen::Vector3d>        Pgs;
    std::vector<Eigen::Vector3d>        Vgs;
    std::vector<Eigen::Matrix3d>        Rs;



    std::vector<IntegrationBase*> pre_integrations;
    std::vector<mea_t*>Rovers;
    std::vector<double*> para_pose;
    std::vector<double*> para_speed_bias;
    std::vector<double>phase_bias_save;
    IMUGNSSBase* imu_gnss_base = 0;
    void LambdaSearch(std::vector<mea_t*>& rovers, std::vector<double*>& parameter_block_addr_,
                      std::vector<int>& parameter_block_global_size_, Eigen::MatrixXd A_, Eigen::VectorXd& b_);
};
