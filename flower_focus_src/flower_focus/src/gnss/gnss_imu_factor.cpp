#include "gnss_imu_factor.h"
#include "../factor/operation.h"
#include "../utility/utility.h"
#include "../parameter/parameters.h"

const int mapindex[O_Size - 1] = {O_Pose0, O_Pose1, O_EXTRINSIC, O_N};

#define MatrixX15d_ Eigen::Matrix<double,Eigen::Dynamic,15,Eigen::RowMajor>



IMUGNSSBase::~IMUGNSSBase() {
    Clear();
    for (int i = 0; i < (int)imu_factors.size(); i++)
        delete imu_factors[i];
    if (gnss_middle_marginfo)delete gnss_middle_marginfo;
}

void IMUGNSSBase::Init() {

    InitHessianRhs(hessian55, hessian_size, rhs5);
    phase_biases_values.resize(hessian_size[O_N]);
    history_flag = false;

    hessian_index[0] = 0;
    global_size = 0;
    for (int i = 0; i < O_Size - 1 - 1; i++)
        hessian_index[i + 1] = hessian_index[i] + hessian_size[mapindex[i]];
    for (int i = 0; i < O_Size - 1; i++)
        global_size += hessian_size[mapindex[i]];


    order2p_local.clear();
    int p_local = 0;
    int p_idx = 0;
    if (ESTIMATE_SCALE) {
        O_SCALE = p_idx; p_idx++; order2p_local.push_back(p_local); p_local += 1;
    }
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC) {
        O_POSGLOBAL = p_idx; p_idx++; order2p_local.push_back(p_local); p_local += 6;
    }
    if (ESTIMATE_EXTRINSIC) {
        O_EXTRINSIC_IC = p_idx; p_idx++;
        order2p_local.push_back(p_local); p_local += 6;
    }

    O_POSi = p_idx; order2p_local.push_back(p_local); p_local += 6; p_idx++;
    O_BIASi = p_idx; order2p_local.push_back(p_local); p_local += 9; p_idx++;
    O_POSj = p_idx; order2p_local.push_back(p_local); p_local += 6; p_idx++;
    O_BIASj = p_idx; order2p_local.push_back(p_local); p_local += 9; p_idx++;

    if (ESTIMATE_TD) {
        O_TD = p_idx; p_idx++;
        order2p_local.push_back(p_local); p_local += 1;
    }
    if (ESTIMATE_ACC_SCALE) {
        O_ACC_S = p_idx; p_idx++;
        order2p_local.push_back(p_local); p_local += 3;
    }
    if (ESTIMATE_PBG) {
        O_PBG = p_idx; p_idx++;
        order2p_local.push_back(p_local); p_local += 3;
    }

    if (gnss_phase_biases.size()) {
        O_PHASEBIAS = p_idx; p_idx ++;
        order2p_local.push_back(p_local); p_local += gnss_phase_biases.size();
    }


    O_FULL = p_idx; p_idx++; order2p_local.push_back(p_local);

    int extrinsic_globalsize = 0;
    if (ESTIMATE_TD)extrinsic_globalsize += 1;
    if (ESTIMATE_ACC_SCALE)extrinsic_globalsize += 3;
    if (ESTIMATE_PBG)extrinsic_globalsize += 3;
    if (extrinsic_globalsize == 0)extrinsic_globalsize = 1;

    extrinsic = Eigen::VectorXd(extrinsic_globalsize);
    phase_biases_values.setZero();

    for (int i = 0; i < (int)imu_factors.size(); i++) {

        if (i == 0) {
            if (imu_factors.size() == 1)
                imu_factors[i]->parameters = std::vector<double*>({PgQbi.data(), Bi.data(), PgQbj.data(), Bj.data(), extrinsic.data()});
            else
                imu_factors[i]->parameters = std::vector<double*>({PgQbi.data(), Bi.data(), gnss_poses[0], gnss_speed_bias[0], extrinsic.data()});
        } else if (i == (int)imu_factors.size() - 1)imu_factors[i]->parameters = std::vector<double*>({gnss_poses[i - 1], gnss_speed_bias[i - 1], PgQbj.data(), Bj.data(), extrinsic.data()});
        else imu_factors[i]->parameters = std::vector<double*>({gnss_poses[i - 1], gnss_speed_bias[i - 1], gnss_poses[i], gnss_speed_bias[i], extrinsic.data()});


        imu_factors[i]->imu_jacobians_in.resize(3);
        imu_factors[i]->imu_jacobians_in[0] = flower_focus::Matrix15d();
        imu_factors[i]->imu_jacobians_in[1] = flower_focus::Matrix15d();
        imu_factors[i]->imu_jacobians_in[2] = flower_focus::MatrixXd(15, hessian_size[O_EXTRINSIC]);
        imu_factors[i]->imu_jacobians_raw_in[0] = imu_factors[i]->imu_jacobians_in[0].data();
        imu_factors[i]->imu_jacobians_raw_in[1] = imu_factors[i]->imu_jacobians_in[1].data();
        imu_factors[i]->imu_jacobians_raw_in[2] = imu_factors[i]->imu_jacobians_in[2].data();
        imu_factors[i]->model_residual_accum.setZero();

    }


}

void IMUGNSSBase::Clear() {
    for (int i = 0; i < O_Size; i++) {
        if (rhs5[i])delete rhs5[i];
        if (gradient5[i])delete gradient5[i];
        if (inc5[i])delete inc5[i];

        if (gauss5[i])delete gauss5[i];
        for (int j = i; j < O_Size; j++) {
            if (hessian55[i * O_Size + j])delete hessian55[i * O_Size + j];
        }
    }
    for (int i = 0; i < O_Size; i++) {
        for (int j = 0; j < (int)hmn_save[i].size(); j++) {
            if (hmn_save[i][j])delete hmn_save[i][j];
        }
    }
    for (int i = 0; i < (int)rhsmn_save.size(); i++) {
        if (rhsmn_save[i])delete rhsmn_save[i];
        if (gradient_save[i])delete gradient_save[i];
        if (gauss_save[i])delete gauss_save[i];
        if (inc_save[i])delete inc_save[i];

    }



}

void IMUGNSSBase::InitHessianRhs(double** hessian55, int* hessian_size, double** rhs5) {
    int extrinsic_localsize = 0;
    if (ESTIMATE_TD)extrinsic_localsize += 1;
    if (ESTIMATE_ACC_SCALE)extrinsic_localsize += 3;
    if (ESTIMATE_PBG)extrinsic_localsize += 3;
    if (extrinsic_localsize == 0)extrinsic_localsize = 1;

    hessian_size[O_Pose1] = 15;
    hessian_size[O_Pose2] = 15;
    hessian_size[O_N] = (int)gnss_phase_biases.size();
    hessian_size[O_Pose0] = 15;
    hessian_size[O_EXTRINSIC] = extrinsic_localsize;
    // if (hessian_size[O_N] == 0)hessian_size[O_N] = 1;


    for (int i = 0; i < O_Size; i++) {
        rhs5[i] = new double[hessian_size[i]];
        gradient5[i] = new double[hessian_size[i]];
        inc5[i] = new double[hessian_size[i]];

        gauss5[i] = new double[hessian_size[i]];
        for (int j = i; j < O_Size; j++) {
            int size = hessian_size[i] * hessian_size[j];
            hessian55[i * O_Size + j] = new double[size];
        }
    }
    rhsmn_save.resize(gnss_speed_bias.size());
    for (int i = 0; i < (int)gnss_speed_bias.size(); i++)
        rhsmn_save[i] = new double [hessian_size[O_Pose1]];

    gradient_save.resize(gnss_speed_bias.size());
    for (int i = 0; i < (int)gnss_speed_bias.size(); i++)
        gradient_save[i] = new double [hessian_size[O_Pose1]];

    gauss_save.resize(gnss_speed_bias.size());
    for (int i = 0; i < (int)gnss_speed_bias.size(); i++)
        gauss_save[i] = new double [hessian_size[O_Pose1]];

    inc_save.resize(gnss_speed_bias.size());
    for (int i = 0; i < (int)gnss_speed_bias.size(); i++)
        inc_save[i] = new double [hessian_size[O_Pose1]];


    for (int i = O_Pose1; i < O_Size; i++) {
        hmn_save[i].resize(gnss_speed_bias.size());
        for (int j = 0; j < (int)gnss_speed_bias.size(); j++)
            hmn_save[i][j] = new double[ hessian_size[O_Pose1] * hessian_size[i] ];
    }
}

void IMUGNSSBase::AddImuFactor(IntegrationBase* pre_integration) {
    imu_factors.push_back(new IMUPreFactor(new IMUFactor(pre_integration)));
    Clear();
    Init();
    history_flag = false;
}

void IMUGNSSBase::AddMidMargInfo(MarginalizationInfo* gnss_middle_marginfo_) {

    // ASSERT(gnss_middle_marginfo == 0);
    // gnss_middle_marginfo = gnss_middle_marginfo_;


    // const flower_focus::MatrixXd& GNSS_IMU_Hessian = gnss_middle_marginfo->A;
    // const flower_focus::VectorXd& GNSS_IMU_rhs = gnss_middle_marginfo->b;
    // gnss_init_cost += gnss_middle_marginfo->init_cost;

    // int set_num = 0;

    // for (int i = 0; i < (int)imu_factors.size(); i++) {
    //     IntegrationBase* pre_integration = imu_factors[i]->pre_integration;
    //     if (pre_integration->sum_dt <= 100000)continue;
    //     set_num++;
    //     gnss_Index = i;
    // }
    // ASSERT(set_num == 1);
    // int N_size_external = 0;
    // int oldNsize = (int)gnss_phase_biases.size();
    // std::map<double*, int>parameter_index;
    // for (int j = 0; j < (int)gnss_middle_marginfo->keep_block_addr.size(); j++) {
    //     double* pointer = gnss_middle_marginfo->keep_block_addr[j];
    //     int size = gnss_middle_marginfo->keep_block_size[j];
    //     if (size == 7)
    //         ASSERT(pointer == gnss_poses[gnss_Index] || pointer == gnss_poses[gnss_Index - 1]);
    //     else if (size == 9)
    //         ASSERT(pointer == gnss_speed_bias[gnss_Index] || pointer == gnss_speed_bias[gnss_Index - 1]);
    //     else if (size == 1) {
    //         if (gnss_phase_biases_set.find(pointer) == gnss_phase_biases_set.end()) {
    //             gnss_phase_biases_set.insert(pointer);
    //             gnss_phase_biases.push_back(pointer);
    //             N_size_external++;
    //         }
    //     } else
    //         ASSERT(0);
    //     parameter_index[pointer] = gnss_middle_marginfo->keep_block_idx[j] - gnss_middle_marginfo->m;
    // }
    // int newNsize = (int)gnss_phase_biases.size();
    // ASSERT(oldNsize + N_size_external == newNsize);
    // //尺寸扩展
    // if (N_size_external) {
    //     phase_biases_hessians.conservativeResize(newNsize, newNsize);
    //     phase_biases_hessians.block(oldNsize, oldNsize, N_size_external, N_size_external).setZero();
    //     phase_biases_hessians.block(0, oldNsize, oldNsize, N_size_external).setZero();
    //     phase_biases_hessians.block(oldNsize, 0, N_size_external, oldNsize).setZero();
    //     phase_biases_rhs.conservativeResize(newNsize);
    //     phase_biases_rhs.segment(oldNsize, N_size_external).setZero();
    //     for (int i = 0; i < (int)pose_phase_biases_hessians.size(); i++) {
    //         pose_phase_biases_hessians[i].conservativeResize(15, newNsize);
    //         pose_phase_biases_hessians[i].block(0, oldNsize, 15, N_size_external).setZero();
    //     }
    // }


    // std::vector<double*>hessian55_m_pointer({gnss_poses[gnss_Index - 1], gnss_speed_bias[gnss_Index - 1], gnss_poses[gnss_Index], gnss_speed_bias[gnss_Index]});
    // std::vector<int>hessian55_m_index({O_Pose1, O_Pose1, O_Pose2, O_Pose2});
    // std::vector<int>hessian55_m_size({6, 9, 6, 9});
    // std::vector<int>hessian55_m_shift({0, 6, 0, 6});

    // for (int i = 0; i < (int)gnss_phase_biases.size(); i++) {
    //     hessian55_m_pointer.push_back(gnss_phase_biases[i]);
    //     hessian55_m_index.push_back(O_N);
    //     hessian55_m_size.push_back(1);
    //     hessian55_m_shift.push_back(i);
    // }

    // for (int i = 0; i < (int)hessian55_m_pointer.size(); i++) {
    //     double* pointer1 = hessian55_m_pointer[i];
    //     if (parameter_index.find(pointer1) == parameter_index.end())continue;
    //     int index1 = parameter_index[pointer1];
    //     int shift1 = hessian55_m_shift[i];
    //     int size1 = hessian55_m_size[i];
    //     int indexm1 = hessian55_m_index[i];
    //     if (indexm1 == O_N)
    //         phase_biases_rhs(shift1) += GNSS_IMU_rhs(index1);
    //     else if (indexm1 == O_Pose1)
    //         pose_rhses[gnss_Index - 1].segment(shift1, size1) += GNSS_IMU_rhs.segment(index1, size1);
    //     else if (indexm1 == O_Pose2)
    //         pose_rhses[gnss_Index].segment(shift1, size1) += GNSS_IMU_rhs.segment(index1, size1);
    //     else
    //         ASSERT(0);
    //     for (int j = i; j < (int)hessian55_m_pointer.size(); j++) {
    //         double* pointer2 = hessian55_m_pointer[j];
    //         if (parameter_index.find(pointer2) == parameter_index.end())continue;
    //         int index2 = parameter_index[pointer2];
    //         int shift2 = hessian55_m_shift[j];
    //         int size2 = hessian55_m_size[j];
    //         int indexm2 = hessian55_m_index[j];
    //         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix = GNSS_IMU_Hessian.block(index1, index2, size1, size2);
    //         if (indexm1 == O_N && indexm2 == O_N)
    //             phase_biases_hessians.block(shift1, shift2, size1, size2) += matrix;
    //         else if (indexm1 == O_Pose1 && indexm2 == O_N)
    //             pose_phase_biases_hessians[gnss_Index - 1].block(shift1, shift2, size1, size2) += matrix;
    //         else if (indexm1 == O_Pose2 && indexm2 == O_N)
    //             pose_phase_biases_hessians[gnss_Index].block(shift1, shift2, size1, size2) += matrix;
    //         else if (indexm1 == O_Pose1 && indexm2 == O_Pose1)
    //             pose_hessians[gnss_Index - 1].block(shift1, shift2, size1, size2) += matrix;
    //         else if (indexm1 == O_Pose2 && indexm2 == O_Pose2)
    //             pose_hessians[gnss_Index].block(shift1, shift2, size1, size2) += matrix;
    //         else if (indexm1 == O_Pose1 && indexm2 == O_Pose2)
    //             pose1_pose2_hessians.block(shift1, shift2, size1, size2) = matrix;
    //         else
    //             ASSERT(0);

    //     }
    // }

    // phase_biases_hessians = phase_biases_hessians.selfadjointView<Eigen::Upper>();
    // // pose1_pose2_hessians=pose1_pose2_hessians.selfadjointView<Eigen::Upper>();
    // // std::cout<<pose1_pose2_hessians<<std::endl;
    // pose_hessians[gnss_Index - 1] = pose_hessians[gnss_Index - 1].selfadjointView<Eigen::Upper>();
    // pose_hessians[gnss_Index] = pose_hessians[gnss_Index].selfadjointView<Eigen::Upper>();
    // Init();




}


void IMUGNSSBase::AddMargInfo(MarginalizationInfo* marg_info_gnss, double* para_Poses, double* para_speed_bias) {
    Clear();
    history_flag = false;
    int N_size_external = 0;
    std::map<double*, int>parameter_index;
    gnss_poses.push_back(para_Poses);
    gnss_speed_bias.push_back(para_speed_bias);
    gnss_speed_bias_lin.push_back(para_speed_bias);
    gnss_poses_lin.push_back(para_Poses);

    gnss_poses_old.resize(gnss_poses.size());
    gnss_speed_bias_old.resize(gnss_speed_bias.size());

    int oldNsize = (int)gnss_phase_biases.size();
    for (int i = 0; i < (int)marg_info_gnss->keep_block_addr.size(); i++) {
        double* pointer = marg_info_gnss->keep_block_addr[i];
        int size = marg_info_gnss->keep_block_size[i];
        if (size == 7)size = 6;
        if (size == 6) {
            gnss_poses_lin[gnss_poses_lin.size() - 1] = marg_info_gnss->keep_block_data[i];
            gnss_poses[gnss_poses.size() - 1] = marg_info_gnss->keep_block_addr[i];
        } else if (size == 9) {
            gnss_speed_bias_lin[gnss_speed_bias_lin.size() - 1] = marg_info_gnss->keep_block_data[i];
            gnss_speed_bias[gnss_speed_bias.size() - 1] = marg_info_gnss->keep_block_addr[i];
        }
        parameter_index[pointer] = marg_info_gnss->keep_block_idx[i] - marg_info_gnss->m;
        if (size == 1) {
            ASSERT(fabs(marg_info_gnss->keep_block_data[i][0]) == 0);
            if (gnss_phase_biases_set.find(pointer) == gnss_phase_biases_set.end()) {
                gnss_phase_biases_set.insert(pointer);
                gnss_phase_biases.push_back(pointer);
                N_size_external++;
            }
        }
    }
    int newNsize = (int)gnss_phase_biases.size();
    ASSERT(oldNsize + N_size_external == newNsize);
    //尺寸扩展
    if (N_size_external) {
        phase_biases_hessians.conservativeResize(newNsize, newNsize);
        phase_biases_hessians.block(oldNsize, oldNsize, N_size_external, N_size_external).setZero();
        phase_biases_hessians.block(0, oldNsize, oldNsize, N_size_external).setZero();
        phase_biases_hessians.block(oldNsize, 0, N_size_external, oldNsize).setZero();
        phase_biases_rhs.conservativeResize(newNsize);
        phase_biases_rhs.segment(oldNsize, N_size_external).setZero();
        for (int i = 0; i < (int)pose_phase_biases_hessians.size(); i++) {
            pose_phase_biases_hessians[i].conservativeResize(15, newNsize);
            pose_phase_biases_hessians[i].block(0, oldNsize, 15, N_size_external).setZero();
        }
    }

    flower_focus::Matrix15d Pose_Hessian;
    flower_focus::MatrixXd Pose_N_Hessian(15, newNsize);
    flower_focus::Vector15d Pose_rhs;

    Pose_Hessian.setZero();
    Pose_N_Hessian.setZero();
    Pose_rhs.setZero();

    const flower_focus::MatrixXd& GNSS_Hessian = marg_info_gnss->A;
    const flower_focus::VectorXd& GNSS_rhs = marg_info_gnss->b;
    gnss_init_cost += marg_info_gnss->init_cost;

    //更新hessian和rhs
    for (int i = 0; i < newNsize; i++) {
        double* pointer1 = gnss_phase_biases[i];
        if (parameter_index.find(pointer1) == parameter_index.end())continue;
        int index1 = parameter_index[pointer1];
        phase_biases_rhs(i) += GNSS_rhs(index1);
        for (int j = 0; j < newNsize; j++) {
            double* pointer2 = gnss_phase_biases[j];
            if (parameter_index.find(pointer2) == parameter_index.end())continue;
            int index2 = parameter_index[pointer2];
            phase_biases_hessians(i, j) += GNSS_Hessian(index1, index2);
        }
    }

    double* Pose = gnss_poses[gnss_poses.size() - 1];
    double* SpeedBias = gnss_speed_bias[gnss_speed_bias.size() - 1];
    for (int i = 0; i < 2; i++) {
        double* pointer1 = i == 0 ? Pose : SpeedBias;
        if (parameter_index.find(pointer1) == parameter_index.end())continue;
        int index1 = parameter_index[pointer1];
        int index1_ = i == 0 ? 0 : 6;
        int size1 = i == 0 ? 6 : 9;
        Pose_rhs.segment(index1_, size1) += GNSS_rhs.segment(index1, size1);
        for (int j = 0; j < 2; j++) {
            double* pointer2 = j == 0 ? Pose : SpeedBias;
            if (parameter_index.find(pointer2) == parameter_index.end())continue;
            int index2 = parameter_index[pointer2];
            int index2_ = j == 0 ? 0 : 6;
            int size2 = j == 0 ? 6 : 9;
            Pose_Hessian.block(index1_, index2_, size1, size2) += GNSS_Hessian.block(index1, index2, size1, size2); //(index2>index1?GNSS_Hessian.block(index1,index2,size1,size2):GNSS_Hessian.block(index2,index1,size1,size2));
        }
        for (int j = 0; j < newNsize; j++) {
            double* pointer2 = gnss_phase_biases[j];
            if (parameter_index.find(pointer2) == parameter_index.end())continue;
            int index2 = parameter_index[pointer2];
            int index2_ = j;
            Pose_N_Hessian.block(index1_, index2_, size1, 1) += GNSS_Hessian.block(index1, index2, size1, 1); //(index2>index1?GNSS_Hessian.block(index1,index2,size1,1):GNSS_Hessian.block(index2,index1,size1,1));
        }
    }


    pose_hessians.push_back(Pose_Hessian);
    pose_phase_biases_hessians.push_back(Pose_N_Hessian);
    pose_rhses.push_back(Pose_rhs);

    Init();
}

void IMUGNSSBase::JacobianResidualUpdateHessianRhs(std::vector<int>index, std::vector<int>mapsize, double** jacobian, double** hessian55, double* residual, double** rhs5, double** gradient5) {
    for (int i = 0; i < (int)index.size(); i++) {
        int sizei = mapsize[i];
        int indexi = index[i];
        Eigen::VectorXd tmp(sizei);
        MatrixTransposeVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 0>
        (jacobian[i], 15, sizei, residual, tmp.data());
        Eigen::Map<Eigen::VectorXd>(rhs5[indexi], sizei) += tmp;
        Eigen::Map<Eigen::VectorXd>(gradient5[indexi], sizei) += tmp;
        if (hessian55) {
            for (int j = 0; j < (int)index.size(); j++) {
                int sizej = mapsize[j];
                int indexj = index[j];
                double* hessian_local = hessian55[indexi * O_Size + indexj];
                if (indexj < indexi)continue;
                MatrixTransposeMatrixMultiply<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, 1>
                (jacobian[i], 15, sizei, jacobian[j], 15, sizej, hessian_local, 0, 0, sizei, sizej);

            }
        }
    }
}

void IMUGNSSBase::HessianUpdateHessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2) {
    ASSERT(block2 >= block1);
    if (size1 == 15 && size2 == 15) {
        Eigen::Map<flower_focus::Matrix15d>(hessian_desc[block1 * O_Size + block2]) +=
            Eigen::Map<flower_focus::Matrix15d>(hessian_src);
    } else {
        Eigen::Map<flower_focus::MatrixXd>(hessian_desc[block1 * O_Size + block2], size1, size2) +=
            Eigen::Map<flower_focus::MatrixXd>(hessian_src, size1, size2);
    }
}

void IMUGNSSBase::RhsUpdateRhs(double** rhs_desc, double* rhs_src, int block1, int size1) {
    if (size1 == 15)
        Eigen::Map<flower_focus::Vector15d>(rhs_desc[block1]) += Eigen::Map<flower_focus::Vector15d>(rhs_src);
    else
        Eigen::Map<flower_focus::VectorXd >(rhs_desc[block1], size1) += Eigen::Map<flower_focus::VectorXd > (rhs_src, size1);
}

void IMUGNSSBase::CopyHessian2Hessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2) {
    ASSERT(block2 >= block1);
    memcpy(hessian_desc[block1 * O_Size + block2], hessian_src, sizeof(double)*size1 * size2);
}

void IMUGNSSBase::CopyRhs2Rhs(double** rhs_desc, double* rhs_src, int block1, int size1) {
    memcpy(rhs_desc[block1], rhs_src, sizeof(double)*size1);
}

void IMUGNSSBase::MargPose1(double** hessian55, int* hessian_size, double** rhs5) {

    flower_focus::Matrix15d Amm_inverse = InvertPSDMatrix(Eigen::Map<flower_focus::Matrix15d>(hessian55[O_Pose1 * O_Size + O_Pose1]) );

    memcpy(hessian55[O_Pose1 * O_Size + O_Pose1], Amm_inverse.data(), sizeof(double) * 15 * 15);

    for (int i = O_Pose1 + 1; i < O_Size; i++) {
        int sizen = hessian_size[i];
        if (sizen == 0)continue;
        flower_focus::MatrixXd Anm_Amm_inverse(sizen, 15);
        MatrixTransposeMatrixMultiply<15, Eigen::Dynamic, 15, 15, 0>(
            hessian55[O_Pose1 * O_Size + i], 15, sizen,
            Amm_inverse.data(), 15, 15,
            Anm_Amm_inverse.data(), 0, 0, sizen, 15);

        MatrixVectorMultiply < Eigen::Dynamic, 15, -1 > (
            Anm_Amm_inverse.data(), sizen, 15,
            rhs5[O_Pose1],
            rhs5[i]);

        for (int j = i; j < O_Size; j++) {
            int sizev = hessian_size[j];
            if (sizev == 0)continue;
            MatrixMatrixMultiply
            < Eigen::Dynamic, 15, 15, Eigen::Dynamic, -1 > (
                Anm_Amm_inverse.data(), sizen, 15,
                hessian55[O_Pose1 * O_Size + j], 15, sizev,
                hessian55[i * O_Size + j], 0, 0, sizen, sizev);
        }
        Eigen::Map<flower_focus::MatrixXd>(hessian55[O_Pose1 * O_Size + i], 15, sizen) = Anm_Amm_inverse.transpose();
    }

}



void IMUGNSSBase::MargPose1rhsOnly(double** hessian55, int* hessian_size, double** rhs5, int index) {


    for (int i = O_Pose1 + 1; i < O_Size; i++) {
        int sizen = hessian_size[i];
        if (sizen == 0)continue;
        MatrixTransposeVectorMultiply < Eigen::Dynamic, Eigen::Dynamic, -1 > (
            hmn_save[i][index], 15, sizen,
            rhs5[O_Pose1],
            rhs5[i]);
    }

}

void IMUGNSSBase::MoveHessianData(int index) {


    for (int i = O_Pose1; i < O_Size; i++)
        memcpy(hmn_save[i][index], hessian55[O_Pose1 * O_Size + i], sizeof(double)*hessian_size[O_Pose1]*hessian_size[i]);
    memcpy(hessian55[O_Pose1 * O_Size + O_Pose1], hessian55[O_Pose2 * O_Size + O_Pose2], sizeof(double)*hessian_size[O_Pose2]*hessian_size[O_Pose2]);
    memset(hessian55[O_Pose2 * O_Size + O_Pose2], 0, sizeof(double)*hessian_size[O_Pose2]*hessian_size[O_Pose2]);

    for (int i = O_Pose2 + 1; i < O_Size; i++) {
        memcpy(hessian55[O_Pose1 * O_Size + i], hessian55[O_Pose2 * O_Size + i], sizeof(double)*hessian_size[O_Pose2]*hessian_size[i]);
        memset(hessian55[O_Pose2 * O_Size + i], 0, sizeof(double)*hessian_size[O_Pose2]*hessian_size[i]);
    }
    memset(hessian55[O_Pose1 * O_Size + O_Pose2], 0, sizeof(double)*hessian_size[O_Pose1]*hessian_size[O_Pose2]);

}

void IMUGNSSBase::MoveRhsData(int index) {

    memcpy(rhsmn_save[index], rhs5[O_Pose1], sizeof(double)*hessian_size[O_Pose1]);
    memcpy(gradient_save[index], gradient5[O_Pose1], sizeof(double)*hessian_size[O_Pose1]);


    memcpy(rhs5[O_Pose1], rhs5[O_Pose2], sizeof(double)*hessian_size[O_Pose2]);
    memcpy(gradient5[O_Pose1], gradient5[O_Pose2], sizeof(double)*hessian_size[O_Pose2]);
    memset(rhs5[O_Pose2], 0, sizeof(double)*hessian_size[O_Pose2]);
    memset(gradient5[O_Pose2], 0, sizeof(double)*hessian_size[O_Pose2]);


}

void IMUGNSSBase::BuildMatrixInfo() {






    matrix_info.clear();
    {
        flower_focus::Matrix6d dP0_dP0 = flower_focus::Matrix6d::Zero();
        flower_focus::Matrix6d dP1_dP1 = flower_focus::Matrix6d::Zero();

        dP0_dP0.block<3, 3>(0, 0) = r_wi_wc * scale_factor;
        dP0_dP0.block<3, 3>(0, 3) = r_wi_wc * Rci * Utility::skewSymmetric(ric.transpose() * (tic - pbg));
        dP0_dP0.block<3, 3>(3, 3) = ric;

        dP1_dP1.block<3, 3>(0, 0) = r_wi_wc * scale_factor;
        dP1_dP1.block<3, 3>(0, 3) = r_wi_wc * Rcj * Utility::skewSymmetric(ric.transpose() * (tic - pbg));
        dP1_dP1.block<3, 3>(3, 3) = ric;

        matrix_info.push_back(MatirxInfo(order2p_local[O_POSi], order2p_local[O_POSi], dP0_dP0 - flower_focus::Matrix6d::Identity(), false, 6, 6));
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSj], order2p_local[O_POSj], dP1_dP1 - flower_focus::Matrix6d::Identity(), false, 6, 6));
    }
    if (ESTIMATE_P_WI_WC || ESTIMATE_R_WI_WC) {
        flower_focus::Matrix6d dP0_dR_WI_WC = flower_focus::Matrix6d::Zero();
        flower_focus::Matrix6d dP1_dR_WI_WC = flower_focus::Matrix6d::Zero();
        if (ESTIMATE_P_WI_WC) {
            dP0_dR_WI_WC.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();//TODO:Matrix3d
            dP1_dR_WI_WC.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
        }
        if (ESTIMATE_R_WI_WC) {
            dP0_dR_WI_WC.block<3, 3>(0, 3) = r_wi_wc * -Utility::skewSymmetric(Pci * scale_factor - Rci * ric.transpose() * (tic - pbg));
            dP0_dR_WI_WC.block<3, 3>(3, 3) = Rbi.transpose() * r_wi_wc;//TODO:Matrix6_3d

            dP1_dR_WI_WC.block<3, 3>(0, 3) = r_wi_wc * -Utility::skewSymmetric(Pcj * scale_factor - Rcj * ric.transpose() * (tic - pbg));
            dP1_dR_WI_WC.block<3, 3>(3, 3) = Rbj.transpose() * r_wi_wc;
        }
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSi], order2p_local[O_POSGLOBAL], dP0_dR_WI_WC, false, 6, 6));
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSj], order2p_local[O_POSGLOBAL], dP1_dR_WI_WC, false, 6, 6));
    }

    if (ESTIMATE_SCALE && fix_scale) {
        flower_focus::Vector6d dP0_dscale = flower_focus::Vector6d::Zero();//TODO:Vector3d
        flower_focus::Vector6d dP1_dscale = flower_focus::Vector6d::Zero();

        dP0_dscale.segment(0, 3) = r_wi_wc * Pci;
        dP1_dscale.segment(0, 3) = r_wi_wc * Pcj;

        matrix_info.push_back(MatirxInfo(order2p_local[O_POSi], order2p_local[O_SCALE], dP0_dscale, false, 6, 1));
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSj], order2p_local[O_SCALE], dP1_dscale, false, 6, 1));
    }

    if (ESTIMATE_EXTRINSIC) {

        flower_focus::Matrix6d dP0_dextrinsic = flower_focus::Matrix6d::Zero();
        flower_focus::Matrix6d dP1_dextrinsic = flower_focus::Matrix6d::Zero();

        dP0_dextrinsic.block<3, 3>(0, 0) = - r_wi_wc * Rci * ric.transpose();
        dP0_dextrinsic.block<3, 3>(0, 3) = - r_wi_wc * Rci * Utility::skewSymmetric(ric.transpose() * (tic - pbg));
        dP0_dextrinsic.block<3, 3>(3, 3) = -ric;

        dP1_dextrinsic.block<3, 3>(0, 0) = - r_wi_wc * Rcj * ric.transpose();
        dP1_dextrinsic.block<3, 3>(0, 3) = - r_wi_wc * Rcj * Utility::skewSymmetric(ric.transpose() * (tic - pbg));
        dP1_dextrinsic.block<3, 3>(3, 3) = -ric;

        matrix_info.push_back(MatirxInfo(order2p_local[O_POSi], order2p_local[O_EXTRINSIC_IC], dP0_dextrinsic, false, 6, 6));
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSj], order2p_local[O_EXTRINSIC_IC], dP1_dextrinsic, false, 6, 6));
    }

    if (ESTIMATE_PBG) {

        flower_focus::Matrix6_3d dP0_dpbg = flower_focus::Matrix6_3d::Zero();//TODO:Matrix3_3d
        flower_focus::Matrix6_3d dP1_dpbg = flower_focus::Matrix6_3d::Zero();

        dP0_dpbg.block<3, 3>(0, 0) = Rbi;
        dP1_dpbg.block<3, 3>(0, 0) = Rbj;
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSi], order2p_local[O_PBG], dP0_dpbg, false, 6, 3));
        matrix_info.push_back(MatirxInfo(order2p_local[O_POSj], order2p_local[O_PBG], dP1_dpbg, false, 6, 3));
    }


}

void IMUGNSSBase::UpdateSchurComponent(double* rhs, double* lhs, double* gradient) {


    if (lhs)
        hessian_dense = flower_focus::MatrixXd (global_size, global_size);
    rhs_dense = flower_focus::VectorXd(global_size);
    gradient_dense = flower_focus::VectorXd(global_size);

    for (int i = 0; i < O_Size - 1; i++) {
        int i2 = mapindex[i];
        rhs_dense.segment(hessian_index[i], hessian_size[i2]) = Eigen::Map<flower_focus::VectorXd >(rhs5[i2], hessian_size[i2]);
        gradient_dense.segment(hessian_index[i], hessian_size[i2]) = Eigen::Map<flower_focus::VectorXd >(gradient5[i2], hessian_size[i2]);
        if (lhs) {
            for (int j = i; j < O_Size - 1; j++) {
                int j2 = mapindex[j];
                if (j2 >= i2) {
                    hessian_dense.block(hessian_index[i], hessian_index[j], hessian_size[i2], hessian_size[j2]) =
                        Eigen::Map<flower_focus::MatrixXd>(hessian55[i2 * O_Size + j2], hessian_size[i2], hessian_size[j2]);
                } else {
                    hessian_dense.block(hessian_index[i], hessian_index[j], hessian_size[i2], hessian_size[j2]) =
                        (Eigen::Map<flower_focus::MatrixXd>(hessian55[j2 * O_Size + i2], hessian_size[j2], hessian_size[i2])).transpose();
                }
            }
        }
    }



    if (lhs) {
        hessian_dense = hessian_dense.selfadjointView<Eigen::Upper>();
        Eigen::Map<flower_focus::MatrixXd>lhs_outside(lhs, order2p_local[O_FULL], order2p_local[O_FULL]);
        lhs_outside.setZero();
        lhs_outside.block(order2p_local[O_POSi], order2p_local[O_POSi], hessian_dense.cols(), hessian_dense.cols()) = hessian_dense;
        flower_focus::MatrixXd A = lhs_outside;
        flower_focus::MatrixXd D = lhs_outside;

        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            matrix_update(A, info.m, D, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
        }
        A = D;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            matrix_update(D, info.m, A, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }
        lhs_outside = A;
        lhs_outside = lhs_outside.selfadjointView<Eigen::Upper>();
        // LOG_OUT<<lhs_outside<<std::endl;
    }

    if (rhs) {
        // LOG_OUT<<rhs_dense<<std::endl;

        Eigen::Map<flower_focus::VectorXd>rhs_outside(rhs, order2p_local[O_FULL]);
        rhs_outside.setZero();
        rhs_outside.segment(order2p_local[O_POSi], rhs_dense.size()) = rhs_dense;

        Eigen::VectorXd b = rhs_outside;
        Eigen::VectorXd b2 = rhs_outside;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            vector_update(b, info.m, b2, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }

        rhs_outside = b2;
        // LOG_OUT<<rhs_outside.transpose()<<std::endl;
    }

    if (gradient) {
        Eigen::Map<flower_focus::VectorXd>gradient_outside(gradient, order2p_local[O_FULL]);
        gradient_outside.setZero();
        gradient_outside.segment(order2p_local[O_POSi], gradient_dense.size()) = gradient_dense;

        Eigen::VectorXd b = gradient_outside;
        Eigen::VectorXd b2 = gradient_outside;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            vector_update(b, info.m, b2, info.idx1, info.idx2,  RIGHT, info.is_identity, info.b_row, info.b_col);
        }

        gradient_outside = b2;
        // LOG_OUT<<gradient_outside.transpose()<<std::endl;
    }



}

void IMUGNSSBase::UpdateRhsN() {
    //+
    Eigen::VectorXd tmp(hessian_size[O_N]);
    MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 0>(
        phase_biases_hessians.data(), hessian_size[O_N], hessian_size[O_N],
        phase_biases_values.data(),
        tmp.data());

    Eigen::Map<Eigen::VectorXd>(rhs5[O_N], hessian_size[O_N]) += tmp;
    Eigen::Map<Eigen::VectorXd>(gradient5[O_N], hessian_size[O_N]) += tmp;

}

void IMUGNSSBase::UpdateRhsPose(int i) {

    flower_focus::Vector15d dx = GetInc(i);

    flower_focus::Vector15d tmp1;
    Eigen::VectorXd tmp2(hessian_size[O_N]);

    MatrixVectorMultiply<15, 15, 0>(
        pose_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_Pose1],
        dx.data(),
        tmp1.data());
    Eigen::Map<flower_focus::Vector15d>(rhs5[O_Pose1]) += tmp1;
    Eigen::Map<flower_focus::Vector15d>(gradient5[O_Pose1]) += tmp1;

    MatrixVectorMultiply<15, Eigen::Dynamic, 0>(
        pose_phase_biases_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_N],
        phase_biases_values.data(),
        tmp1.data());
    Eigen::Map<flower_focus::Vector15d>(rhs5[O_Pose1]) += tmp1;
    Eigen::Map<flower_focus::Vector15d>(gradient5[O_Pose1]) += tmp1;

    MatrixTransposeVectorMultiply<15, Eigen::Dynamic, 0>(
        pose_phase_biases_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_N],
        dx.data(),
        tmp2.data());
    Eigen::Map<Eigen::VectorXd>(rhs5[O_N], hessian_size[O_N]) += tmp2;
    Eigen::Map<Eigen::VectorXd>(gradient5[O_N], hessian_size[O_N]) += tmp2;


}

void IMUGNSSBase::EvaluateAlpha(double* alpha1) {
    Eigen::VectorXd model_residuals;
    if (gnss_speed_bias.size()) {
        model_residuals = imu_factors[0]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(gradient5[O_Pose0]) +
                          imu_factors[0]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(gradient_save[0]) +
                          imu_factors[0]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(gradient5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
        alpha1[0] += model_residuals.squaredNorm();
    } else {
        model_residuals = imu_factors[0]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(gradient5[O_Pose0]) +
                          imu_factors[0]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(gradient5[O_Pose2]) +
                          imu_factors[0]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(gradient5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
        alpha1[0] += model_residuals.squaredNorm();
    }

    Eigen::VectorXd gradientN;
    if (gnss_speed_bias.size()) {
        gradientN = Eigen::Map<flower_focus::VectorXd>(gradient5[O_N], hessian_size[O_N]);
        alpha1[0] += gradientN.transpose() * phase_biases_hessians * gradientN;
    }
    for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {
        if (i != (int)gnss_speed_bias.size() - 1) {
            model_residuals = imu_factors[i + 1]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(gradient_save[i]) +
                              imu_factors[i + 1]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(gradient_save[i + 1]) +
                              imu_factors[i + 1]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(gradient5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
            alpha1[0] += model_residuals.squaredNorm();
        } else {
            model_residuals = imu_factors[i + 1]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(gradient_save[i]) +
                              imu_factors[i + 1]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(gradient5[O_Pose2]) +
                              imu_factors[i + 1]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(gradient5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
            alpha1[0] += model_residuals.squaredNorm();
        }
        alpha1[0] += 2 * Eigen::Map<flower_focus::Vector15d>(gradient_save[i]).transpose() * pose_phase_biases_hessians[i] * gradientN;
        alpha1[0] += Eigen::Map<flower_focus::Vector15d>(gradient_save[i]).transpose() * pose_hessians[i] * Eigen::Map<flower_focus::Vector15d>(gradient_save[i]);

    }

}


void IMUGNSSBase::UpdateHiddenGaussStep(double* gauss_squared_norm, double* gradient_squared_norm, double* gradient_dot_gauss) {
    flower_focus::Vector15d tmp;
    memcpy(tmp.data(), gauss5[O_Pose2], 15 * sizeof(double));

    for (int i = (int)hmn_save[0].size() - 1; i >= 0; i--) {

        MatrixVectorMultiply<15, 15, 0>(
            hmn_save[O_Pose1][i], hessian_size[O_Pose1], hessian_size[O_Pose1],
            rhsmn_save[i],
            gauss5[O_Pose1]);

        for (int j = O_Pose2; j < O_Size; j++) {
            MatrixVectorMultiply < 15, Eigen::Dynamic, -1 > (
                hmn_save[j][i], hessian_size[O_Pose1], hessian_size[j],
                gauss5[j],
                gauss5[O_Pose1]);
        }



        memcpy(gauss_save[i], gauss5[O_Pose1], sizeof(double) * 15);
        memcpy(gauss5[O_Pose2], gauss5[O_Pose1], sizeof(double) * 15);

        gauss_squared_norm[0] += Eigen::Map<flower_focus::Vector15d>(gauss_save[i]).squaredNorm();
        gradient_squared_norm[0] += Eigen::Map<flower_focus::Vector15d>(gradient_save[i]).squaredNorm();
        gradient_dot_gauss[0] += Eigen::Map<flower_focus::Vector15d>(gradient_save[i]).dot(Eigen::Map<flower_focus::Vector15d>(gauss_save[i]));


    }
    ASSERT(!isnan(gauss_squared_norm[0]));

    memcpy(gauss5[O_Pose2], tmp.data(), 15 * sizeof(double));


}


void IMUGNSSBase::ResetMem() {
    for (int i = 0; i < O_Size; i++) {
        memset(rhs5[i], 0, sizeof(double)*hessian_size[i]);
        memset(gradient5[i], 0, sizeof(double)*hessian_size[i]);
        for (int j = i; j < O_Size; j++)
            memset(hessian55[i * O_Size + j], 0, sizeof(double) * (hessian_size[i]*hessian_size[j]));
    }
}

void IMUGNSSBase::UpdateTrustRegion(double mu) {
    Eigen::Map<flower_focus::Matrix15d>(hessian55[O_Pose1 * O_Size + O_Pose1]).diagonal().array() += mu;//diag.array().matrix();
}



bool IMUGNSSBase::EvaluateRhs(double* rhs, double* gradient) {
    ASSERT(history_flag);
#if BTOC_UPDATE
    p_jacobians.clear();
#endif

    for (int i = 0; i < O_Size; i++) {
        memset(rhs5[i], 0, sizeof(double)*hessian_size[i]);
        memset(gradient5[i], 0, sizeof(double)*hessian_size[i]);
    }
    if (gnss_speed_bias.size()) {
        CopyRhs2Rhs(rhs5, phase_biases_rhs.data(), O_N, hessian_size[O_N]);
        CopyRhs2Rhs(gradient5, phase_biases_rhs.data(), O_N, hessian_size[O_N]);
        UpdateRhsN();
    }


    JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose0, O_Pose1, O_EXTRINSIC}), std::vector<int>({15, 15, hessian_size[O_EXTRINSIC]}), imu_factors[0]->imu_jacobians_raw_in, 0, imu_factors[0]->imu_residual_in.data(), rhs5, gradient5);

    for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {

        JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose1, O_Pose2, O_EXTRINSIC}), std::vector<int>({15, 15, hessian_size[O_EXTRINSIC]}), imu_factors[i + 1]->imu_jacobians_raw_in, 0, imu_factors[i + 1]->imu_residual_in.data(), rhs5, gradient5);

        //update information relate to gnss factor
        UpdateRhsPose(i);
        RhsUpdateRhs(rhs5, pose_rhses[i].data(), O_Pose1, hessian_size[O_Pose1]);
        RhsUpdateRhs(gradient5, pose_rhses[i].data(), O_Pose1, hessian_size[O_Pose1]);
#if BTOC_UPDATE
        p_jacobians.push_back(PJacobian(r_wi_wc, p_wi_wc, scale_factor, Eigen::Quaterniond(gnss_poses[i] + 3).toRotationMatrix(), flower_focus::Vector3d(gnss_poses[i]), ric, tic, pbg));
#endif
        //marginalization
        MargPose1rhsOnly(hessian55, hessian_size, rhs5, i);
        MoveRhsData(i);

    }
    //update schur jacobian of selected states.
    UpdateSchurComponent(rhs, 0, gradient);
    last_gnss_cost = gnss_cost;



    return true;
}



flower_focus::Vector15d IMUGNSSBase::GetInc(int i) {
    flower_focus::Vector15d dx;

    Eigen::Map<flower_focus::Vector7d>x(gnss_poses[i], 7);
    Eigen::Map<flower_focus::Vector7d>x0(gnss_poses_lin[i], 7);
    dx.segment<3>(0) = x.head<3>() - x0.head<3>();
    dx.segment<3>(3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
    if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0))
        dx.segment<3>(3) = 2.0 * -(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
    Eigen::Map<flower_focus::Vector9d >x_(gnss_speed_bias[i], 9);
    Eigen::Map<flower_focus::Vector9d >x0_(gnss_speed_bias_lin[i], 9);
    dx.segment<9>(6) = x_ - x0_;

    return dx;
}




void IMUGNSSBase::Reparameterization(double const* const* parameters) {

    int p_idx = 0;
    int p_idx2 = 0;

    if (ESTIMATE_SCALE) {
        scale_factor = parameters[p_idx][0];
        p_idx++;
    } else scale_factor = 1;
    if (ESTIMATE_P_WI_WC || ESTIMATE_R_WI_WC) {
        p_wi_wc = Eigen::Vector3d(parameters[p_idx][0], parameters[p_idx][1], parameters[p_idx][2]);
        r_wi_wc = Eigen::Quaterniond(parameters[p_idx][6], parameters[p_idx][3], parameters[p_idx][4], parameters[p_idx][5]).toRotationMatrix();
        p_idx++;
    } else {
        r_wi_wc = R_WI_WC;
        p_wi_wc = P_WI_WC;
    }
    if (ESTIMATE_EXTRINSIC) {
        tic = Eigen::Vector3d(parameters[p_idx][0], parameters[p_idx][1], parameters[p_idx][2]);
        Eigen::Quaterniond qic = Eigen::Quaterniond(parameters[p_idx][6], parameters[p_idx][3], parameters[p_idx][4], parameters[p_idx][5]);
        ric = qic.toRotationMatrix();
        p_idx++;
    } else {
        tic = TIC[0];
        ric = RIC[0];
    }

    Pci = Eigen::Vector3d(parameters[p_idx][0], parameters[p_idx][1], parameters[p_idx][2]);
    Rci = Eigen::Quaterniond(parameters[p_idx][6], parameters[p_idx][3], parameters[p_idx][4], parameters[p_idx][5]).toRotationMatrix();
    p_idx++;
    for (int i = 0; i < 9; i++)Bi(i) = parameters[p_idx][i];
    p_idx++;


    Pcj = Eigen::Vector3d(parameters[p_idx][0], parameters[p_idx][1], parameters[p_idx][2]);
    Rcj = Eigen::Quaterniond(parameters[p_idx][6], parameters[p_idx][3], parameters[p_idx][4], parameters[p_idx][5]).toRotationMatrix();
    p_idx++;
    for (int i = 0; i < 9; i++)Bj(i) = parameters[p_idx][i];
    p_idx++;




    if (ESTIMATE_TD) {
        extrinsic(p_idx2) = parameters[p_idx][0]; p_idx2 += 1;
        p_idx++;
    }
    if (ESTIMATE_ACC_SCALE) {
        Eigen::Vector3d acc_scale;
        acc_scale = flower_focus::Vector3d(parameters[p_idx]); extrinsic.segment(p_idx2, 3) = acc_scale; p_idx2 += 3;
        p_idx++;
    }
    pbg = Pbg;
    if (ESTIMATE_PBG) {
        pbg = flower_focus::Vector3d(parameters[p_idx]); extrinsic.segment(p_idx2, 3) = pbg; p_idx2 += 3;
        p_idx++;
    }

    ASSERT(p_idx2 == extrinsic.size());

    for (int i = 0; i < phase_biases_values.size(); i++)phase_biases_values(i) = parameters[p_idx++][0];



    Pgi = r_wi_wc * Pci * scale_factor + p_wi_wc - r_wi_wc * Rci * ric.transpose() * (tic - pbg);
    Rbi = r_wi_wc * Rci * ric.transpose();
    Qbi = Eigen::Quaterniond(Rbi);



    Pgj = r_wi_wc * Pcj * scale_factor + p_wi_wc - r_wi_wc * Rcj * ric.transpose() * (tic - pbg);
    Rbj = r_wi_wc * Rcj * ric.transpose();
    Qbj = Eigen::Quaterniond(Rbj);

    PgQbi(0) = Pgi.x(); PgQbi(1) = Pgi.y(); PgQbi(2) = Pgi.z();
    PgQbi(3) = Qbi.x(); PgQbi(4) = Qbi.y(); PgQbi(5) = Qbi.z(); PgQbi(6) = Qbi.w();

    PgQbj(0) = Pgj.x(); PgQbj(1) = Pgj.y(); PgQbj(2) = Pgj.z();
    PgQbj(3) = Qbj.x(); PgQbj(4) = Qbj.y(); PgQbj(5) = Qbj.z(); PgQbj(6) = Qbj.w();







}

#if BTOC_UPDATE
void IMUGNSSBase::ReparameterizationNew(double const* const* parameters) {


    int p_idx = 0;

    if (ESTIMATE_SCALE) {
        new_scale = parameters[p_idx][0];
        p_idx++;
    } else
        new_scale = 1;
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC) {
        new_p_wi_wc = Eigen::Vector3d(parameters[p_idx][0], parameters[p_idx][1], parameters[p_idx][2]);
        new_r_wi_wc = Eigen::Quaterniond(parameters[p_idx][6], parameters[p_idx][3], parameters[p_idx][4], parameters[p_idx][5]).toRotationMatrix();
        p_idx++;
    } else {
        new_p_wi_wc = P_WI_WC;
        new_r_wi_wc = R_WI_WC;
    }
    if (ESTIMATE_EXTRINSIC) {
        new_tic = Eigen::Vector3d(parameters[p_idx][0], parameters[p_idx][1], parameters[p_idx][2]);
        Eigen::Quaterniond qic = Eigen::Quaterniond(parameters[p_idx][6], parameters[p_idx][3], parameters[p_idx][4], parameters[p_idx][5]);
        new_ric = qic.toRotationMatrix();
        p_idx++;
    } else {
        new_tic = TIC[0];
        new_ric = RIC[0];
    }
    p_idx += 4;

    if (ESTIMATE_TD)
        p_idx++;
    if (ESTIMATE_ACC_SCALE)
        p_idx++;
    new_pbg = Pbg;
    if (ESTIMATE_PBG) {
        new_pbg = flower_focus::Vector3d(parameters[p_idx]);
        p_idx++;
    }





}
#endif


void IMUGNSSBase::UpdateInsideStateUseGaussAndGradient(double const* const* parameters, double gradient_scale, double gauss_scale) {

    for (int order = O_Pose1; order < O_Size; order++) {
        Eigen::Map<flower_focus::VectorXd>dx(inc5[order], hessian_size[order]);
        dx = Eigen::Map<flower_focus::VectorXd>(gauss5[order], hessian_size[order]) * gauss_scale - Eigen::Map< flower_focus::VectorXd>(gradient5[order], hessian_size[order]) * gradient_scale;
    }
#if BTOC_UPDATE
    ASSERT(p_jacobians.size() == hmn_save[0].size());
    if (ESTIMATE_PBG)
        pbg_inc = pbg_gauss * gauss_scale - pbg_gradient * gradient_scale;
    if (ESTIMATE_SCALE)
        scale_inc = scale_gauss * gauss_scale - scale_gradient * gradient_scale;
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
        posgloabl_inc = posgloabl_gauss * gauss_scale - posgloabl_gradient * gradient_scale;
    if (ESTIMATE_EXTRINSIC)
        extrinsic_inc = extrinsic_gauss * gauss_scale - extrinsic_gradient * gradient_scale;

#endif
    for (int i = (int)hmn_save[0].size() - 1; i >= 0; i--) {

        Eigen::Map<flower_focus::Vector15d>dx(inc_save[i]);
        dx = flower_focus::Vector15d(gauss_save[i]) * gauss_scale - flower_focus::Vector15d(gradient_save[i]) * gradient_scale;

#if BTOC_UPDATE
        flower_focus::Vector6d new_dx = p_jacobians[i].dP0_dP0 * dx.segment(0, 6);
        if (ESTIMATE_EXTRINSIC)new_dx += p_jacobians[i].dP0_dextrinsic * extrinsic_inc;
        if (ESTIMATE_PBG)new_dx += p_jacobians[i].dP0_dpbg * pbg_inc;
        if (ESTIMATE_SCALE)new_dx += p_jacobians[i].dP0_dscale * scale_inc;
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)new_dx += p_jacobians[i].dP0_dR_WI_WC * posgloabl_inc;


        Eigen::Map<flower_focus::Vector3d> Pgi(gnss_poses[i]);
        Eigen::Map<Eigen::Quaterniond> Qbi(gnss_poses[i] + 3);


        Eigen::Quaterniond Qci = Eigen::Quaterniond(p_jacobians[i].r_wi_wc.transpose() * Qbi.toRotationMatrix() * p_jacobians[i].ric);
        Eigen::Vector3d Pci = (p_jacobians[i].r_wi_wc.transpose() * (Pgi - p_jacobians[i].p_wi_wc) + Qci.toRotationMatrix() * p_jacobians[i].ric.transpose() * (p_jacobians[i].tic - p_jacobians[i].pbg)) / p_jacobians[i].scale_factor;


        Pci -= Eigen::Map<flower_focus::Vector3d>(new_dx.data());
        Qci = (  Qci *  Utility::deltaQ(-Eigen::Map<flower_focus::Vector3d>(new_dx.data() + 3))  ).normalized();


        Pgi = new_r_wi_wc * Pci * new_scale + new_p_wi_wc - new_r_wi_wc * Qci.toRotationMatrix() * new_ric.transpose() * (new_tic - new_pbg);
        Rbi = new_r_wi_wc * Qci.toRotationMatrix() * new_ric.transpose();
        Qbi = Eigen::Quaterniond(Rbi);

#else
        Eigen::Map<flower_focus::Vector3d> P(gnss_poses[i]);
        Eigen::Map<Eigen::Quaterniond> Q(gnss_poses[i] + 3);

        P -= Eigen::Map<flower_focus::Vector3d>(dx.data());
        Q = (  Q *  Utility::deltaQ(-Eigen::Map<flower_focus::Vector3d>(dx.data() + 3))  ).normalized();
#endif



        Eigen::Map<flower_focus::Vector9d> B(gnss_speed_bias[i]);
        B -=  Eigen::Map<flower_focus::Vector9d>(dx.data() + 6);
    }



}

bool IMUGNSSBase::EvaluateCost(double* cost, double* model_cost_change) {

    if (cost) {
        if (gnss_speed_bias.size()) {
            gnss_cost = gnss_init_cost;
            gnss_cost += 0.5 * phase_biases_values.transpose() * phase_biases_hessians * phase_biases_values;

            gnss_cost += phase_biases_values.transpose() * phase_biases_rhs;
            // ASSERT(gnss_cost == 0);
            for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {
                flower_focus::Vector15d dx = GetInc(i);
                gnss_cost += 0.5 * dx.transpose() * pose_hessians[i] * dx;
                gnss_cost += dx.transpose() * pose_phase_biases_hessians[i] * phase_biases_values;

                gnss_cost += dx.transpose() * pose_rhses[i];
            }
            ASSERT(gnss_cost >= 0);
            cost[0] += gnss_cost;
        }


        for (int i = 0; i < (int)imu_factors.size(); i++) {
            imu_factors[i]->factor->Evaluate(imu_factors[i]->parameters.data(), imu_factors[i]->imu_residual_in.data(), 0);
            cost[0] += 0.5 * imu_factors[i]->imu_residual_in.squaredNorm();
            // if (imu_factors.size())LOG_OUT << "3:" << 0.5 * imu_factors[i]->imu_residual_in.squaredNorm() << std::endl;
        }

    }

    if (model_cost_change) {
        if (gnss_speed_bias.size()) {
            for (int i = 0; i < (int)gnss_phase_biases.size(); i++)
                phase_biases_values(i) = gnss_phase_biases[i][0];
            double gnss_cost = gnss_init_cost;
            gnss_cost += 0.5 * phase_biases_values.transpose() * phase_biases_hessians * phase_biases_values;

            gnss_cost += phase_biases_values.transpose() * phase_biases_rhs;
            // ASSERT(gnss_cost == 0);
            for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {
                flower_focus::Vector15d dx = GetInc(i);
                gnss_cost += 0.5 * dx.transpose() * pose_hessians[i] * dx;
                gnss_cost += dx.transpose() * pose_phase_biases_hessians[i] * phase_biases_values;

                gnss_cost += dx.transpose() * pose_rhses[i];
            }
            ASSERT(gnss_cost >= 0);
            if(!USE_GNSS)ASSERT(last_gnss_cost - gnss_cost==0);
            model_cost_change[0] += last_gnss_cost - gnss_cost;
            // LOG_OUT << "a:" << last_gnss_cost - gnss_cost << std::endl;
        }

        Eigen::VectorXd model_residuals;
        if (gnss_speed_bias.size()) {
            model_residuals = imu_factors[0]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(inc5[O_Pose0]) +
                              imu_factors[0]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(inc_save[0]) +
                              imu_factors[0]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(inc5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
            imu_factors[0]->model_residual_accum += model_residuals;
            model_cost_change[0] += model_residuals.transpose() * (imu_factors[0]->imu_residual_in_save - model_residuals / 2);
            // LOG_OUT << "b:" << model_residuals.transpose() * (imu_factors[0]->imu_residual_in_save - model_residuals / 2) << std::endl;
        } else {
            model_residuals = imu_factors[0]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(inc5[O_Pose0]) +
                              imu_factors[0]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(inc5[O_Pose2]) +
                              imu_factors[0]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(inc5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
            imu_factors[0]->model_residual_accum += model_residuals;
            model_cost_change[0] += model_residuals.transpose() * (imu_factors[0]->imu_residual_in_save - model_residuals / 2);
            ASSERT(!isnan(model_residuals.norm()));
            ASSERT(!isnan(imu_factors[0]->model_residual_accum.norm()));
        }




        for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {
            if (i != (int)gnss_speed_bias.size() - 1) {
                model_residuals = imu_factors[i + 1]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(inc_save[i]) +
                                  imu_factors[i + 1]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(inc_save[i + 1]) +
                                  imu_factors[i + 1]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(inc5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
                imu_factors[i + 1]->model_residual_accum += model_residuals;
                model_cost_change[0] += model_residuals.transpose() * (imu_factors[i + 1]->imu_residual_in_save - model_residuals / 2);
                // LOG_OUT << "d:" << model_residuals.transpose() * (imu_factors[i + 1]->imu_residual_in_save - model_residuals / 2) << std::endl;
            } else {
                model_residuals = imu_factors[i + 1]->imu_jacobians_in[0] * Eigen::Map<flower_focus::Vector15d>(inc_save[i]) +
                                  imu_factors[i + 1]->imu_jacobians_in[1] * Eigen::Map<flower_focus::Vector15d>(inc5[O_Pose2]) +
                                  imu_factors[i + 1]->imu_jacobians_in[2] * Eigen::Map<flower_focus::VectorXd>(inc5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
                imu_factors[i + 1]->model_residual_accum += model_residuals;
                model_cost_change[0] += model_residuals.transpose() * (imu_factors[i + 1]->imu_residual_in_save - model_residuals / 2);
                // LOG_OUT << "e:" << model_residuals.transpose() * (imu_factors[i + 1]->imu_residual_in_save - model_residuals / 2) << std::endl;
            }
            ASSERT(!isnan(model_residuals.norm()));
            ASSERT(!isnan(imu_factors[i + 1]->model_residual_accum.norm()));


        }
    }




    return true;


}

void IMUGNSSBase::SaveLoadVariable(bool is_save) {
    if (is_save) {
        for (int i = 0; i < (int)gnss_poses_old.size(); i++)memcpy(gnss_poses_old[i].data(), gnss_poses[i], sizeof(double) * 7);
        for (int i = 0; i < (int)gnss_speed_bias.size(); i++)memcpy(gnss_speed_bias_old[i].data(), gnss_speed_bias[i], sizeof(double) * 9);
        for (int i = 0; i < (int)imu_factors.size(); i++) {
            imu_factors[i]->model_residual_accum_save = imu_factors[i]->model_residual_accum;
            imu_factors[i]->imu_residual_in_save = imu_factors[i]->imu_residual_in;
        }
    } else {
        for (int i = 0; i < (int)gnss_poses_old.size(); i++)memcpy(gnss_poses[i], gnss_poses_old[i].data(), sizeof(double) * 7);
        for (int i = 0; i < (int)gnss_speed_bias.size(); i++)memcpy(gnss_speed_bias[i], gnss_speed_bias_old[i].data(), sizeof(double) * 9);

        for (int i = 0; i < (int)imu_factors.size(); i++) {
            imu_factors[i]->model_residual_accum = imu_factors[i]->model_residual_accum_save;
            imu_factors[i]->imu_residual_in = imu_factors[i]->imu_residual_in_save;
        }
    }
}


bool IMUGNSSBase::EvaluateLhsRhs( double* rhs, double* lhs, double* gradient, double mu) {

    BuildMatrixInfo();
    ResetMem();
    history_flag = true;
#if BTOC_UPDATE
    p_jacobians.clear();
#endif
    if (gnss_speed_bias.size()) {
        CopyHessian2Hessian(hessian55, phase_biases_hessians.data(), O_N, O_N, hessian_size[O_N], hessian_size[O_N]);
        CopyRhs2Rhs(rhs5, phase_biases_rhs.data(), O_N, hessian_size[O_N]);
        CopyRhs2Rhs(gradient5, phase_biases_rhs.data(), O_N, hessian_size[O_N]);
        UpdateRhsN();
    }


    //update  information about first imu factor

    imu_factors[0]->factor->Evaluate(imu_factors[0]->parameters.data(), 0, imu_factors[0]->imu_jacobians_raw_in);
    JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose0, O_Pose1, O_EXTRINSIC}),
                                     std::vector<int>({15, 15, hessian_size[O_EXTRINSIC]}),
                                     imu_factors[0]->imu_jacobians_raw_in, hessian55, imu_factors[0]->imu_residual_in.data(), rhs5, gradient5);

    for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {

        imu_factors[i + 1]->factor->Evaluate(imu_factors[i + 1]->parameters.data(), 0, imu_factors[i + 1]->imu_jacobians_raw_in);
        JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose1, O_Pose2, O_EXTRINSIC}),
                                         std::vector<int>({15, 15, hessian_size[O_EXTRINSIC]}),
                                         imu_factors[i + 1]->imu_jacobians_raw_in, hessian55, imu_factors[i + 1]->imu_residual_in.data(), rhs5, gradient5);

        //update information relate to gnss factor
        UpdateRhsPose(i);
        HessianUpdateHessian(hessian55, pose_hessians[i].data(), O_Pose1, O_Pose1, hessian_size[O_Pose1], hessian_size[O_Pose1]);
        HessianUpdateHessian(hessian55, pose_phase_biases_hessians[i].data(), O_Pose1, O_N, hessian_size[O_Pose1], hessian_size[O_N]);
        RhsUpdateRhs(rhs5, pose_rhses[i].data(), O_Pose1, hessian_size[O_Pose1]);
        RhsUpdateRhs(gradient5, pose_rhses[i].data(), O_Pose1, hessian_size[O_Pose1]);
        UpdateTrustRegion(mu);

#if BTOC_UPDATE
        p_jacobians.push_back(PJacobian(r_wi_wc, p_wi_wc, scale_factor, Eigen::Quaterniond(gnss_poses[i] + 3).toRotationMatrix(), flower_focus::Vector3d(gnss_poses[i]), ric, tic, pbg));
#endif
        //marginalization
        MargPose1(hessian55, hessian_size, rhs5);
        MoveRhsData(i);
        MoveHessianData(i);

    }
    //update schur jacobian of selected states.
    UpdateSchurComponent(rhs, lhs, gradient);
    // ASSERT(Eigen::Map<flower_focus::VectorXd>(rhs5[O_N], hessian_size[O_N]).norm() == 0);
    for (int i = 0; i < (int)imu_factors.size(); i++) {
        imu_factors[i]->imu_residual_in_linerized = imu_factors[i]->imu_residual_in;
        imu_factors[i]->model_residual_accum.setZero();
    }
    last_gnss_cost = gnss_cost;




    return true;
}


void IMUGNSSBase::EvaluateGaussAndGradientStep( const double* gradient, const double* gauss, double* gauss_squared_norm, double* gradient_squared_norm, double* gradient_dot_gauss, double* alpha) {


    Eigen::VectorXd gauss_inside;
    Eigen::VectorXd gradient_inside;

    {
        Eigen::VectorXd gauss_outside = Eigen::Map<const Eigen::VectorXd>(gauss, order2p_local[O_FULL]);
        // LOG_OUT<<gauss_outside.transpose()<<std::endl;
#if BTOC_UPDATE
        if (ESTIMATE_PBG)
            pbg_gauss = gauss_outside.segment(order2p_local[O_PBG], 3);
        if (ESTIMATE_SCALE)
            scale_gauss = gauss_outside.segment(order2p_local[O_SCALE], 1);
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
            posgloabl_gauss = gauss_outside.segment(order2p_local[O_POSGLOBAL], 6);
        if (ESTIMATE_EXTRINSIC)
            extrinsic_gauss = gauss_outside.segment(order2p_local[O_EXTRINSIC_IC], 6);
#endif

        ASSERT(!isnan(gauss_outside.norm()));
        Eigen::VectorXd b2 = gauss_outside;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            vector_update(gauss_outside, info.m, b2, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
        }
        gauss_inside = b2;
        ASSERT(!isnan(gauss_inside.norm()));

    }
    {
        Eigen::VectorXd gradient_outside = Eigen::Map<const Eigen::VectorXd>(gradient, order2p_local[O_FULL]);

#if BTOC_UPDATE
        if (ESTIMATE_PBG)
            pbg_gradient = gradient_outside.segment(order2p_local[O_PBG], 3);
        if (ESTIMATE_SCALE)
            scale_gradient = gradient_outside.segment(order2p_local[O_SCALE], 1);
        if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
            posgloabl_gradient = gradient_outside.segment(order2p_local[O_POSGLOBAL], 6);
        if (ESTIMATE_EXTRINSIC)
            extrinsic_gradient = gradient_outside.segment(order2p_local[O_EXTRINSIC_IC], 6);

#endif

        Eigen::VectorXd b2 = gradient_outside;
        for (int i = 0; i < (int)matrix_info.size(); ++i) {
            auto& info = matrix_info[i];
            vector_update(gradient_outside, info.m, b2, info.idx1, info.idx2,  LEFT, info.is_identity, info.b_row, info.b_col);
        }
        gradient_inside = b2;
    }

    int delta_count = 0;
    if (ESTIMATE_SCALE)
        delta_count += 1;
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
        delta_count += 6;
    if (ESTIMATE_EXTRINSIC)
        delta_count += 6;
    {
        Eigen::Map<flower_focus::Vector15d>dx(gauss5[O_Pose0]);
        dx = gauss_inside.segment(delta_count, hessian_size[O_Pose0]);
        delta_count += hessian_size[O_Pose0];
    }
    {
        Eigen::Map<flower_focus::Vector15d>dx(gauss5[O_Pose2]);
        dx = gauss_inside.segment(delta_count, hessian_size[O_Pose2]);
        delta_count += hessian_size[O_Pose2];
    }
    {
        Eigen::Map<flower_focus::VectorXd>dx(gauss5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
        dx = gauss_inside.segment(delta_count, hessian_size[O_EXTRINSIC]);
        delta_count += hessian_size[O_EXTRINSIC];
    }
    {
        Eigen::Map<flower_focus::VectorXd>dx(gauss5[O_N], hessian_size[O_N]);
        dx = gauss_inside.segment(delta_count, hessian_size[O_N]);
        delta_count += hessian_size[O_N];
    }

    delta_count = 0;
    if (ESTIMATE_SCALE)
        delta_count += 1;
    if (ESTIMATE_R_WI_WC || ESTIMATE_P_WI_WC)
        delta_count += 6;
    if (ESTIMATE_EXTRINSIC)
        delta_count += 6;
    {
        Eigen::Map<flower_focus::Vector15d>dx(gradient5[O_Pose0]);
        dx = gradient_inside.segment(delta_count, hessian_size[O_Pose0]);
        delta_count += hessian_size[O_Pose0];
    }
    {
        Eigen::Map<flower_focus::Vector15d>dx(gradient5[O_Pose2]);
        dx = gradient_inside.segment(delta_count, hessian_size[O_Pose2]);
        delta_count += hessian_size[O_Pose2];
    }
    {
        Eigen::Map<flower_focus::VectorXd>dx(gradient5[O_EXTRINSIC], hessian_size[O_EXTRINSIC]);
        dx = gradient_inside.segment(delta_count, hessian_size[O_EXTRINSIC]);
        delta_count += hessian_size[O_EXTRINSIC];
    }
    {
        Eigen::Map<flower_focus::VectorXd>dx(gradient5[O_N], hessian_size[O_N]);
        dx = gradient_inside.segment(delta_count, hessian_size[O_N]);
        delta_count += hessian_size[O_N];
    }


    UpdateHiddenGaussStep(gauss_squared_norm, gradient_squared_norm, gradient_dot_gauss);
    EvaluateAlpha(alpha);


}


IMUGNSSFactor::IMUGNSSFactor(IMUGNSSBase* IMUGNSS_info_): IMUGNSS_info(IMUGNSS_info_) {


}



bool IMUGNSSFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians)  const {



    return 1;

}















