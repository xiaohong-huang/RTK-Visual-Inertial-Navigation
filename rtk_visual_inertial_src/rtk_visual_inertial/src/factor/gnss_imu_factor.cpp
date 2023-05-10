#include "gnss_imu_factor.h"
#include "ceres/invert_psd_matrix.h"
#include "ceres/small_blas.h"
#include "../utility/utility.h"
#include "../parameter/parameters.h"

const double min_diagonal = 1e-6;
const double max_diagonal = 1e32;
const double eps = 1e-8;
const int mapindex[O_Size - 1] = {O_Pose0, O_Pose1, O_N};

#define Matrix15d_ Eigen::Matrix<double,15,15,Eigen::RowMajor>
#define Vector15d_ Eigen::Matrix<double,15,1,Eigen::ColMajor>
#define MatrixXd_ Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
#define VectorXd_ Eigen::Matrix<double,Eigen::Dynamic,1,Eigen::ColMajor>
#define MatrixX15d_ Eigen::Matrix<double,Eigen::Dynamic,15,Eigen::RowMajor>
#define Matrix15Xd_ Eigen::Matrix<double,15,Eigen::Dynamic,Eigen::RowMajor>
#define Vector7d_ Eigen::Matrix<double,7,1,Eigen::ColMajor>
#define Vector9d_ Eigen::Matrix<double,9,1,Eigen::ColMajor>
#define Vector3d_ Eigen::Matrix<double,3,1,Eigen::ColMajor>



IMUGNSSBase::~IMUGNSSBase() {
    Clear();
    for (int i = 0; i < (int)imu_factors.size(); i++) {
        delete imu_factors[i];
    }
    delete last_imu_factor;
    if (gnss_middle_marginfo)delete gnss_middle_marginfo;
}

void IMUGNSSBase::Init() {
    InitHessianRhs(hessian55, hessian_size, rhs5);
    phase_biases_values.resize(hessian_size[O_N]);
    history_flag = false;
    INC = VectorXd_(hessian_size[O_N] + 30);

    hessian_index[0] = 0;
    global_size = 0;
    for (int i = 0; i < O_Size - 1 - 1; i++) {
        hessian_index[i + 1] = hessian_index[i] + hessian_size[mapindex[i]];
    }
    for (int i = 0; i < O_Size - 1; i++) {
        global_size += hessian_size[mapindex[i]];
    }
}

void IMUGNSSBase::Clear() {
    for (int i = 0; i < O_Size; i++) {
        if (rhs5[i])delete rhs5[i];
        if (delta5[i])delete delta5[i];
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
    }



}

void IMUGNSSBase::InitHessianRhs(double** hessian55, int* hessian_size, double** rhs5) {
    hessian_size[0] = 15;
    hessian_size[1] = 15;
    hessian_size[2] = (int)gnss_phase_biases.size();
    hessian_size[3] = 15;


    for (int i = 0; i < O_Size; i++) {
        rhs5[i] = new double[hessian_size[i]];
        delta5[i] = new double[hessian_size[i]];
        for (int j = i; j < O_Size; j++) {
            int size = hessian_size[i] * hessian_size[j];
            hessian55[i * O_Size + j] = new double[size];
        }
    }
    rhsmn_save.resize(gnss_speed_bias.size());
    for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {
        rhsmn_save[i] = new double [hessian_size[O_Pose1]];
    }

    for (int i = O_Pose1; i < O_Size; i++) {
        hmn_save[i].resize(gnss_speed_bias.size());
        for (int j = 0; j < (int)gnss_speed_bias.size(); j++) {
            hmn_save[i][j] = new double[ hessian_size[O_Pose1] * hessian_size[i] ];
        }
    }

}

void IMUGNSSBase::SetLastImuFactor(IntegrationBase* pre_integration, double* para_Pose2, double* para_SpeedBias2) {
    if (last_imu_factor)delete last_imu_factor;
    last_imu_factor = new IMUFactor(pre_integration);

    param.clear();
    param = std::vector<double*>({para_pose0, para_speed_bias0, para_Pose2, para_SpeedBias2});
    for (int i = 0; i < (int)gnss_phase_biases.size(); i++) {
        param.push_back(gnss_phase_biases[i]);
    }

    if (USE_GLOBAL_OPTIMIZATION) {
        for (int i = 0; i < (int)gnss_poses.size(); i++) {
            if (my_problem->HasParameterBlock(gnss_poses[i]))my_problem->RemoveParameterBlock(gnss_poses[i]);
            if (my_problem->HasParameterBlock(gnss_speed_bias[i]))my_problem->RemoveParameterBlock(gnss_speed_bias[i]);
        }
        IMUGNSSFactor* factor = new IMUGNSSFactor(this);
        residualBlockId = my_problem->AddResidualBlock(factor, 0, param);
        // factor->Evaluate(param.data(),0,0);
    }

}

void IMUGNSSBase::AddMidMargInfo(MarginalizationInfo* gnss_middle_marginfo_) {

    assert(gnss_middle_marginfo == 0);
    gnss_middle_marginfo = gnss_middle_marginfo_;


    const MatrixXd_ &GNSS_IMU_Hessian = gnss_middle_marginfo->A;
    const VectorXd_ &GNSS_IMU_rhs = gnss_middle_marginfo->b;

    int set_num = 0;

    for (int i = 0; i < (int)imu_factors.size(); i++) {
        IntegrationBase* pre_integration = imu_factors[i]->pre_integration;
        if (pre_integration->sum_dt <= IMUINTEGRATIONTIMEMAX + 100000)continue;
        set_num++;
        gnss_Index = i;
    }
    assert(set_num == 1);
    int N_size_external = 0;
    int oldNsize = (int)gnss_phase_biases.size();
    std::map<double*, int>parameter_index;
    for (int j = 0; j < (int)gnss_middle_marginfo->keep_block_addr.size(); j++) {
        double* pointer = gnss_middle_marginfo->keep_block_addr[j];
        int size = gnss_middle_marginfo->keep_block_size[j];
        if (size == 7) {
            assert(pointer == gnss_poses[gnss_Index] || pointer == gnss_poses[gnss_Index - 1]);
        } else if (size == 9) {
            assert(pointer == gnss_speed_bias[gnss_Index] || pointer == gnss_speed_bias[gnss_Index - 1]);
        } else if (size == 1) {
            if (gnss_phase_biases_set.find(pointer) == gnss_phase_biases_set.end()) {
                gnss_phase_biases_set.insert(pointer);
                gnss_phase_biases.push_back(pointer);
                N_size_external++;
            }
        } else {
            assert(0);
        }
        parameter_index[pointer] = gnss_middle_marginfo->keep_block_idx[j] - gnss_middle_marginfo->m;
    }
    int newNsize = (int)gnss_phase_biases.size();
    assert(oldNsize + N_size_external == newNsize);
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


    std::vector<double*>hessian55_m_pointer({gnss_poses[gnss_Index - 1], gnss_speed_bias[gnss_Index - 1], gnss_poses[gnss_Index], gnss_speed_bias[gnss_Index]});
    std::vector<int>hessian55_m_index({O_Pose1, O_Pose1, O_Pose2, O_Pose2});
    std::vector<int>hessian55_m_size({6, 9, 6, 9});
    std::vector<int>hessian55_m_shift({0, 6, 0, 6});

    for (int i = 0; i < (int)gnss_phase_biases.size(); i++) {
        hessian55_m_pointer.push_back(gnss_phase_biases[i]);
        hessian55_m_index.push_back(O_N);
        hessian55_m_size.push_back(1);
        hessian55_m_shift.push_back(i);
    }

    for (int i = 0; i < (int)hessian55_m_pointer.size(); i++) {
        double* pointer1 = hessian55_m_pointer[i];
        if (parameter_index.find(pointer1) == parameter_index.end())continue;
        int index1 = parameter_index[pointer1];
        int shift1 = hessian55_m_shift[i];
        int size1 = hessian55_m_size[i];
        int indexm1 = hessian55_m_index[i];
        if (indexm1 == O_N) {
            phase_biases_rhs(shift1) += GNSS_IMU_rhs(index1);
        } else if (indexm1 == O_Pose1) {
            pose_rhses[gnss_Index - 1].segment(shift1, size1) += GNSS_IMU_rhs.segment(index1, size1);
        } else if (indexm1 == O_Pose2) {
            pose_rhses[gnss_Index].segment(shift1, size1) += GNSS_IMU_rhs.segment(index1, size1);
        } else {
            assert(0);
        }
        for (int j = i; j < (int)hessian55_m_pointer.size(); j++) {
            double* pointer2 = hessian55_m_pointer[j];
            if (parameter_index.find(pointer2) == parameter_index.end())continue;
            int index2 = parameter_index[pointer2];
            int shift2 = hessian55_m_shift[j];
            int size2 = hessian55_m_size[j];
            int indexm2 = hessian55_m_index[j];
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix = GNSS_IMU_Hessian.block(index1, index2, size1, size2);
            if (indexm1 == O_N && indexm2 == O_N) {
                phase_biases_hessians.block(shift1, shift2, size1, size2) += matrix;
            } else if (indexm1 == O_Pose1 && indexm2 == O_N) {
                pose_phase_biases_hessians[gnss_Index - 1].block(shift1, shift2, size1, size2) += matrix;
            } else if (indexm1 == O_Pose2 && indexm2 == O_N) {
                pose_phase_biases_hessians[gnss_Index].block(shift1, shift2, size1, size2) += matrix;
            } else if (indexm1 == O_Pose1 && indexm2 == O_Pose1) {
                pose_hessians[gnss_Index - 1].block(shift1, shift2, size1, size2) += matrix;
            } else if (indexm1 == O_Pose2 && indexm2 == O_Pose2) {
                pose_hessians[gnss_Index].block(shift1, shift2, size1, size2) += matrix;
            } else if (indexm1 == O_Pose1 && indexm2 == O_Pose2) {
                pose1_pose2_hessians.block(shift1, shift2, size1, size2) = matrix;
            } else {
                assert(0);
            }

        }
    }

    phase_biases_hessians = phase_biases_hessians.selfadjointView<Eigen::Upper>();
    // pose1_pose2_hessians=pose1_pose2_hessians.selfadjointView<Eigen::Upper>();
    // std::cout<<pose1_pose2_hessians<<std::endl;
    pose_hessians[gnss_Index - 1] = pose_hessians[gnss_Index - 1].selfadjointView<Eigen::Upper>();
    pose_hessians[gnss_Index] = pose_hessians[gnss_Index].selfadjointView<Eigen::Upper>();
    Init();




}


void IMUGNSSBase::AddMargInfo(MarginalizationInfo* marg_info_gnss, IntegrationBase* pre_integration,
                              double* para_Poses, double* para_speed_bias) {
    Clear();
    residualBlockId = 0;
    imu_factors.push_back(new IMUFactor(pre_integration));
    int N_size_external = 0;
    std::map<double*, int>parameter_index;
    gnss_poses.push_back(para_Poses);
    gnss_speed_bias.push_back(para_speed_bias);
    gnss_speed_bias_lin.push_back(para_speed_bias);
    gnss_poses_lin.push_back(para_Poses);

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
            assert(fabs(marg_info_gnss->keep_block_data[i][0]) == 0);
            if (gnss_phase_biases_set.find(pointer) == gnss_phase_biases_set.end()) {
                gnss_phase_biases_set.insert(pointer);
                gnss_phase_biases.push_back(pointer);
                N_size_external++;
            }
        }
    }
    int newNsize = (int)gnss_phase_biases.size();
    assert(oldNsize + N_size_external == newNsize);
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

    Matrix15d_ Pose_Hessian;
    MatrixXd_ Pose_N_Hessian(15, newNsize);
    Vector15d_ Pose_rhs;

    Pose_Hessian.setZero();
    Pose_N_Hessian.setZero();
    Pose_rhs.setZero();

    const MatrixXd_ &GNSS_Hessian = marg_info_gnss->A;
    const VectorXd_ &GNSS_rhs = marg_info_gnss->b;

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

void IMUGNSSBase::JacobianResidualUpdateHessianRhs(std::vector<int>index, std::vector<int>mapsize, double** jacobian, double** hessian55, double* residual, double** rhs5) {
    for (int i = 0; i < (int)index.size(); i++) {
        int sizei = mapsize[i];
        int indexi = index[i];
        ceres::internal::MatrixTransposeVectorMultiply<15, Eigen::Dynamic, 1>
        (jacobian[i], 15, sizei, residual, rhs5[indexi]);
        for (int j = 0; j < (int)index.size(); j++) {
            int sizej = mapsize[j];
            int indexj = index[j];
            double* hessian_local = hessian55[indexi * O_Size + indexj];
            if (indexj < indexi)continue;
            if (sizei == 15 && sizej == 15) {
                ceres::internal::MatrixTransposeMatrixMultiply<15, 15, 15, 15, 1>
                (jacobian[i], 15, 15, jacobian[j], 15, 15, hessian_local, 0, 0, 15, 15);
            } else {
                assert(0);
            }
        }
    }
}

void IMUGNSSBase::HessianUpdateHessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2) {
    assert(block2 >= block1);
    if (size1 == 15 && size2 == 15) {
        Eigen::Map<Matrix15d_>(hessian_desc[block1 * O_Size + block2]) +=
            Eigen::Map<Matrix15d_>(hessian_src);
    } else {
        Eigen::Map<MatrixXd_>(hessian_desc[block1 * O_Size + block2], size1, size2) +=
            Eigen::Map<MatrixXd_>(hessian_src, size1, size2);
    }
}

void IMUGNSSBase::RhsUpdateRhs(double** rhs_desc, double* rhs_src, int block1, int size1) {
    if (size1 == 15) {
        Eigen::Map<Vector15d_>(rhs_desc[block1]) += Eigen::Map<Vector15d_>(rhs_src);
    } else {
        Eigen::Map<VectorXd_ >(rhs_desc[block1], size1) += Eigen::Map<VectorXd_ > (rhs_src, size1);
    }
}

void IMUGNSSBase::CopyHessian2Hessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2) {
    assert(block2 >= block1);
    memcpy(hessian_desc[block1 * O_Size + block2], hessian_src, sizeof(double)*size1 * size2);
}

void IMUGNSSBase::CopyRhs2Rhs(double** rhs_desc, double* rhs_src, int block1, int size1) {
    memcpy(rhs_desc[block1], rhs_src, sizeof(double)*size1);
}

void IMUGNSSBase::MargPose1(double** hessian55, int* hessian_size, double** rhs5) {

    Matrix15d_ Amm_inverse = ceres::internal::InvertPSDMatrix<15>(1,
                             Eigen::Map<Matrix15d_>(hessian55[O_Pose1 * O_Size + O_Pose1]) );

    memcpy(hessian55[O_Pose1 * O_Size + O_Pose1], Amm_inverse.data(), sizeof(double) * 15 * 15);

    for (int i = O_Pose1 + 1; i < O_Size; i++) {
        int sizen = hessian_size[i];
        MatrixX15d_ Anm_Amm_inverse(sizen, 15);
        ceres::internal::MatrixTransposeMatrixMultiply<15, Eigen::Dynamic, 15, 15, 0>(
            hessian55[O_Pose1 * O_Size + i], 15, sizen,
            Amm_inverse.data(), 15, 15,
            Anm_Amm_inverse.data(), 0, 0, sizen, 15);

        ceres::internal::MatrixVectorMultiply < Eigen::Dynamic, 15, -1 > (
            Anm_Amm_inverse.data(), sizen, 15,
            rhs5[O_Pose1],
            rhs5[i]);

        for (int j = i; j < O_Size; j++) {
            int sizev = hessian_size[j];
            ceres::internal::MatrixMatrixMultiply
            < Eigen::Dynamic, 15, 15, Eigen::Dynamic, -1 > (
                Anm_Amm_inverse.data(), sizen, 15,
                hessian55[O_Pose1 * O_Size + j], 15, sizev,
                hessian55[i * O_Size + j], 0, 0, sizen, sizev);
        }
    }

}

void IMUGNSSBase::MoveHessianData(int index) {

    memcpy(rhsmn_save[index], rhs5[O_Pose1], sizeof(double)*hessian_size[O_Pose1]);
    for (int i = O_Pose1; i < O_Size; i++) {
        memcpy(hmn_save[i][index], hessian55[O_Pose1 * O_Size + i], sizeof(double)*hessian_size[O_Pose1]*hessian_size[i]);
    }
    memcpy(hessian55[O_Pose1 * O_Size + O_Pose1], hessian55[O_Pose2 * O_Size + O_Pose2], sizeof(double)*hessian_size[O_Pose2]*hessian_size[O_Pose2]);
    memset(hessian55[O_Pose2 * O_Size + O_Pose2], 0, sizeof(double)*hessian_size[O_Pose2]*hessian_size[O_Pose2]);
    memcpy(rhs5[O_Pose1], rhs5[O_Pose2], sizeof(double)*hessian_size[O_Pose2]);
    memset(rhs5[O_Pose2], 0, sizeof(double)*hessian_size[O_Pose2]);

    for (int i = O_Pose2 + 1; i < O_Size; i++) {
        memcpy(hessian55[O_Pose1 * O_Size + i], hessian55[O_Pose2 * O_Size + i], sizeof(double)*hessian_size[O_Pose2]*hessian_size[i]);
        memset(hessian55[O_Pose2 * O_Size + i], 0, sizeof(double)*hessian_size[O_Pose2]*hessian_size[i]);
    }
    memset(hessian55[O_Pose1 * O_Size + O_Pose2], 0, sizeof(double)*hessian_size[O_Pose1]*hessian_size[O_Pose2]);

}

void IMUGNSSBase::UpdateSchurComponent() {



    hessian_dense = MatrixXd_ (global_size, global_size);
    rhs_dense = VectorXd_(global_size);

    for (int i = 0; i < O_Size - 1; i++) {
        int i2 = mapindex[i];
        rhs_dense.segment(hessian_index[i], hessian_size[i2]) = Eigen::Map<VectorXd_ >(rhs5[i2], hessian_size[i2]);
        for (int j = i; j < O_Size - 1; j++) {
            int j2 = mapindex[j];
            if (j2 >= i2) {
                hessian_dense.block(hessian_index[i], hessian_index[j], hessian_size[i2], hessian_size[j2]) =
                    Eigen::Map<MatrixXd_>(hessian55[i2 * O_Size + j2], hessian_size[i2], hessian_size[j2]);
            } else {
                hessian_dense.block(hessian_index[i], hessian_index[j], hessian_size[i2], hessian_size[j2]) =
                    (Eigen::Map<MatrixXd_>(hessian55[j2 * O_Size + i2], hessian_size[j2], hessian_size[i2])).transpose();
            }
        }
    }

    hessian_dense = hessian_dense.selfadjointView<Eigen::Upper>();


    Eigen::SelfAdjointEigenSolver<MatrixXd_> saes2(hessian_dense);
    VectorXd_ S = VectorXd_((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
    VectorXd_ S_sqrt = S.cwiseSqrt();
    schur_jacobian = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();


    VectorXd_ S_inv = VectorXd_((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));
    VectorXd_ S_inv_sqrt = S_inv.cwiseSqrt();
    schur_residual = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * rhs_dense;

}
void IMUGNSSBase::UpdateJacobResidual(double** jacobians, double* residuals) {
    int size = schur_residual.rows();
    if (residuals) {
        if (!history_flag || update_flag)
            Eigen::Map<VectorXd_>(residuals, size) = schur_residual;
        else
            Eigen::Map<VectorXd_>(residuals, size) = schur_residual - schur_jacobian * INC;
    }


    // std::cout<<"schur_residual2:"<<Eigen::Map<VectorXd_>(residuals,size).norm()<<std::endl;
    if (jacobians) {
        if (jacobians[0]) {
            Eigen::Map<MatrixXd_>tmp(jacobians[0], size, 7);
            tmp.setZero();
            tmp.block(0, 0, size, 6) = schur_jacobian.block(0, hessian_index[0], size, 6);
        }
        if (jacobians[1]) {
            Eigen::Map<MatrixXd_>(jacobians[1], size, 9) = schur_jacobian.block(0, hessian_index[0] + 6, size, 9);
        }
        if (jacobians[2]) {
            Eigen::Map<MatrixXd_>tmp(jacobians[2], size, 7);
            tmp.setZero();
            tmp.block(0, 0, size, 6) = schur_jacobian.block(0, hessian_index[1], size, 6);
        }
        if (jacobians[3]) {
            Eigen::Map<MatrixXd_>(jacobians[3], size, 9) = schur_jacobian.block(0, hessian_index[1] + 6, size, 9);
        }
        for (int i = 0; i < hessian_size[O_N]; i++) {
            if (jacobians[4 + i]) {
                Eigen::Map<MatrixXd_>(jacobians[4 + i], size, 1) = schur_jacobian.block(0, hessian_index[2] + i, size, 1);
            }
        }
    }

}

void IMUGNSSBase::UpdateRhsN() {
    //+
    ceres::internal::MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 1>(
        phase_biases_hessians.data(), hessian_size[O_N], hessian_size[O_N],
        phase_biases_values.data(),
        rhs5[O_N]);
}

void IMUGNSSBase::UpdateRhsPose(int i) {

    Vector15d_ dx = GetInc(i);


    ceres::internal::MatrixVectorMultiply<15, 15, 1>(
        pose_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_Pose1],
        dx.data(),
        rhs5[O_Pose1]);

    ceres::internal::MatrixVectorMultiply<15, Eigen::Dynamic, 1>(
        pose_phase_biases_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_N],
        phase_biases_values.data(),
        rhs5[O_Pose1]);

    ceres::internal::MatrixTransposeVectorMultiply<15, Eigen::Dynamic, 1>(
        pose_phase_biases_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_N],
        dx.data(),
        rhs5[O_N]);

    // VectorXd_ rt(phase_biases_values.rows());
    // ceres::internal::MatrixTransposeVectorMultiply<15,Eigen::Dynamic, 0>(
    //     pose_phase_biases_hessians[i].data(), hessian_size[O_Pose1], hessian_size[O_N],
    //     dx.data(),
    //     rt.data());



}

void IMUGNSSBase::UpdateDeltaValues(Vector7d_ &Pi, Vector9d_ &Bi, Vector7d_ &Pj, Vector9d_ &Bj, VectorXd_ &phase_biases_values) {
    double delta_sum = 0;
    int delta_count = 0;
    {
        Eigen::Map<Vector15d_>dx(delta5[O_Pose2]);
        Vector7d_ x = Pj;
        Vector7d_ x0 = Pj_old;
        dx.segment<3>(0) =  x0.head<3>() - x.head<3>();
        dx.segment<3>(3) = -2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
        if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0))
            dx.segment<3>(3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
        dx.segment<9>(6) = Bj_old - Bj;
        INC.segment(15, 15) = dx;
        delta_sum += dx.norm();
        delta_count += 15;
    }

    {
        Eigen::Map<VectorXd_ >dx(delta5[O_N], hessian_size[O_N]);
        dx = phase_biases_values_old - phase_biases_values;
        INC.segment(30, hessian_size[O_N]) = dx;
        delta_sum += dx.norm();
        delta_count += hessian_size[O_N];
    }

    {
        Eigen::Map<Vector15d_>dx(delta5[O_Pose0]);
        Vector7d_ x = Pi;
        Vector7d_ x0 = Pi_old;
        dx.segment<3>(0) = x0.head<3>() - x.head<3>() ;
        dx.segment<3>(3) = -2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
        if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0))
            dx.segment<3>(3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
        dx.segment<9>(6) = Bi_old - Bi;
        INC.segment(0, 15) = dx;
        delta_sum += dx.norm();
        delta_count += 15;
    }


}

void IMUGNSSBase::UpdateHiddenState() {

    for (int i = (int)hmn_save[0].size() - 1; i >= 0; i--) {

        for (int j = O_Pose2; j < O_Size; j++) {
            ceres::internal::MatrixVectorMultiply < 15, Eigen::Dynamic, -1 > (
                hmn_save[j][i], hessian_size[O_Pose1], hessian_size[j],
                delta5[j],
                rhsmn_save[i]);
        }

        ceres::internal::MatrixVectorMultiply<15, 15, 0>(
            hmn_save[O_Pose1][i], hessian_size[O_Pose1], hessian_size[O_Pose1],
            rhsmn_save[i],
            delta5[O_Pose2]);


        Eigen::Map<Vector3d_> P(gnss_poses[i]);
        Eigen::Map<Eigen::Quaterniond> Q(gnss_poses[i] + 3);

        P -= Eigen::Map<Vector3d_>(delta5[O_Pose2]);
        Q = (  Q *  Utility::deltaQ(-Eigen::Map<Vector3d_>(delta5[O_Pose2] + 3))  ).normalized();


        Eigen::Map<Vector9d_> B(gnss_speed_bias[i]);
        B -=  Eigen::Map<Vector9d_>(delta5[O_Pose2] + 6);
    }




}

void IMUGNSSBase::ResetMem() {
    for (int i = 0; i < O_Size; i++) {
        memset(rhs5[i], 0, sizeof(double)*hessian_size[i]);
        for (int j = i; j < O_Size; j++) {
            memset(hessian55[i * O_Size + j], 0, sizeof(double) * (hessian_size[i]*hessian_size[j]));
        }
    }
}

void IMUGNSSBase::UpdateTrustRegion() {
    Vector15d_ diag = Eigen::Map<Matrix15d_>(hessian55[O_Pose1 * O_Size + O_Pose1]).diagonal().array();
    for (int i = 0; i < 15; ++i) {
        diag(i) = std::min(std::max(diag(i), min_diagonal), max_diagonal);
    }
    diag = (diag / MAX_TRUST_REGION_RADIUS).array();

    Eigen::Map<Matrix15d_>(hessian55[O_Pose1 * O_Size + O_Pose1]).diagonal() += diag.array().matrix();
}


void IMUGNSSBase::SaveLastStates(Vector7d_ &Pi, Vector9d_ &Bi, Vector7d_ &Pj, Vector9d_ &Bj) {
    Pi_old = Pi;
    Bi_old = Bi;
    Pj_old = Pj;
    Bj_old = Bj;
    phase_biases_values_old = phase_biases_values;
}

Vector15d_ IMUGNSSBase::GetInc(int i) {
    Vector15d_ dx;

    Eigen::Map<Vector7d_>x(gnss_poses[i], 7);
    Eigen::Map<Vector7d_>x0(gnss_poses_lin[i], 7);
    dx.segment<3>(0) = x.head<3>() - x0.head<3>();
    dx.segment<3>(3) = 2.0 * (Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
    if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0)) {
        dx.segment<3>(3) = 2.0 * -(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
    }
    Eigen::Map<Vector9d_ >x_(gnss_speed_bias[i], 9);
    Eigen::Map<Vector9d_ >x0_(gnss_speed_bias_lin[i], 9);
    dx.segment<9>(6) = x_ - x0_;

    return dx;
}
bool IMUGNSSBase::Evaluate(double const* const* parameters, double* residuals, double** jacobians) {

    Matrix15d_ J1;
    Matrix15d_ J2;
    Vector7d_ Pi;
    Vector9d_ Bi;
    Vector7d_ Pj;
    Vector9d_ Bj;

    double* jacobian[3] = {J1.data(), J2.data(), 0};
    double residual[15];
    evaluate_round++;
    for (int i = 0; i < hessian_size[O_N]; i++)phase_biases_values(i) = parameters[4 + i][0];
    for (int i = 0; i < 7; i++) {
        Pi(i) = parameters[0][i];
        Pj(i) = parameters[2][i];
    }
    for (int i = 0; i < 9; i++) {
        Bi(i) = parameters[1][i];
        Bj(i) = parameters[3][i];
    }
    if (!history_flag)SaveLastStates(Pi, Bi, Pj, Bj);

    //update hedden states
    UpdateDeltaValues(Pi, Bi, Pj, Bj, phase_biases_values);

    if (jacobians) update_flag = true;
    else update_flag = false;

    if (history_flag && update_flag)UpdateHiddenState();

    if (!history_flag || update_flag) {

        update_round++;
        // std::cout<<"update jacobian\r\n";
        history_flag = true;
        SaveLastStates(Pi, Bi, Pj, Bj);
        ResetMem();
        //update information about phase biases
        CopyHessian2Hessian(hessian55, phase_biases_hessians.data(), O_N, O_N, hessian_size[O_N], hessian_size[O_N]);
        CopyRhs2Rhs(rhs5, phase_biases_rhs.data(), O_N, hessian_size[O_N]);
        UpdateRhsN();

        //update  information about first imu factor
        const double* parameter[4] = {parameters[0], parameters[1], gnss_poses[0], gnss_speed_bias[0]};
        imu_factors[0]->Evaluate2(parameter, residual, jacobian);

        JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose0, O_Pose1}),
                                         std::vector<int>({15, 15}), jacobian, hessian55, residual, rhs5);

        for (int i = 0; i < (int)gnss_speed_bias.size(); i++) {

            //update information relate to imu factor
            const double* parameter[4];
            if (i != (int)gnss_speed_bias.size() - 1) {
                parameter[0] = gnss_poses[i];
                parameter[1] = gnss_speed_bias[i];
                parameter[2] = gnss_poses[i + 1];
                parameter[3] = gnss_speed_bias[i + 1];


                if (imu_factors[i + 1]->pre_integration->sum_dt <= IMUINTEGRATIONTIMEMAX + 100000) {
                    imu_factors[i + 1]->Evaluate2(parameter, residual, jacobian);
                    JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose1, O_Pose2}),
                                                     std::vector<int>({15, 15}), jacobian, hessian55, residual, rhs5);
                } else {
                    assert(i == gnss_Index - 1);
                    Vector15d_ inc1 = GetInc(i);
                    Vector15d_ inc2 = GetInc(i + 1);

                    ceres::internal::MatrixVectorMultiply<15, 15, 1>(
                        pose1_pose2_hessians.data(), hessian_size[O_Pose1], hessian_size[O_Pose2],
                        inc2.data(),
                        rhs5[O_Pose1]);


                    ceres::internal::MatrixTransposeVectorMultiply<15, 15, 1>(
                        pose1_pose2_hessians.data(), hessian_size[O_Pose1], hessian_size[O_Pose2],
                        inc1.data(),
                        rhs5[O_Pose2]);

                    HessianUpdateHessian(hessian55, pose1_pose2_hessians.data(), O_Pose1, O_Pose2, hessian_size[O_Pose1], hessian_size[O_Pose2]);
                }

            } else {
                parameter[0] = gnss_poses[i];
                parameter[1] = gnss_speed_bias[i];
                parameter[2] = parameters[2];
                parameter[3] = parameters[3];
                last_imu_factor->Evaluate2(parameter, residual, jacobian);
                JacobianResidualUpdateHessianRhs(std::vector<int>({O_Pose1, O_Pose2}),
                                                 std::vector<int>({15, 15}), jacobian, hessian55, residual, rhs5);

            }



            //update information relate to gnss factor
            UpdateRhsPose(i);
            HessianUpdateHessian(hessian55, pose_hessians[i].data(), O_Pose1, O_Pose1, hessian_size[O_Pose1], hessian_size[O_Pose1]);
            HessianUpdateHessian(hessian55, pose_phase_biases_hessians[i].data(), O_Pose1, O_N, hessian_size[O_Pose1], hessian_size[O_N]);
            RhsUpdateRhs(rhs5, pose_rhses[i].data(), O_Pose1, hessian_size[O_Pose1]);
            // UpdateTrustRegion();

            //marginalization
            MargPose1(hessian55, hessian_size, rhs5);
            MoveHessianData(i);

        }
        //update schur jacobian of selected states.
        UpdateSchurComponent();

    }

    //update jacobian of selected states.
    UpdateJacobResidual(jacobians, residuals);
    // std::cout<<"schur_residual"<<schur_residual.norm()<<std::endl;

    // std::cout<<"t.toc:"<<t.toc()<<std::endl;


    return true;
}


IMUGNSSFactor::IMUGNSSFactor(IMUGNSSBase* IMUGNSS_info_): IMUGNSS_info(IMUGNSS_info_) {

    int cnt = 0;
    mutable_parameter_block_sizes()->push_back(7);
    mutable_parameter_block_sizes()->push_back(9);
    mutable_parameter_block_sizes()->push_back(7);
    mutable_parameter_block_sizes()->push_back(9);
    cnt += 30;
    for (int i = 0; i < IMUGNSS_info->hessian_size[O_N]; i++) {
        mutable_parameter_block_sizes()->push_back(1);
        cnt++;
    }
    set_num_residuals(cnt);

}
bool IMUGNSSFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {

    return IMUGNSS_info->Evaluate(parameters, residuals, jacobians);
}















