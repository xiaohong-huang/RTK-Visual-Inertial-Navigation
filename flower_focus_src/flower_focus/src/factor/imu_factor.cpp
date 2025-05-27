#include"imu_factor.h"
#include "../utility/utility.h"
#include "../parameter/parameters.h"
#include "integration_base.h"



// ----------------------------------------------------------------------------------------------------------------------Pb,Vg
bool IMUFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {

#define LEN 6
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);

    Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
    Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
    Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);


    int p_idx = 0;

    double td = 0;
    Eigen::Vector3d acc_scale = Eigen::Vector3d({1, 1, 1});
    Eigen::Vector3d pbg = Pbg;
    Eigen::Matrix3d Rwgw2 = Rwgw;

    if (ESTIMATE_TD) {
        td = parameters[4][p_idx];
        p_idx++;
    }
    if (ESTIMATE_ACC_SCALE) {
        acc_scale = Eigen::Vector3d(parameters[4][p_idx + 0], parameters[4][p_idx + 1], parameters[4][p_idx + 2]);
        p_idx += 3;
    }
    if (ESTIMATE_PBG) {
        pbg = Eigen::Vector3d(parameters[4][p_idx + 0], parameters[4][p_idx + 1], parameters[4][p_idx + 2]);
        p_idx += 3;
    }

    // LOG_OUT<<Rwgw2<<std::endl;


    ASSERT(pre_integration-> dt_i != 0 && pre_integration->dt_j != 0);


    Eigen::Matrix3d dp_dba = pre_integration->jacobian.block<3, 3>(O_P, O_BA);
    Eigen::Matrix3d dp_dbg = pre_integration->jacobian.block<3, 3>(O_P, O_BG);
    Eigen::Matrix3d dq_dbg = pre_integration->jacobian.block<3, 3>(O_R, O_BG);
    Eigen::Matrix3d dv_dba = pre_integration->jacobian.block<3, 3>(O_V, O_BA);
    Eigen::Matrix3d dv_dbg = pre_integration->jacobian.block<3, 3>(O_V, O_BG);
    Eigen::Matrix<double, 15, 15> sqrt_info = pre_integration->get_sqrtinfo();
    double sum_dt = pre_integration->sum_dt;


    Eigen::Matrix3d Ri = Qi.toRotationMatrix();
    Eigen::Matrix3d Rj = Qj.toRotationMatrix();
    Eigen::Matrix3d Ri_inv = Ri.transpose();
    Eigen::Quaterniond Qij = Qi.inverse() * Qj;

    Eigen::Vector3d newG =  Rwgw2 * G;
    Eigen::Vector3d newVi =  Vi;
    Eigen::Vector3d newVj =  Vj;

    Eigen::Vector3d un_gyri = pre_integration->gyri - Bgi;
    Eigen::Vector3d un_gyrj = pre_integration->gyrj - Bgj;


    double new_dt_i;
    double new_dt_j;

    if (ESTIMATE_TD) {
        new_dt_i = pre_integration->dt_i + td;
        new_dt_j = pre_integration->dt_j + td;
    } else {
        new_dt_i = pre_integration->dt_i;
        new_dt_j = pre_integration->dt_j;
    }


    Eigen::Quaterniond dQti(1, un_gyri(0) * new_dt_i / 2, un_gyri(1) *new_dt_i / 2, un_gyri(2) * new_dt_i / 2);
    Eigen::Quaterniond dQtj(1, un_gyrj(0) * new_dt_j / 2, un_gyrj(1) *new_dt_j / 2, un_gyrj(2) * new_dt_j / 2);
    Eigen::Matrix3d dRti = dQti.toRotationMatrix();


    Eigen::Vector3d dP = 0.5 * newG * sum_dt * sum_dt +  Pj - Pi - Rj * pbg + newVj * new_dt_j - newVi * (new_dt_i + sum_dt)  - Rj * (Utility::skewSymmetric(un_gyrj) * pbg) * new_dt_j;
    Eigen::Vector3d dV = newG * sum_dt + newVj - Rj * (Utility::skewSymmetric(un_gyrj) * pbg)  - newVi;

    Eigen::Vector3d dba = Bai - pre_integration->linearized_ba;
    Eigen::Vector3d dbg = Bgi - pre_integration->linearized_bg;

    Eigen::Vector3d dacc_scale = acc_scale - pre_integration->acc_scale;

    Eigen::Matrix3d dp_dacc_scale = pre_integration->jacobian_s.block<3, 3>(O_P, 0);
    Eigen::Matrix3d dv_dacc_scale = pre_integration->jacobian_s.block<3, 3>(O_V, 0);

    Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * dbg);
    Eigen::Vector3d corrected_delta_v = pre_integration->delta_v + dv_dba * dba + dv_dbg * dbg + dv_dacc_scale * dacc_scale;
    Eigen::Vector3d corrected_delta_p = pre_integration->delta_p + dp_dba * dba + dp_dbg * dbg + dp_dacc_scale * dacc_scale;


    if (residuals) {

        Eigen::Map<Eigen::Matrix<double, 15, 1>> residual(residuals);
        residual.block<3, 1>(O_P, 0) = dRti.transpose() * ( Ri_inv * dP + pbg + (Utility::skewSymmetric(un_gyri) * pbg) * (sum_dt + new_dt_i) ) - corrected_delta_p;
        residual.block<3, 1>(O_R, 0) = 2 * (corrected_delta_q.inverse() * dQti.inverse()  * Qij * dQtj).vec();
        residual.block<3, 1>(O_V, 0) = dRti.transpose() * (Ri_inv * dV + Utility::skewSymmetric(un_gyri) * pbg)  - corrected_delta_v;
        residual.block<3, 1>(O_BA, 0) = Baj - Bai;
        residual.block<3, 1>(O_BG, 0) = Bgj - Bgi;
        // LOG_OUT<<residual.transpose()<<std::endl;
        residual = sqrt_info * residual;
    }


    if (jacobians) {

        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
            jacobian_pose_i.setZero();
            jacobian_pose_i.block<3, 3>(O_P, O_P) = -(dRti.transpose() * Ri_inv);
            jacobian_pose_i.block<3, 3>(O_P, O_R) = dRti.transpose() * Utility::skewSymmetric(Ri_inv * dP);
            jacobian_pose_i.block<3, 3>(O_P, O_V) = -(dRti.transpose() * Ri_inv * (sum_dt + new_dt_i));
            jacobian_pose_i.block<3, 3>(O_P, O_BA) = -dp_dba;
            jacobian_pose_i.block<3, 3>(O_P, O_BG) = -dp_dbg
                                                     - Utility::skewSymmetric(dRti.transpose() * ( Ri_inv * dP + pbg + (Utility::skewSymmetric(un_gyri) * pbg) * (sum_dt + new_dt_i) ) * new_dt_i)
                                                     + dRti.transpose() * (Utility::skewSymmetric(pbg)) * (sum_dt + new_dt_i);

            jacobian_pose_i.block<3, 3>(O_R, O_R) =
                -(Utility::Qleft(dQtj.inverse() * Qij.inverse()) * Utility::Qright(dQti * corrected_delta_q)).bottomRightCorner<3, 3>();
            jacobian_pose_i.block<3, 3>(O_R, O_BG) =
                -(Utility::Qleft(dQtj.inverse() * Qij.inverse() * dQti * pre_integration->delta_q).bottomRightCorner<3, 3>() * dq_dbg)
                + (new_dt_i * Utility::Qleft(corrected_delta_q.inverse()) * Utility::Qright(dQti.inverse() * Qij * dQtj)).bottomRightCorner<3, 3>();

            jacobian_pose_i.block<3, 3>(O_V, O_R) = dRti.transpose() * Utility::skewSymmetric(Ri_inv * dV);
            jacobian_pose_i.block<3, 3>(O_V, O_V) = -(dRti.transpose() * Ri_inv);
            jacobian_pose_i.block<3, 3>(O_V, O_BA) = -dv_dba;
            jacobian_pose_i.block<3, 3>(O_V, O_BG) = -dv_dbg
                                                     - Utility::skewSymmetric(dRti.transpose() * (Ri_inv * dV + Utility::skewSymmetric(un_gyri) * pbg) * new_dt_i)
                                                     + dRti.transpose() * Utility::skewSymmetric(pbg) ;

            jacobian_pose_i.block<3, 3>(O_BA, O_BA) = -Eigen::Matrix3d::Identity();
            jacobian_pose_i.block<3, 3>(O_BG, O_BG) = -Eigen::Matrix3d::Identity();
            jacobian_pose_i = sqrt_info * jacobian_pose_i;
        }
        if (jacobians[1]) {
            Eigen::Map<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
            jacobian_pose_j.setZero();
            jacobian_pose_j.block<3, 3>(O_P, O_P) = dRti.transpose() * Ri_inv ;
            jacobian_pose_j.block<3, 3>(O_P, O_R) = dRti.transpose() * Ri_inv * Rj * Utility::skewSymmetric(pbg + (Utility::skewSymmetric(un_gyrj) * pbg) * new_dt_j) ;
            jacobian_pose_j.block<3, 3>(O_P, O_V) = dRti.transpose() * Ri_inv * new_dt_j;
            jacobian_pose_j.block<3, 3>(O_P, O_BG) = - dRti.transpose() * Ri_inv * Rj * (Utility::skewSymmetric(pbg)) * new_dt_j;

            jacobian_pose_j.block<3, 3>(O_R, O_R) = (Utility::Qright(dQtj) * Utility::Qleft(corrected_delta_q.inverse() * dQti.inverse() * Qij)).bottomRightCorner<3, 3>();
            jacobian_pose_j.block<3, 3>(O_R, O_BG) = -(new_dt_j * Utility::Qleft(corrected_delta_q.inverse()  * dQti.inverse() * Qij * dQtj)).bottomRightCorner<3, 3>(); //待定

            jacobian_pose_j.block<3, 3>(O_V, O_V) = dRti.transpose() * Ri_inv;
            jacobian_pose_j.block<3, 3>(O_V, O_R) = dRti.transpose() * Ri_inv * Rj * Utility::skewSymmetric((Utility::skewSymmetric(un_gyrj) * pbg));
            jacobian_pose_j.block<3, 3>(O_V, O_BG) = - dRti.transpose() * Ri_inv * Rj * Utility::skewSymmetric(pbg) ;

            jacobian_pose_j.block<3, 3>(O_BA, O_BA) = Eigen::Matrix3d::Identity();
            jacobian_pose_j.block<3, 3>(O_BG, O_BG) = Eigen::Matrix3d::Identity();
            jacobian_pose_j = sqrt_info * jacobian_pose_j;
        }


        if (jacobians[2]) {
            int extrinsic_localsize = 0;
            if (ESTIMATE_TD)extrinsic_localsize += 1;
            if (ESTIMATE_ACC_SCALE)extrinsic_localsize += 3;
            if (ESTIMATE_PBG)extrinsic_localsize += 3;
            if (extrinsic_localsize) {
                Eigen::Map<flower_focus::MatrixXd> jacobian_extrinsic(jacobians[2], 15, extrinsic_localsize);
                jacobian_extrinsic.setZero();
                int p_idx = 0;
                if (ESTIMATE_TD) {
                    jacobian_extrinsic.block<3, 1>(O_P, p_idx) = dRti.transpose() * (Ri_inv * (newVj - newVi - Rj * (Utility::skewSymmetric(un_gyrj) * pbg)) + Utility::skewSymmetric(un_gyri) * pbg)
                                                                 + Utility::skewSymmetric(dRti.transpose() * ( Ri_inv * dP + pbg + (Utility::skewSymmetric(un_gyri) * pbg) * (sum_dt + new_dt_i) )) * un_gyri;
                    Eigen::Quaterniond a(0, -un_gyri(0), -un_gyri(1), -un_gyri(2) );
                    Eigen::Quaterniond b(0, un_gyrj(0), un_gyrj(1), un_gyrj(2) );
                    Eigen::Quaterniond c = dQti.inverse() * Qij *  dQtj;
                    jacobian_extrinsic.block<3, 1>(O_R, p_idx) =  (corrected_delta_q.inverse() * a * c).vec() + (corrected_delta_q.inverse() * c * b).vec();
                    jacobian_extrinsic.block<3, 1>(O_V, p_idx) = Utility::skewSymmetric(dRti.transpose() * (Ri_inv * dV + Utility::skewSymmetric(un_gyri) * pbg)) * un_gyri;
                    jacobian_extrinsic.block<15, 1>(0, p_idx) = sqrt_info * jacobian_extrinsic.block<15, 1>(0, p_idx);
                    p_idx += 1;
                }
                if (ESTIMATE_ACC_SCALE) {
                    if ( has_excitation) {
                        ASSERT((pre_integration->jacobian_s.block<3, 3>(O_R, 0)).norm() == 0);
                        ASSERT((pre_integration->jacobian_s.block<3, 3>(O_BA, 0)).norm() == 0);
                        ASSERT((pre_integration->jacobian_s.block<3, 3>(O_BG, 0)).norm() == 0);
                        jacobian_extrinsic.block<3, 3>(O_P, p_idx) = -dp_dacc_scale;
                        jacobian_extrinsic.block<3, 3>(O_V, p_idx) = -dv_dacc_scale;
                        ASSERT((jacobian_extrinsic.block<15, 3>(0, p_idx) + (pre_integration->jacobian_s.block<15, 3>(0, 0))).norm() == 0);
                        jacobian_extrinsic.block<15, 3>(0, p_idx) = sqrt_info * jacobian_extrinsic.block<15, 3>(0, p_idx);
                    }
                    p_idx += 3;
                }
                if (ESTIMATE_PBG) {
                    jacobian_extrinsic.block<3, 3>(O_P, p_idx) = dRti.transpose() * (
                                                                     Ri_inv * -Rj * (flower_focus::Matrix3d::Identity() + Utility::skewSymmetric(un_gyrj * new_dt_j)) + flower_focus::Matrix3d::Identity() + Utility::skewSymmetric(un_gyri * (sum_dt + new_dt_i))
                                                                 );
                    jacobian_extrinsic.block<3, 3>(O_V, p_idx) = dRti.transpose() * (
                                                                     Ri_inv * -Rj * Utility::skewSymmetric(un_gyrj) + Utility::skewSymmetric(un_gyri)
                                                                 );
                    jacobian_extrinsic.block<15, 3>(0, p_idx) = sqrt_info * jacobian_extrinsic.block<15, 3>(0, p_idx);


                    p_idx += 3;
                }


            }

        }
    }
#undef LEN

    return true;


}




