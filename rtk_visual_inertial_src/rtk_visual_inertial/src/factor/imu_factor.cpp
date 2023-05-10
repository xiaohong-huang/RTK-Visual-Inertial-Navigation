#include"imu_factor.h"
#include "../utility/utility.h"
#include "../parameter/parameters.h"

bool IMUFactor::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {

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



    Eigen::Map<Eigen::Matrix<double, 15, 1>> residual(residuals);
    residual = pre_integration->evaluate(Pi, Qi, Vi, Bai, Bgi, Pj, Qj, Vj, Baj, Bgj, Pbg);

    Eigen::Matrix<double, 15, 15> sqrt_info = pre_integration->get_sqrtinfo();
    residual = sqrt_info * residual;
    Eigen::Vector3d newG = Rwgw * G;
    if (jacobians) {
        double sum_dt = pre_integration->sum_dt;
        Eigen::Matrix3d dp_dba = pre_integration->jacobian.template block<3, 3>(O_P, O_BA);
        Eigen::Matrix3d dp_dbg = pre_integration->jacobian.template block<3, 3>(O_P, O_BG);

        Eigen::Matrix3d dq_dbg = pre_integration->jacobian.template block<3, 3>(O_R, O_BG);

        Eigen::Matrix3d dv_dba = pre_integration->jacobian.template block<3, 3>(O_V, O_BA);
        Eigen::Matrix3d dv_dbg = pre_integration->jacobian.template block<3, 3>(O_V, O_BG);

        if (pre_integration->jacobian.maxCoeff() > 1e8 || pre_integration->jacobian.minCoeff() < -1e8) {
            printf("numerical unstable in preintegration");
        }
        Eigen::Matrix3d Ri_inv = Qi.inverse().toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();



        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
            jacobian_pose_i.setZero();
            jacobian_pose_i.block<3, 3>(O_P, O_P) = -Ri_inv;
            jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * newG * sum_dt * sum_dt + ((Pj - Pi) - Qj * Pbg) - Vi * sum_dt));
            Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
            jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();
            jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (newG * sum_dt + (Vj - Qj * (Utility::skewSymmetric(pre_integration->gyrj - Bgj) * Pbg)) - Vi));
            jacobian_pose_i = sqrt_info * jacobian_pose_i;
            if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8) {
                printf("numerical unstable in preintegration");

            }
        }
        if (jacobians[1]) {
            Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
            jacobian_speedbias_i.setZero();
            jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Ri_inv * sum_dt;
            jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
            jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg + Utility::skewSymmetric(Pbg) * sum_dt;
            jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration->delta_q).bottomRightCorner<3, 3>() * dq_dbg;
            jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Ri_inv;
            jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
            jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg + Utility::skewSymmetric(Pbg);
            jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();
            jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();
            jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

        }
        if (jacobians[2]) {
            Eigen::Map<Eigen::Matrix<double, 15, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
            jacobian_pose_j.setZero();
            jacobian_pose_j.block<3, 3>(O_P, O_P) = Ri_inv;
            jacobian_pose_j.block<3, 3>(O_P, O_R) = Ri_inv * Rj * Utility::skewSymmetric(Pbg);
            Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
            jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();
            jacobian_pose_j.block<3, 3>(O_V, O_R) = Ri_inv * Rj * Utility::skewSymmetric(Utility::skewSymmetric(pre_integration->gyrj - Bgj) * Pbg);
            jacobian_pose_j = sqrt_info * jacobian_pose_j;

        }
        if (jacobians[3]) {
            Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
            jacobian_speedbias_j.setZero();
            jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Ri_inv;
            jacobian_speedbias_j.block<3, 3>(O_V, O_BG - O_V) = Ri_inv * Rj * (-Utility::skewSymmetric(Pbg));
            jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();
            jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();
            jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;

        }

    }

    return true;
}

bool IMUFactor::Evaluate2(double const* const* parameters, double* residuals, double** jacobians) {

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



    Eigen::Map<Eigen::Matrix<double, 15, 1>> residual(residuals);
    residual = pre_integration->evaluate(Pi, Qi, Vi, Bai, Bgi, Pj, Qj, Vj, Baj, Bgj, Pbg);

    Eigen::Matrix<double, 15, 15> sqrt_info = pre_integration->get_sqrtinfo();
    residual = sqrt_info * residual;
    Eigen::Vector3d newG = Rwgw * G;
    if (jacobians) {
        double sum_dt = pre_integration->sum_dt;
        Eigen::Matrix3d dp_dba = pre_integration->jacobian.template block<3, 3>(O_P, O_BA);
        Eigen::Matrix3d dp_dbg = pre_integration->jacobian.template block<3, 3>(O_P, O_BG);

        Eigen::Matrix3d dq_dbg = pre_integration->jacobian.template block<3, 3>(O_R, O_BG);

        Eigen::Matrix3d dv_dba = pre_integration->jacobian.template block<3, 3>(O_V, O_BA);
        Eigen::Matrix3d dv_dbg = pre_integration->jacobian.template block<3, 3>(O_V, O_BG);

        if (pre_integration->jacobian.maxCoeff() > 1e8 || pre_integration->jacobian.minCoeff() < -1e8) {
            printf("numerical unstable in preintegration");
        }
        Eigen::Matrix3d Ri_inv = Qi.inverse().toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();



        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
            jacobian_pose_i.setZero();
            jacobian_pose_i.block<3, 3>(O_P, O_P) = -Ri_inv;
            jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * newG * sum_dt * sum_dt + ((Pj - Pi) - Qj * Pbg) - Vi * sum_dt));
            Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
            jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();
            jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (newG * sum_dt + (Vj - Qj * (Utility::skewSymmetric(pre_integration->gyrj - Bgj) * Pbg)) - Vi));
            if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8) {
                printf("numerical unstable in preintegration");

            }
            jacobian_pose_i.block<3, 3>(O_P, O_V) = -Ri_inv * sum_dt;
            jacobian_pose_i.block<3, 3>(O_P, O_BA) = -dp_dba;
            jacobian_pose_i.block<3, 3>(O_P, O_BG) = -dp_dbg + Utility::skewSymmetric(Pbg) * sum_dt;
            jacobian_pose_i.block<3, 3>(O_R, O_BG) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration->delta_q).bottomRightCorner<3, 3>() * dq_dbg;
            jacobian_pose_i.block<3, 3>(O_V, O_V) = -Ri_inv;
            jacobian_pose_i.block<3, 3>(O_V, O_BA) = -dv_dba;
            jacobian_pose_i.block<3, 3>(O_V, O_BG) = -dv_dbg + Utility::skewSymmetric(Pbg);
            jacobian_pose_i.block<3, 3>(O_BA, O_BA) = -Eigen::Matrix3d::Identity();
            jacobian_pose_i.block<3, 3>(O_BG, O_BG) = -Eigen::Matrix3d::Identity();
            jacobian_pose_i = sqrt_info * jacobian_pose_i;

        }
        if (jacobians[1]) {
            Eigen::Map<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
            jacobian_pose_j.setZero();
            jacobian_pose_j.block<3, 3>(O_P, O_P) = Ri_inv;
            jacobian_pose_j.block<3, 3>(O_P, O_R) = Ri_inv * Rj * Utility::skewSymmetric(Pbg);
            Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
            jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();
            jacobian_pose_j.block<3, 3>(O_V, O_R) = Ri_inv * Rj * Utility::skewSymmetric(Utility::skewSymmetric(pre_integration->gyrj - Bgj) * Pbg);
            jacobian_pose_j.block<3, 3>(O_V, O_V) = Ri_inv;
            jacobian_pose_j.block<3, 3>(O_V, O_BG) = Ri_inv * Rj * (-Utility::skewSymmetric(Pbg));
            jacobian_pose_j.block<3, 3>(O_BA, O_BA) = Eigen::Matrix3d::Identity();
            jacobian_pose_j.block<3, 3>(O_BG, O_BG) = Eigen::Matrix3d::Identity();
            jacobian_pose_j = sqrt_info * jacobian_pose_j;

        }
        if (jacobians[2]) {
            Eigen::Map<Eigen::Matrix<double, 15, 3, Eigen::RowMajor>> jacobian_speedbias_t(jacobians[2]);
            jacobian_speedbias_t.setZero();
            jacobian_speedbias_t.block<3, 3>(O_P, O_P) = -Ri_inv * Rj + Eigen::Matrix<double, 3, 3>::Identity() + Utility::skewSymmetric(pre_integration->gyri - Bgi) * sum_dt;
            jacobian_speedbias_t.block<3, 3>(O_V, O_P) = -Ri_inv * Rj * Utility::skewSymmetric(pre_integration->gyrj - Bgj) + Utility::skewSymmetric(pre_integration->gyri - Bgi);
            jacobian_speedbias_t = sqrt_info * jacobian_speedbias_t;

        }
    }

    return true;
}





