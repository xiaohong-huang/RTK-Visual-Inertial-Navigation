#include"integration_base.h"
#include "../parameter/parameters.h"
#include "../utility/utility.h"

IntegrationBase::IntegrationBase(const Eigen::Vector3d& _acc_0, const Eigen::Vector3d& _gyr_0,
                                 const Eigen::Vector3d& _linearized_ba, const Eigen::Vector3d& _linearized_bg
                                 , const Eigen::Vector3d& _acc_scale, const Eigen::Vector3d& _gyr_scale
                                )
    : acc_0{_acc_0}, gyr_0{_gyr_0}, linearized_acc{_acc_0}, linearized_gyr{_gyr_0},
      linearized_ba{_linearized_ba}, linearized_bg{_linearized_bg},
      acc_scale {_acc_scale}, gyr_scale {_gyr_scale},
      jacobian {Eigen::Matrix<double, 15, 15>::Identity()}, covariance {Eigen::Matrix<double, 15, 15>::Zero()},
      sum_dt{0.0}, delta_p{Eigen::Vector3d::Zero()}, delta_q{Eigen::Quaterniond::Identity()}, delta_v{Eigen::Vector3d::Zero()},
      covariance_update{true}

{
    jacobian_s.setZero();
    noise = Eigen::Matrix<double, 18, 18>::Zero();
    noise.block<3, 3>(0, 0) =  (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
    noise.block<3, 3>(3, 3) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
    noise.block<3, 3>(6, 6) =  (ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
    noise.block<3, 3>(9, 9) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
    noise.block<3, 3>(12, 12) =  (ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
    noise.block<3, 3>(15, 15) =  (GYR_W * GYR_W) * Eigen::Matrix3d::Identity();
    gyri = gyr_0;
    gyrj = gyr_0;
    dt_i=0;dt_j=0;
}

void IntegrationBase::push_back(double dt, const Eigen::Vector3d& acc, const Eigen::Vector3d& gyr) {
    gyrj = gyr;
    dt_buf.push_back(dt);
    acc_buf.push_back(acc);
    gyr_buf.push_back(gyr);
    propagate(dt, acc, gyr);
}



void IntegrationBase::midPointIntegration(double _dt,
                                          const Eigen::Vector3d& _acc_0, const Eigen::Vector3d& _gyr_0,
                                          const Eigen::Vector3d& _acc_1, const Eigen::Vector3d& _gyr_1,
                                          const Eigen::Vector3d& delta_p, const Eigen::Quaterniond& delta_q, const Eigen::Vector3d& delta_v,
                                          const Eigen::Vector3d& linearized_ba, const Eigen::Vector3d& linearized_bg,
                                          Eigen::Vector3d& result_delta_p, Eigen::Quaterniond& result_delta_q, Eigen::Vector3d& result_delta_v,
                                          Eigen::Vector3d& result_linearized_ba, Eigen::Vector3d& result_linearized_bg, bool update_jacobian) {


    Eigen::Vector3d un_acc_0 = delta_q * ((_acc_0.array() * acc_scale.array()).matrix() - linearized_ba);
    Eigen::Vector3d un_gyr = 0.5 * ((_gyr_0 + _gyr_1).array() * gyr_scale.array()).matrix() - linearized_bg;

    result_delta_q = delta_q * Eigen::Quaterniond(1, un_gyr(0) * _dt / 2, un_gyr(1) * _dt / 2, un_gyr(2) * _dt / 2);

    Eigen::Vector3d un_acc_1 = result_delta_q * ((_acc_1.array() * acc_scale.array()).matrix() - linearized_ba);

    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    result_delta_p = delta_p + delta_v * _dt + 0.5 * un_acc * _dt * _dt;
    result_delta_v = delta_v + un_acc * _dt;
    result_linearized_ba = linearized_ba;
    result_linearized_bg = linearized_bg;

    if (update_jacobian) {
        Eigen::Vector3d w_x = 0.5 * ((_gyr_0 + _gyr_1).array() * gyr_scale.array()).matrix() - linearized_bg;
        Eigen::Vector3d a_0_x = (_acc_0.array() * acc_scale.array()).matrix() - linearized_ba;
        Eigen::Vector3d a_1_x = (_acc_1.array() * acc_scale.array()).matrix() - linearized_ba;

        Eigen::Matrix3d R_w_x, R_a_0_x, R_a_1_x;

        R_w_x << 0, -w_x(2), w_x(1),
              w_x(2), 0, -w_x(0),
              -w_x(1), w_x(0), 0;
        R_a_0_x << 0, -a_0_x(2), a_0_x(1),
                a_0_x(2), 0, -a_0_x(0),
                -a_0_x(1), a_0_x(0), 0;
        R_a_1_x << 0, -a_1_x(2), a_1_x(1),
                a_1_x(2), 0, -a_1_x(0),
                -a_1_x(1), a_1_x(0), 0;

        Eigen::MatrixXd F = Eigen::MatrixXd::Zero(15, 15);
        F.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
        F.block<3, 3>(0, 3) = -0.25 * delta_q.toRotationMatrix() * R_a_0_x * _dt * _dt +
                              -0.25 * result_delta_q.toRotationMatrix() * R_a_1_x * (Eigen::Matrix3d::Identity() - R_w_x * _dt) * _dt * _dt;
        F.block<3, 3>(0, 6) = Eigen::MatrixXd::Identity(3, 3) * _dt;
        F.block<3, 3>(0, 9) = -0.25 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * _dt * _dt;
        F.block<3, 3>(0, 12) = -0.25 * result_delta_q.toRotationMatrix() * R_a_1_x * _dt * _dt * -_dt;

        F.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() - R_w_x * _dt;
        F.block<3, 3>(3, 12) = -1.0 * Eigen::MatrixXd::Identity(3, 3) * _dt;

        F.block<3, 3>(6, 3) = -0.5 * delta_q.toRotationMatrix() * R_a_0_x * _dt +
                              -0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * (Eigen::Matrix3d::Identity() - R_w_x * _dt) * _dt;
        F.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity();
        F.block<3, 3>(6, 9) = -0.5 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * _dt;
        F.block<3, 3>(6, 12) = -0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * _dt * -_dt;

        F.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity();

        F.block<3, 3>(12, 12) = Eigen::Matrix3d::Identity();



        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(15, 18);
        V.block<3, 3>(0, 0) =  0.25 * delta_q.toRotationMatrix() * _dt * _dt;
        V.block<3, 3>(0, 3) =  0.25 * -result_delta_q.toRotationMatrix() * R_a_1_x  * _dt * _dt * 0.5 * _dt;
        V.block<3, 3>(0, 6) =  0.25 * result_delta_q.toRotationMatrix() * _dt * _dt;
        V.block<3, 3>(0, 9) =  V.block<3, 3>(0, 3);

        V.block<3, 3>(3, 3) =  0.5 * Eigen::MatrixXd::Identity(3, 3) * _dt;
        V.block<3, 3>(3, 9) =  0.5 * Eigen::MatrixXd::Identity(3, 3) * _dt;

        V.block<3, 3>(6, 0) =  0.5 * delta_q.toRotationMatrix() * _dt;
        V.block<3, 3>(6, 3) =  0.5 * -result_delta_q.toRotationMatrix() * R_a_1_x  * _dt * 0.5 * _dt;
        V.block<3, 3>(6, 6) =  0.5 * result_delta_q.toRotationMatrix() * _dt;
        V.block<3, 3>(6, 9) =  V.block<3, 3>(6, 3);

        V.block<3, 3>(9, 12) = Eigen::MatrixXd::Identity(3, 3) * _dt;
        V.block<3, 3>(12, 15) = Eigen::MatrixXd::Identity(3, 3) * _dt;

        noise.block<3, 3>(3, 3) +=  (GYR_NONLINEARITY * GYR_NONLINEARITY * _gyr_0.cwiseAbs2() ).asDiagonal();
        noise.block<3, 3>(9, 9) +=  (GYR_NONLINEARITY * GYR_NONLINEARITY * _gyr_1.cwiseAbs2() ).asDiagonal();

        Eigen::Matrix<double, 15, 6> jacobian_sk = Eigen::Matrix<double, 15, 6>::Zero();

        jacobian_sk.block<3, 3>(0, 0) = 0.25 * (delta_q.toRotationMatrix() * _acc_0.asDiagonal() + result_delta_q.toRotationMatrix() * _acc_1.asDiagonal()) * _dt * _dt;
        jacobian_sk.block<3, 3>(6, 0) = 0.5 * (delta_q.toRotationMatrix() * _acc_0.asDiagonal() + result_delta_q.toRotationMatrix() * _acc_1.asDiagonal()) * _dt;
        jacobian_sk.block<3, 3>(0, 3) = -F.block<3, 3>(0, 12) * (0.5 * (_gyr_0 + _gyr_1)).asDiagonal();
        jacobian_sk.block<3, 3>(3, 3) = -F.block<3, 3>(3, 12) * (0.5 * (_gyr_0 + _gyr_1)).asDiagonal();
        jacobian_sk.block<3, 3>(6, 3) = -F.block<3, 3>(6, 12) * (0.5 * (_gyr_0 + _gyr_1)).asDiagonal();
        jacobian_s = F * jacobian_s + jacobian_sk;


        jacobian = F * jacobian;
        covariance = F * covariance * F.transpose() + V * noise * V.transpose();
        covariance.block<3, 3>(0, 0) += Eigen::Matrix3d::Identity() * INT_N * INT_N * dt * dt;

        noise.block<3, 3>(3, 3) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
        noise.block<3, 3>(9, 9) =  (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();


        covariance_update = true;
    }
}



Eigen::Matrix<double, 15, 15> IntegrationBase::get_sqrtinfo() {

    if (covariance_update) {
        covariance_update = false;
        sqrt_info = Eigen::LLT<Eigen::Matrix<double, 15, 15>>(covariance.inverse()).matrixL().transpose();
    }
    return sqrt_info;

}

void IntegrationBase::propagate(double _dt, const Eigen::Vector3d& _acc_1, const Eigen::Vector3d& _gyr_1) {
    dt = _dt;
    acc_1 = _acc_1;
    gyr_1 = _gyr_1;
    Eigen::Vector3d result_delta_p;
    Eigen::Quaterniond result_delta_q;
    Eigen::Vector3d result_delta_v;
    Eigen::Vector3d result_linearized_ba;
    Eigen::Vector3d result_linearized_bg;
    midPointIntegration(_dt, acc_0, gyr_0, _acc_1, _gyr_1, delta_p, delta_q, delta_v,
                        linearized_ba, linearized_bg,
                        result_delta_p, result_delta_q, result_delta_v,
                        result_linearized_ba, result_linearized_bg, 1);

    delta_p = result_delta_p;
    delta_q = result_delta_q;
    delta_v = result_delta_v;
    linearized_ba = result_linearized_ba;
    linearized_bg = result_linearized_bg;
    delta_q.normalize();
    sum_dt += dt;
    acc_0 = acc_1;
    gyr_0 = gyr_1;

}


void IntegrationBase::repropagate(const Eigen::Vector3d& _linearized_ba, const Eigen::Vector3d& _linearized_bg
                                  , const Eigen::Vector3d& _acc_scale, const Eigen::Vector3d& _gyr_scale
                                 ) {
    sum_dt = 0.0;
    acc_0 = linearized_acc;
    gyr_0 = linearized_gyr;
    delta_p.setZero();
    delta_q.setIdentity();
    delta_v.setZero();
    linearized_ba = _linearized_ba;
    linearized_bg = _linearized_bg;
    jacobian.setIdentity();
    covariance.setZero();
    jacobian_s.setZero();
    acc_scale = _acc_scale;
    gyr_scale = _gyr_scale;

    covariance_update = true;
    for (int i = 0; i < static_cast<int>(dt_buf.size()); i++)
        propagate(dt_buf[i], acc_buf[i], gyr_buf[i]);
}
