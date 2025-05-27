

#pragma once

#include <Eigen/Dense>
#include "../parameter/parameters.h"
class IntegrationBase {
  public:
    IntegrationBase() = delete;
    IntegrationBase(const Eigen::Vector3d& _acc_0, const Eigen::Vector3d& _gyr_0,
                    const Eigen::Vector3d& _linearized_ba, const Eigen::Vector3d& _linearized_bg
                    , const Eigen::Vector3d& _acc_scale, const Eigen::Vector3d& _gyr_scale );

    void push_back(double dt, const Eigen::Vector3d& acc, const Eigen::Vector3d& gyr);

    void midPointIntegration(double _dt,
                             const Eigen::Vector3d& _acc_0, const Eigen::Vector3d& _gyr_0,
                             const Eigen::Vector3d& _acc_1, const Eigen::Vector3d& _gyr_1,
                             const Eigen::Vector3d& delta_p, const Eigen::Quaterniond& delta_q, const Eigen::Vector3d& delta_v,
                             const Eigen::Vector3d& linearized_ba, const Eigen::Vector3d& linearized_bg,
                             Eigen::Vector3d& result_delta_p, Eigen::Quaterniond& result_delta_q, Eigen::Vector3d& result_delta_v,
                             Eigen::Vector3d& result_linearized_ba, Eigen::Vector3d& result_linearized_bg, bool update_jacobian);

    void propagate(double _dt, const Eigen::Vector3d& _acc_1, const Eigen::Vector3d& _gyr_1);

    Eigen::Matrix<double, 15, 15> get_sqrtinfo();
    void repropagate(const Eigen::Vector3d& _linearized_ba, const Eigen::Vector3d& _linearized_bg
                     , const Eigen::Vector3d& _acc_scale, const Eigen::Vector3d& _gyr_scale);
                     
    double dt;
    Eigen::Vector3d acc_0, gyr_0;
    Eigen::Vector3d acc_1, gyr_1;


    const Eigen::Vector3d linearized_acc, linearized_gyr;
    Eigen::Vector3d linearized_ba, linearized_bg;

    Eigen::Vector3d acc_scale;
    Eigen::Vector3d gyr_scale;
    
    std::vector<double> dt_buf;
    std::vector<Eigen::Vector3d> acc_buf;
    std::vector<Eigen::Vector3d> gyr_buf;

    Eigen::Matrix<double, 15, 15> jacobian, covariance;
    Eigen::Matrix<double, 15, 15> step_jacobian;
    Eigen::Matrix<double, 15, 18> step_V;
    Eigen::Matrix<double, 18, 18> noise;

    double sum_dt;
    Eigen::Vector3d delta_p;
    Eigen::Quaterniond delta_q;
    Eigen::Vector3d delta_v;

    Eigen::Matrix<double, 15, 15> sqrt_info;
    bool covariance_update;

    Eigen::Vector3d gyri, gyrj;
    Eigen::Matrix<double, 15, 6> jacobian_s;
    double dt_i,dt_j;


};
