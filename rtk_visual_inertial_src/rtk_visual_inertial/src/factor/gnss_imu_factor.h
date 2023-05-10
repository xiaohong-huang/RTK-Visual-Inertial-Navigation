
#pragma once
#include <eigen3/Eigen/Dense>
#include "integration_base.h"
#include"imu_factor.h"
#include "marginalization_factor.h"

enum HessianOrder {
    O_Pose1 = 0,
    O_Pose2,
    O_N,
    O_Pose0,
    O_Size

};



class IMUGNSSBase {
  public:


    IMUGNSSBase(double* para_Pose0_, double* para_SpeedBias0_, ceres::Problem* Myproblem_):
        para_pose0(para_Pose0_), para_speed_bias0(para_SpeedBias0_), my_problem(Myproblem_) {};
    ~IMUGNSSBase();
    void Init();
    void Clear();
    void InitHessianRhs(double** hessian55, int* hessian_size, double** rhs5);
    void SetLastImuFactor(IntegrationBase* pre_integration, double* para_Pose2, double* para_SpeedBias2);
    void AddMargInfo(MarginalizationInfo* marg_info_gnss, IntegrationBase* pre_integration, double* para_Poses, double* para_speed_bias);
    void JacobianResidualUpdateHessianRhs(std::vector<int>index, std::vector<int>mapsize, double** jacobian, double** hessian55, double* residual, double** rhs5);
    void HessianUpdateHessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2);
    void RhsUpdateRhs(double** rhs_desc, double* rhs_src, int block1, int size1);
    void CopyHessian2Hessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2);
    void CopyRhs2Rhs(double** rhs_desc, double* rhs_src, int block1, int size1);
    void MargPose1(double** hessian55, int* hessian_size, double** rhs5);
    void MoveHessianData(int index);
    void UpdateSchurComponent();
    void UpdateJacobResidual(double** jacobians, double* residuals);
    void UpdateRhsN();
    void UpdateRhsPose(int i);
    void UpdateDeltaValues(Eigen::Matrix<double, 7, 1, Eigen::ColMajor>& Pi, Eigen::Matrix<double, 9, 1, Eigen::ColMajor>& Bi, Eigen::Matrix<double, 7, 1, Eigen::ColMajor>& Pj, Eigen::Matrix<double, 9, 1, Eigen::ColMajor>& Bj, Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>& phase_biases_values);
    void UpdateHiddenState();
    void ResetMem();
    void UpdateTrustRegion();
    void SaveLastStates(Eigen::Matrix<double, 7, 1, Eigen::ColMajor>& Pi, Eigen::Matrix<double, 9, 1, Eigen::ColMajor>& Bi, Eigen::Matrix<double, 7, 1, Eigen::ColMajor>& Pj, Eigen::Matrix<double, 9, 1, Eigen::ColMajor>& Bj);
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) ;
    void AddMidMargInfo(MarginalizationInfo* marginfo_);
    Eigen::Matrix<double, 15, 1, Eigen::ColMajor> GetInc(int i);

    //pointer to the phase biases
    std::vector<double*>gnss_phase_biases;
    //pointer to the phase biases
    std::set<double*>gnss_phase_biases_set;
    //pointer to the velocity, gyr bias and accerametor bias
    std::vector<double*>gnss_speed_bias;
    //linerization points of the velocity, gyr bias and accerametor bias
    std::vector<double*>gnss_speed_bias_lin;
    //pointer to the position and orientation
    std::vector<double*>gnss_poses;
    //linerization points of the position and orientation
    std::vector<double*>gnss_poses_lin;
    //phase biases values
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> phase_biases_values;
    //hessian matrix related to the phase biases
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> phase_biases_hessians;
    //right hand side vector related to the phase biases.
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> phase_biases_rhs;
    //hessian matrix related to the GNSS frame(position, orientation, velocity, gyr bias and acc bias)
    std::vector<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> pose_hessians;
    //hessian matrix related to two consecutive GNSS frame, not zero only when the middle marginalization is performed.
    Eigen::Matrix<double, 15, 15, Eigen::RowMajor> pose1_pose2_hessians;
    //hessian matrix related to the GNSS frrame and phase biases.
    std::vector<Eigen::Matrix<double, 15, Eigen::Dynamic, Eigen::RowMajor>> pose_phase_biases_hessians;
    //right hand side vector related to the GNSS frame.
    std::vector<Eigen::Matrix<double, 15, 1, Eigen::ColMajor>> pose_rhses;
    //incremental vector related to the outside states (phase biases+states of the consecutive visual frames).
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> INC;
    //imu factor remained in the current GNSS-IMU factor
    std::vector<IMUFactor*>imu_factors;
    //the last imu factor related the newest gnss frame and his consecutive visual frame
    IMUFactor* last_imu_factor = 0;
    //store the hessian matirx related to the marginalized states and its related states.
    std::vector<double*>hmn_save[O_Size];
    //store the right hand side vector related to the marginalized states and its related states.
    std::vector<double*>rhsmn_save;
    //size of the hessian matrix related to each component
    int hessian_size[O_Size] = {-1};
    //hessian matrix related to the current elimination step
    double* hessian55[O_Size * O_Size] = {0};
    //right hand side vector related to the current elimination step
    double* rhs5[O_Size] = {0};
    //incremental vector related to the back-ward operation
    double* delta5[O_Size] = {0};
    //states estimation in last iteration.
    Eigen::Matrix<double, 7, 1, Eigen::ColMajor> Pi_old;
    //states estimation in last iteration.
    Eigen::Matrix<double, 9, 1, Eigen::ColMajor>  Bi_old;
    //states estimation in last iteration.
    Eigen::Matrix<double, 7, 1, Eigen::ColMajor> Pj_old;
    //states estimation in last iteration.
    Eigen::Matrix<double, 9, 1, Eigen::ColMajor>  Bj_old;
    //states estimation in last iteration.
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> phase_biases_values_old;
    //1 if the Jacobian updation required
    bool update_flag;
    //1 if it's the first iteration
    bool history_flag;
    //jacobian related to the outside states (phase biases+states of the consecutive visual frames)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> schur_jacobian;
    //residual related to the outside states (phase biases+states of the consecutive visual frames)
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> schur_residual;
    //index states in the hessian matrix
    int hessian_index[O_Size - 1];
    //size of the outside states.
    int global_size = 0;
    //number of round for updating the Jacobian
    int update_round = 0;
    //number of round for updating the residual
    int evaluate_round = 0;
    //residual block id of the current GNSS-IMU factor in Ceres Solver
    ceres::ResidualBlockId residualBlockId = 0;
    //pointer to the first position and orientation
    double* para_pose0;
    //pointer to the first velocity, gyr bias and acc bias
    double* para_speed_bias0;
    //outsite states
    std::vector<double*>param;
    //global ceres problem;avaliable only when USE_GLOBAL_OPTIMIZATION is 1
    ceres::Problem* my_problem;
    //threshold for re-eliminating. not use for current version
    double update_threshold = 1e-3;
    //middle marginalization factor
    MarginalizationInfo* gnss_middle_marginfo = 0;
    //number of the avaliable gnss frames//need to fixed
    int gnss_Index = 0;
    //hessian matrix related to the outside states (phase biases+states of the consecutive visual frames)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> hessian_dense;
    //right hand side vector related to the outside states (phase biases+states of the consecutive visual frames)
    Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> rhs_dense;

};


class IMUGNSSFactor: public ceres::CostFunction {
  public:
    IMUGNSSFactor(IMUGNSSBase* IMUGNSS_info_);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)const;

    IMUGNSSBase* IMUGNSS_info;
};
