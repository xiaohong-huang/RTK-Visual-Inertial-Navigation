
#pragma once
#include <eigen3/Eigen/Dense>
#include "../factor/integration_base.h"
#include"../factor/imu_factor.h"
#include "../factor/marginalization_factor.h"
#include "../utility/utility.h"
#include "../factor/sparse_matrix.h"

#define BTOC_UPDATE 1

enum HessianOrder {
    O_Pose1 = 0,
    O_Pose2,
    O_N,
    O_Pose0,
    O_EXTRINSIC,
    O_Size

};
#if BTOC_UPDATE
class PJacobian {
  public:
    inline  PJacobian(const Eigen::Matrix3d& R_WI_WC_, const Eigen::Vector3d& P_WI_WC_, const double& scale_factor_, const Eigen::Matrix3d& Rbi, const Eigen::Vector3d& Pgi, const Eigen::Matrix3d& ric_, const Eigen::Vector3d tic_, const Eigen::Vector3d& pbg_) {
        dP0_dP0.setZero();
        dP0_dscale.setZero();
        dP0_dR_WI_WC.setZero();
        dP0_dextrinsic.setZero();
        dP0_dpbg.setZero();
        r_wi_wc = R_WI_WC_;
        p_wi_wc = P_WI_WC_;
        scale_factor = scale_factor_;
        ric = ric_;
        tic = tic_;
        pbg = pbg_;
        Eigen::Matrix3d Rci = r_wi_wc.transpose() * Rbi * ric;
        Eigen::Vector3d Pci = (r_wi_wc.transpose() * (Pgi - p_wi_wc) +  r_wi_wc.transpose() * Rbi *  (tic - pbg)) / scale_factor;

        ASSERT((r_wi_wc * Pci * scale_factor + p_wi_wc - r_wi_wc * Rci * ric.transpose() * (tic - pbg) - Pgi).norm() < 1e-4);
        ASSERT((r_wi_wc * Rci * ric.transpose() - Rbi).norm() < 1e-4);

        dP0_dP0.block<3, 3>(0, 0) = r_wi_wc.transpose() / scale_factor;
        dP0_dP0.block<3, 3>(0, 3) = r_wi_wc.transpose() * Rbi * -Utility::skewSymmetric( tic - pbg) / scale_factor;
        dP0_dP0.block<3, 3>(3, 3) = ric.transpose();
        if (ESTIMATE_SCALE)
            dP0_dscale.segment(0, 3) = -Pci / scale_factor;
        if (ESTIMATE_P_WI_WC)
            dP0_dR_WI_WC.block<3, 3>(0, 0) = -r_wi_wc.transpose() / scale_factor;
        if (ESTIMATE_R_WI_WC) {
            dP0_dR_WI_WC.block<3, 3>(0, 3) = Utility::skewSymmetric(Pci);
            dP0_dR_WI_WC.block<3, 3>(3, 3) = -Rci.transpose();
        }
        if (ESTIMATE_EXTRINSIC) {
            dP0_dextrinsic.block<3, 3>(0, 0) = r_wi_wc.transpose() * Rbi / scale_factor;
            dP0_dextrinsic.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
        }
        if (ESTIMATE_PBG)
            dP0_dpbg.block<3, 3>(0, 0) = -r_wi_wc.transpose() * Rbi / scale_factor;

    }
    flower_focus::Matrix6d dP0_dP0;
    flower_focus::Vector6d dP0_dscale;
    flower_focus::Matrix6d dP0_dR_WI_WC;
    flower_focus::Matrix6d dP0_dextrinsic;
    flower_focus::Matrix6_3d dP0_dpbg;

    Eigen::Matrix3d r_wi_wc;
    Eigen::Vector3d p_wi_wc;
    double scale_factor;
    Eigen::Matrix3d ric;
    Eigen::Vector3d tic;
    Eigen::Vector3d pbg;

};
#endif


class IMUPreFactor {
  public:
    IMUPreFactor() = delete;
    IMUPreFactor(IMUFactor* _factor) {
        factor = _factor;
    }
    ~IMUPreFactor() {
        delete factor;
    }
    std::vector<double*>parameters;
    IMUFactor* factor = 0;
    std::vector<flower_focus::MatrixXd > imu_jacobians_in;

    flower_focus::Vector15d  imu_residual_in;
    flower_focus::Vector15d  imu_residual_in_save;
    flower_focus::Vector15d  model_residual_accum;
    flower_focus::Vector15d  model_residual_accum_save;
    flower_focus::Vector15d  imu_residual_in_linerized;
    double* imu_jacobians_raw_in[10] = {0};

};


class IMUGNSSBase {
  public:


    IMUGNSSBase() {

    };
    ~IMUGNSSBase();
    void Init();
    void Clear();
    void InitHessianRhs(double** hessian55, int* hessian_size, double** rhs5);
    void AddImuFactor(IntegrationBase* pre_integration);
    void AddMargInfo(MarginalizationInfo* marg_info_gnss, double* para_Poses, double* para_speed_bias);
    void JacobianResidualUpdateHessianRhs(std::vector<int>index, std::vector<int>mapsize, double** jacobian, double** hessian55, double* residual, double** rhs5, double** gradient5);
    void HessianUpdateHessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2);
    void RhsUpdateRhs(double** rhs_desc, double* rhs_src, int block1, int size1);
    void CopyHessian2Hessian(double** hessian_desc, double* hessian_src, int block1, int block2, int size1, int size2);
    void CopyRhs2Rhs(double** rhs_desc, double* rhs_src, int block1, int size1);
    void MargPose1(double** hessian55, int* hessian_size, double** rhs5);
    void MargPose1rhsOnly(double** hessian55, int* hessian_size, double** rhs5, int index);
    void MoveHessianData(int index);
    void MoveRhsData(int index);
    void UpdateSchurComponent(double* rhs, double* lhs, double* gradient);
    void UpdateRhsN();
    void UpdateRhsPose(int i);
    void UpdateHiddenGaussStep(double*, double*, double*);
    void ResetMem();
    void UpdateTrustRegion(double mu);
    bool EvaluateLhsRhs( double* rhs, double* lhs, double* gradient, double mu);
    bool EvaluateCost(double* cost, double* model_cost_change);
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
    flower_focus::VectorXd phase_biases_values;
    //hessian matrix related to the phase biases
    flower_focus::MatrixXd phase_biases_hessians;
    //right hand side vector related to the phase biases.
    flower_focus::VectorXd phase_biases_rhs;
    //hessian matrix related to the GNSS frame(position, orientation, velocity, gyr bias and acc bias)
    std::vector<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> pose_hessians;
    //hessian matrix related to two consecutive GNSS frame, not zero only when the middle marginalization is performed.
    Eigen::Matrix<double, 15, 15, Eigen::RowMajor> pose1_pose2_hessians;
    //hessian matrix related to the GNSS frrame and phase biases.
    std::vector<Eigen::Matrix<double, 15, Eigen::Dynamic, Eigen::RowMajor>> pose_phase_biases_hessians;
    //right hand side vector related to the GNSS frame.
    std::vector<Eigen::Matrix<double, 15, 1, Eigen::ColMajor>> pose_rhses;
    //store the hessian matirx related to the marginalized states and its related states.
    std::vector<double*>hmn_save[O_Size];
    //store the right hand side vector related to the marginalized states and its related states.
    std::vector<double*>rhsmn_save;
    std::vector<double*>gauss_save;
    std::vector<double*>gradient_save;
    //size of the hessian matrix related to each component
    int hessian_size[O_Size] = {-1};
    //hessian matrix related to the current elimination step
    double* hessian55[O_Size * O_Size] = {0};
    //right hand side vector related to the current elimination step
    double* rhs5[O_Size] = {0};
    double* gradient5[O_Size] = {0};
    //incremental vector related to the back-ward operation
    double* gauss5[O_Size] = {0};
    double* inc5[O_Size] = {0};
    std::vector<double*>inc_save;
    //1 if the Jacobian updation required
    bool update_flag;
    //1 if it's the first iteration
    bool history_flag;
    //jacobian related to the outside states (phase biases+states of the consecutive visual frames)
    flower_focus::MatrixXd schur_jacobian;
    //residual related to the outside states (phase biases+states of the consecutive visual frames)
    flower_focus::VectorXd schur_residual;
    //index states in the hessian matrix
    int hessian_index[O_Size - 1];
    //size of the outside states.
    int global_size = 0;
    //number of round for updating the Jacobian
    int update_round = 0;
    //number of round for updating the residual
    int evaluate_round = 0;
    //threshold for re-eliminating. not use for current version
    double update_threshold = 1e-3;
    //middle marginalization factor
    MarginalizationInfo* gnss_middle_marginfo = 0;
    //number of the avaliable gnss frames//need to fixed
    int gnss_Index = 0;
    //hessian matrix related to the outside states (phase biases+states of the consecutive visual frames)
    flower_focus::MatrixXd hessian_dense;
    //right hand side vector related to the outside states (phase biases+states of the consecutive visual frames)
    flower_focus::VectorXd rhs_dense;
    flower_focus::VectorXd gradient_dense;

    double gnss_init_cost = 0;

    flower_focus::Vector3d Pgi;
    flower_focus::Vector3d Pci;
    flower_focus::Matrix3d Rbi;
    Eigen::Quaterniond Qbi;
    flower_focus::Vector9d Bi;
    Eigen::Matrix3d Rci;
    Eigen::Matrix3d r_wi_wc;
    Eigen::Vector3d p_wi_wc;
    flower_focus::Vector7d PgQbi;

    flower_focus::Vector3d Pgj;
    flower_focus::Vector3d Pcj;
    flower_focus::Matrix3d Rbj;
    Eigen::Quaterniond Qbj;
    flower_focus::Vector9d Bj;
    Eigen::Matrix3d Rcj;
    flower_focus::Vector7d PgQbj;

    Eigen::Vector3d tic;
    Eigen::Matrix3d ric;

    Eigen::VectorXd extrinsic;

    double scale_factor;
    void Reparameterization(double const* const* parameters);
    void ReparameterizationNew(double const* const* parameters);

    int O_POSi;
    int O_BIASi;
    int O_POSj;
    int O_BIASj;
    int O_TD;
    int O_ACC_S;
    int O_PBG;
    int O_RWGW;
    int O_PHASEBIAS;
    int O_SCALE;
    int O_POSGLOBAL;
    int O_EXTRINSIC_IC;
    int O_FULL;
    std::vector<int>order2p_local;

    void EvaluateGaussAndGradientStep(const double* gradient, const double* gauss, double* gauss_squared_norm, double* gradient_squared_norm, double* gradient_dot_gauss, double* alpha);

    void EvaluateAlpha(double* alpha);
    std::vector<MatirxInfo>matrix_info;
    std::vector<IMUPreFactor*>imu_factors;
    void UpdateInsideStateUseGaussAndGradient(double const* const* parameters, double gradient_scale, double gauss_scale);
    double gnss_cost;
    double last_gnss_cost;
    void SaveLoadVariable(bool);
    std::vector<flower_focus::Vector7d>gnss_poses_old;
    std::vector<flower_focus::Vector9d>gnss_speed_bias_old;
    bool EvaluateRhs(double* rhs, double* gradient);
    void BuildMatrixInfo();
    double nonlinear_quantity;
    bool evaluate_rhs_only;
    Eigen::Vector3d pbg;
#if BTOC_UPDATE
    flower_focus::Vector3d pbg_gauss;
    flower_focus::Vector1d scale_gauss;
    flower_focus::Vector6d posgloabl_gauss;
    flower_focus::Vector6d extrinsic_gauss;

    flower_focus::Vector3d pbg_gradient;
    flower_focus::Vector1d scale_gradient;
    flower_focus::Vector6d posgloabl_gradient;
    flower_focus::Vector6d extrinsic_gradient;

    flower_focus::Vector3d pbg_inc;
    flower_focus::Vector1d scale_inc;
    flower_focus::Vector6d posgloabl_inc;
    flower_focus::Vector6d extrinsic_inc;


    double new_scale;
    flower_focus::Matrix3d new_r_wi_wc;
    flower_focus::Vector3d new_p_wi_wc;
    Eigen::Vector3d new_tic;
    Eigen::Matrix3d new_ric;
    flower_focus::Vector3d new_pbg;
    std::vector<PJacobian>p_jacobians;
#endif




    void EvaluateNonlinearQuantity() {
        if (history_flag) {
            nonlinear_quantity = 0;
            for (int i = 0; i < (int)imu_factors.size(); i++)
                nonlinear_quantity += (imu_factors[i]->imu_residual_in + imu_factors[i]->model_residual_accum - imu_factors[i]->imu_residual_in_linerized).squaredNorm();

            ASSERT(!isnan(nonlinear_quantity));
        } else
            nonlinear_quantity = NLQ_THRESHOLD + 100;

    }

};


class IMUGNSSFactor: public CostFunction {
  public:
    IMUGNSSFactor(IMUGNSSBase* IMUGNSS_info_);
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)const;


    bool EvaluateLhsRhs(double const* const* parameters, double* rhs, double* lhs, double* gradient, double mu) {
        IMUGNSS_info->Reparameterization(parameters);
        IMUGNSS_info->EvaluateLhsRhs(rhs, lhs, gradient, mu);
        return true;
    }

    bool EvaluateRhs(double const* const* parameters, double* rhs, double* gradient) {
        IMUGNSS_info->Reparameterization(parameters);
        IMUGNSS_info->EvaluateRhs(rhs, gradient);
        return true;
    }

    bool EvaluateCost(double const* const* parameters, double* cost) {
        IMUGNSS_info->Reparameterization(parameters);
        IMUGNSS_info->EvaluateCost(cost, 0);
        return true;
    }
    bool SaveOrRestoreHidenStates(bool is_save) {
        IMUGNSS_info->SaveLoadVariable(is_save);
        if (is_save)IMUGNSS_info->EvaluateNonlinearQuantity();
        return true;
    }

    bool  UpdateInsideGaussStep(double* gradient, double* gauss,
                                double* gauss_newton_step_inside_square_norm, double* alpha1, double* gradient_squared_norm, double* gradient_dot_gauss_newton_inside)  {
        IMUGNSS_info->EvaluateGaussAndGradientStep(gradient, gauss, gauss_newton_step_inside_square_norm, gradient_squared_norm, gradient_dot_gauss_newton_inside, alpha1);
        return true;
    }

    bool  UpdateInsideStateUseGaussAndGradient(double const* const* parameters, double gradient_scale, double gauss_scale)  {
#if BTOC_UPDATE
        IMUGNSS_info->ReparameterizationNew(parameters);
#endif
        IMUGNSS_info->UpdateInsideStateUseGaussAndGradient(parameters, gradient_scale, gauss_scale);
        return true;
    }

    bool  EvaluateModelCostChange(double* model_cost_change)  {
        IMUGNSS_info->EvaluateCost(0, model_cost_change);
        return true;
    }

    IMUGNSSBase* IMUGNSS_info;
};
