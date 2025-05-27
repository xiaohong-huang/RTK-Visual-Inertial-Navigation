
#pragma once
#include "../parameter/parameters.h"
#include <eigen3/Eigen/Dense>
#include "integration_base.h"

#include"imu_factor.h"
#include "marginalization_factor.h"
#include "loss_function.h"
#include <omp.h>
#include "sparse_matrix.h"
#include "../gnss/gnss_imu_factor.h"

#if ESTIMATE_EXTRINSIC2&&USE_STEREO
#define VRNUM (6+6+1+1+6)
#elif USE_STEREO&&!ESTIMATE_EXTRINSIC2
#define VRNUM (6+6+1+1)
#else
#define VRNUM (6+6+1)
#endif


class IMUGNSSPreFactor {
  public:
    IMUGNSSPreFactor() = delete;
    IMUGNSSPreFactor(IMUGNSSFactor* _factor) {
        factor = _factor;
    }
    ~IMUGNSSPreFactor() {
        delete factor;
    }
    std::vector<double*>parameters;
    int pindexi, pindexj;
    IMUGNSSFactor* factor = 0;
    std::vector<int>evaluate_idx;
    std::vector<int>orders;
    std::vector<int>order_idx;
    std::vector<int>parameter_sizes;
    int parameter_size_sum;
    std::vector<double*>phase_biases;

};


class VisualFactor {
  public:

    std::set<int>pindexs_set;
    std::vector<int> pindexs;
    std::vector<Eigen::Vector3d> ptss;
    double* lhs_pos_idepth;
    double rhs_idepth[1];
    double lhs_idepth[1];
    double rhs_idepth_save[1];
    double* p_idepth;
    int min_pindex;
#if USE_STEREO
    std::vector<bool>is_stereo;
    double lhs_scale_idepth=0;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    double lhs_extrinsic_idepth[6];
#endif
#endif
};

class VisualInertialBase  {
  public:
    VisualInertialBase() {
        cauchy_loss_function = new CauchyLoss(1.0);
    };
    ~VisualInertialBase();
    void AddVisualFactorShort(int feature_id, double*, double*, double*, int, int, Eigen::Vector3d, Eigen::Vector3d
#if USE_STEREO
                              , bool is_stereo_
#endif
                             );
    void AddVisualFactorLong(int feature_id, double*, double*, double*, int, int, Eigen::Vector3d, Eigen::Vector3d
#if USE_STEREO
                             , bool is_stereo_
#endif
                            );
    void AddIMUFactor(IMUGNSSPreFactor* imu_pre_factor);

    void ConstructBiasConnections();
    void ConstructOrderingMap();
    void ResetMem();
    void PrepareMatrix();
    void FactorUpdates();

    void VisualJacobianResidualUpdatelhsRhs(int pindexi, int pindexj,
                                            Eigen::Vector3d pts_i, Eigen::Vector3d pts_j,
                                            double* lhs_pos_idepth,
                                            double* rhs_idepth,
                                            double* lhs_idepth,
                                            LossFunction* loss_function0,
                                            double* jacobian_pointer,
                                            double* visual_residual_raw_in, double visual_cost,
                                            int min_pindex,
                                            double* p_idepth
#if USE_STEREO
                                            , double* lhs_scale_idepth, bool is_stereo
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                                            , double* lhs_extrinsic_idepth
#endif
#endif
                                           );

    void EvaluateIdepthsShort();
    void EvaluateIdepthsLong();

    void IMUJacobianResidualUpdatelhsRhs();

    void MargIdepth(const std::set<int>& pindexs, double* lhs_pos_idepth, double* rhs_idepth, double* lhs_idepth, int
#if USE_STEREO
                    , double* lhs_scale_idepth
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                    , double* lhs_extrinsic_idepth
#endif
#endif
                   );
    void MargeIdepthsShort();
    void MargInsideBias();
    void MargInsidePos();

    void UpdateMargposUseGradientAndGauss(double gradient_scale, double gauss_scale);
    void UpdateMargPhaseBiasUseGradientAndGauss(double gradient_scale, double gauss_scale);
    void UpdateMargBiasUseGradientAndGauss(double gradient_scale, double gauss_scale);
    void UpdateMargFeatureUseGradientAndGauss(double gradient_scale, double gauss_scale);
    void UpdateInsideStateUseGradientAndGauss(double gradient_scale, double gauss_scale);

    void UpdateInsideGaussStep();
    void UpdateMargposGaussStep();
    void UpdateMargBiasGaussStep();
    void UpdateMargPhaseBiasGaussStep();
    void UpdateMargFeatureGaussStep();

    void SaveOrRestoreHidenStates(bool );

    void EvaluateCost(double* cost, double* model_cost_change);
    void EvaluateAlpha(double* alpha1);

    void SaveGradientOutside();
    void SaveGaussStepOutside();
    void SaveGaussStep();

    int GetBIndex(int);
    int GetPIndex(int);

    bool EvaluateLhsRhs(double const* const* parameters, double* rhs, double* lhs, double* gradient, double mu);

    void SaveOutsidePointerRaw();
    void setRPwc();

    void VisualEvaluate(int pindexi, int pindexj,
                        double* parameters, double* residuals, double** jacobians,
                        const Eigen::Vector3d pts_i, const Eigen::Vector3d pts_j
#if USE_STEREO
                        , bool is_stereo
#endif
                       );
    void EvaluateLastMargInfo();
    void AddLastMargeInfo(MarginalizationInfo* last_marg_info_);
    void ConstructPriorIdx();
    void EvaluatePriorAlpha(double* alpha1);
    void EvaluatePriorModelCostChange(double* model_cost_change);
    void EvaluatePriorCost(double* cost);

    void UpdateGaussAndGradientOutside();
    void UpdateLhsRhsGradientNew();
    void ForwardDeliverIdepthValues();
    void UpdateStateNew(double gauss_scale, double gradient_scale);

    void RemoveFeature(int feature_id);
    void SaveLoadCandidateResidual(bool is_save);
    void ComputeRhs();
    void VisualEvaluateRhs(int pindexi, int pindexj,
                           double* rhs_idepth,
                           double* jacobian_pointer,
                           double* visual_residual_raw_in
#if USE_STEREO
                           , double* lhs_scale_idepth, bool is_stereo
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
                           , double* lhs_extrinsic_idepth
#endif
#endif
                          );
    void MargeRhs();
    void EvaluateNonlinearQuantity();
    void MargRhsNew();

    void BuildForwardMatrixInfo();
    void BuildCurPrevIdx();
    void GetMatrixNext(flower_focus::MatrixXd& lhs, Eigen::VectorXd& rhs, std::vector<double*>& new_parameters, std::vector<int>& new_sizes);
    void Reset() {
        if (memory) {
            delete memory; memory = 0;
        }
    }
    void ResetInit() {
        init = false;
        if (next_vibase)next_vibase->ResetInit();
    }

    int visual_residual_short_num;
    int visual_residual_long_num;
    int visual_lhs_count;

    flower_focus::VectorXd prior_residual;
    flower_focus::VectorXd prior_residual_save;

    std::map<int, VisualFactor*, less<int>>idepth_map_factors_short;
    std::map<int, VisualFactor*, less<int>>idepth_map_factors_long;
    std::map<int, IMUGNSSPreFactor*, less<int>>imu_map_factors;
    std::vector<VisualFactor*>long_feature_factors;


    LossFunction* cauchy_loss_function = 0;

    double* memory = 0;

    double* lhs_pos_idepth_global;
    double* old_estimations_pointer;

    double* visual_jacobians;
    double* visual_residuals;
    double* visual_costs;
    double* visual_residuals_save;

    double* visual_model_residual_accum;
    double* visual_model_residual_accum_save;
    double* visual_residuals_linerized;

    std::vector<double*> lhs_posbias;
    double* rhs_posbias_pointer;

    double* rhs_posbias_save_pointer;

    double* lhs_margpos_outside_pointer;
    double* lhs_margpos_pointer;
    double* rhs_margpos_pointer;

    double* inc_posbias_pointer;
    double* inc_idepth_short_pointer;
    double* inc_outside_pointer;
    double* gauss_posbias_pointer;
    double* gauss_idepth_short_pointer;

    double* idepth_short_save_pointer;
    double* para_posebias_inside_save_pointer;

    std::vector<double*>outside_pointer_raw;

    std::unordered_map<int, double*>pos_index2pointers;
    std::unordered_map<double*, int>pointer2orders;
    std::unordered_map<int, double*>order2pointers;
    std::unordered_map<double*, int>pos_pointer2indexs;
    std::unordered_map<int, double*>bias_index2pointers;
    std::vector<Eigen::Matrix3d>Rwcs;
    std::vector<Eigen::Vector3d>Pwcs;


    std::vector<int> pindex2order;
    std::vector<int> bindex2order;
    std::vector<int> order2size;
    std::vector<int> order2p_local;
    std::vector<int> order2p_global;
    std::vector<std::vector<int>> connections_bias;

    double gradient_squared_norm_inside;
    double gradient_dot_gauss_newton_inside;
    double gauss_newton_step_inside_square_norm;


    int marg_count;

    int bias_para_count;
    int pos_para_count;
    int para_count;

    int idepth_long_count;

    double mu = 0;


    int O_MARG_BIAS;
    int E_O_MARG_BIAS;
    int O_PHASE_BIAS_SHORT;
    int O_MARG_POS;
    int E_O_MARG_POS;
    int O_POSR;
    int O_BIAS0;
    int O_POS0;
    int O_PHASE_BIAS_LONG;
    int O_POSK;
    int O_BIASK;
    int O_POSR_NEXT;
    int O_SCALE;
    int O_POSGLOBAL;
    int O_TD;
    int O_EXTRINSIC;
    int O_ACC_S;
    int O_EXTRINSIC_RIGHT;
    int O_PBG;
    int O_IDEPTHR;
    int O_FULL;

    MarginalizationInfo* last_marg_info = 0;
    flower_focus::VectorXd prior_rhs;
    flower_focus::VectorXd prior_rhs_save;

    Eigen::Matrix<double, 1, 6, Eigen::RowMajor> yaw_constraint_jacobian;
    flower_focus::Vector1d yaw_constraint_residual;
    flower_focus::Vector1d yaw_constraint_residual_save;

    flower_focus::Vector1d yaw_constraint_model_residual_accum;
    flower_focus::Vector1d yaw_constraint_residual_linerized;
    flower_focus::Vector1d yaw_constraint_model_residual_accum_save;

    Eigen::Vector3d InitMag;

    std::vector<int>prior_idx;

    VisualInertialBase* prev_vibase = 0;
    VisualInertialBase* next_vibase = 0;
    std::unordered_map<double*, double*>deliver_idepth_pointer;//1
    std::unordered_map<double*, int>pointer2p_local;

    Eigen::VectorXd gauss_outside;
    flower_focus::MatrixXd lhs_outside;
    Eigen::VectorXd rhs_outside;
    Eigen::VectorXd gradient_outside;

    int num_alone = 0;
    std::vector<MatirxInfo>matrix_info;
    std::vector<int>cur_idxs, prev_idxs, cur_sizes;
    int outside_threshold;
    bool ready_for_sequential_forward = false;
    bool ready_for_sequential_backward = false;

    double* global_pos_pointer = 0;
    double* rwgw_pointer = 0;
    double* scale_pointer = 0;
    double distance0;
    Eigen::VectorXd init_distance_constraint_jacobian;
    double init_distance_constraint_residual;
    double init_distance_constraint_residual_save;

    double init_distance_constraint_model_residual_accum;
    double init_distance_constraint_residual_linerized;
    double init_distance_constraint_model_residual_accum_save;
    double last_cost = 0;
    double* td_pointer = 0;
    double* extrinsic_pointer = 0;
    bool init = false;


    double nonlinear_quantity;
    bool history_flag;
    double nonlinear_quantity_deliver_accum = 1e8;
    double* acc_scale_pointer = 0;
    double* pbg_pointer = 0;


    std::set<double*>phase_bias_pointers_set;
    std::vector<double*>long_phase_bias;
    std::vector<double*>short_phase_bias;
    std::map<double*, int>long_phase_bias_idx_map;
    std::map<double*, int>short_phase_bias_idx_map;
    int phase_bias_alone_num;

    void MargeShortPhaseBias();
    void BuildPhaseBias();


#if USE_STEREO
    Eigen::Vector3d tC2C1;
    Eigen::Matrix3d rC2C1;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    double* extrinsic2_pointer;
#endif
#endif
};


extern VisualInertialBase* VI_info_pointer;
extern VisualInertialBase* VI_info_pointer_head;
extern VisualInertialBase* VI_info_pointer_tail;
