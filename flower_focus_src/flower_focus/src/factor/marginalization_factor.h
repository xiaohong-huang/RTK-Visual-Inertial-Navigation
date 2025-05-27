

#pragma once

#include <cstdlib>
#include <pthread.h>
#include <unordered_map>
#include "../parameter/parameters.h"
//very poor
#include "../solver/cost_function.h"
#include "loss_function.h"

#include "sparse_matrix.h"
struct ResidualBlockInfo {
    ResidualBlockInfo(CostFunction* _cost_function, LossFunction* _loss_function, std::vector<double*> _parameter_blocks, std::vector<int> _drop_set, std::vector<int> _sparse_set)
        : cost_function(_cost_function), loss_function(_loss_function), parameter_blocks(_parameter_blocks), drop_set(_drop_set), sparse_set(_sparse_set) {}

    void Evaluate();

    CostFunction* cost_function;
    LossFunction* loss_function;
    std::vector<double*> parameter_blocks;
    std::vector<int> drop_set;
    std::vector<int> sparse_set;

    double** raw_jacobians;
    std::vector<flower_focus::MatrixXd > jacobians;
    Eigen::VectorXd residuals;


};

struct ThreadsStruct {
    std::vector<ResidualBlockInfo*> sub_factors;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    std::unordered_map<double*, int> parameter_block_size; //global size
    std::unordered_map<double*, int> parameter_block_idx; //local size
};

class MarginalizationInfo {
  public:
    MarginalizationInfo() {};
    MarginalizationInfo(MarginalizationInfo* old) {
        keep_block_size = old->keep_block_size;
        keep_block_idx = old->keep_block_idx;
        keep_block_data = old->keep_block_data;
        keep_block_addr = old->keep_block_addr;
        A = old->A;
        b = old->b;
        n = old->n;
        m = old->m;
    }
    void setmarginalizeinfo(std::vector<double*>& parameter_block_addr,
                            std::vector<int>& parameter_block_global_size, Eigen::MatrixXd A_, Eigen::VectorXd& b_,
                            bool Sqrt);
    ~MarginalizationInfo();
    void addResidualBlockInfo(ResidualBlockInfo* residual_block_info);
    void marginalize(bool initialinformation);
    void getParameterBlocks(bool);
    void resetLinerizationPoint();
    void GetJacobianAndResidual();
    void marginalize_pointers(std::set<double*>feature_marge_pointer,bool);

    int m, n;
    double cost0;

    std::vector<ResidualBlockInfo*> factors;
    std::unordered_map<double*, int> parameter_block_size; //global size
    std::unordered_map<double*, int> parameter_block_idx; //local size
    std::unordered_map<double*, int> parameter_block_drop_idx; //local size
    std::vector<double*>parameters;
    std::vector<double*>drop_parameters;
    std::unordered_map<double*, double*> parameter_block_data;
    std::vector<int> keep_block_size; //global size
    std::vector<int> keep_block_idx;  //local size
    std::vector<double*> keep_block_data;
    std::vector<double*> keep_block_addr;
    std::unordered_map<double*, int>keep_block_addr_set;
    


    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::MatrixXd linearized_jacobians;
    Eigen::VectorXd linearized_residuals;
    double init_cost;



    void BuildMatrixF(std::vector<double*>pointers);
    void BuildMatrixG(std::vector<double*>pointers, Eigen::VectorXd dx, Eigen::VectorXd& dy);

    Eigen::MatrixXd FTAF;
    Eigen::VectorXd FTb;
    Eigen::MatrixXd linearized_jacobiansF;
    
    int v_i;
    std::vector<int>p_is;
    std::vector<int>lambda_is;
    int extrinsic_i;
    Eigen::Vector3d gyr0;
    void compute_residual_rhs(Eigen::VectorXd& residual, Eigen::VectorXd& rhs);
    bool fix = false;
    std::vector<MatirxInfo>matrix_infoF;
    std::vector<MatirxInfo>matrix_infoG;
    void UpdateMatix(const Eigen::MatrixXd& A, Eigen::MatrixXd& A2);
    void UpdateVector(const Eigen::VectorXd& b, Eigen::VectorXd& b2, bool onlyF);
    void UpdateVectorTranspose(const Eigen::VectorXd& b, Eigen::VectorXd& b2);
    int R_WI_WC_i;
    int scale_factor_i;
    void Getestimate(Eigen::Vector3d& v_c0_c0, Eigen::Vector3d& P_I_ck0, Eigen::Matrix3d& R_I_ck, double& sum_distance, std::vector<double>& idepth0, Eigen::Matrix3d& R_WI_WC, std::vector<double*>pointer);
    void UpdateF();
    std::vector<Eigen::MatrixXd> idepth_pos;
    int max_idepth_idx, min_idepth_idx;
    
    
};


class MarginalizationFactor : public CostFunction {
  public:
    MarginalizationFactor(MarginalizationInfo* _marginalization_info);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    MarginalizationInfo* marginalization_info;
};
