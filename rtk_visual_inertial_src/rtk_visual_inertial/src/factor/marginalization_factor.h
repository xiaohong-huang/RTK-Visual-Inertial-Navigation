

#pragma once

#include <cstdlib>
#include <pthread.h>
#include <ceres/ceres.h>
#include <unordered_map>

//very poor
const int NUM_THREADS = 8;

struct ResidualBlockInfo {
    ResidualBlockInfo(ceres::CostFunction* _cost_function, ceres::LossFunction* _loss_function, std::vector<double*> _parameter_blocks, std::vector<int> _drop_set, std::vector<int> _sparse_set)
        : cost_function(_cost_function), loss_function(_loss_function), parameter_blocks(_parameter_blocks), drop_set(_drop_set), sparse_set(_sparse_set) {}

    void Evaluate();

    ceres::CostFunction* cost_function;
    ceres::LossFunction* loss_function;
    std::vector<double*> parameter_blocks;
    std::vector<int> drop_set;
    std::vector<int> sparse_set;

    double** raw_jacobians;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobians;
    Eigen::VectorXd residuals;

    int localSize(int size) {
        return size == 7 ? 6 : size;
    }
};

struct ThreadsStruct {
    std::vector<ResidualBlockInfo*> sub_factors;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    std::unordered_map<long, int> parameter_block_size; //global size
    std::unordered_map<long, int> parameter_block_idx; //local size
};

class MarginalizationInfo {
  public:
    void checksparse(Eigen::MatrixXd& Am2m2);
    MarginalizationInfo() {};
    MarginalizationInfo(MarginalizationInfo* old) {
        keep_block_size = old->keep_block_size;
        keep_block_idx = old->keep_block_idx;
        keep_block_data = old->keep_block_data;
        keep_block_addr = old->keep_block_addr;
        linearized_jacobians = old->linearized_jacobians;
        linearized_residuals = old->linearized_residuals;
        n = old->n;
        m = old->m;
    }
    void setmarginalizeinfo(std::vector<const double*>& parameter_block_addr,
                            std::vector<int>& parameter_block_global_size, Eigen::MatrixXd& A_, Eigen::VectorXd& b_,
                            bool Sqrt);
    ~MarginalizationInfo();
    int localSize(int size) const;
    int globalSize(int size) const;
    void addResidualBlockInfo(ResidualBlockInfo* residual_block_info);
    //marginalizing the selected states
    void marginalize(bool initialinformation, bool multithread);
    //marginalizing except the selected states
    void marginalize_exclude();
    std::vector<double*> getParameterBlocks();
    //reset the linearization point of prior
    void ResetLinearizationPoint(std::vector<double*>parameters) ;

    std::vector<ResidualBlockInfo*> factors;
    int m, n;
    std::unordered_map<long, int> parameter_block_size; //global size
    int sum_block_size;
    std::unordered_map<long, int> parameter_block_idx; //local size
    std::unordered_map<long, int> parameter_block_sparse_idx; //local size
    std::unordered_map<long, int> parameter_block_drop_idx; //local size

    std::unordered_map<long, double*> parameter_block_data;

    std::vector<int> keep_block_size; //global size
    std::vector<int> keep_block_idx;  //local size
    //linerization points for the remained states
    std::vector<double*> keep_block_data;
    //prior jacobian
    Eigen::MatrixXd linearized_jacobians;
    //prior residual
    Eigen::VectorXd linearized_residuals;

    const double eps = 1e-8;
    //hessian matrix
    Eigen::MatrixXd A;
    //right hand side matrix
    Eigen::VectorXd b;
    //pointer to the remained states
    std::vector<double*> keep_block_addr;
    //for trust region updation
    const double min_diagonal = 1e-6;
    const double max_diagonal = 1e32;


};

class MarginalizationFactor : public ceres::CostFunction {
  public:
    MarginalizationFactor(MarginalizationInfo* _marginalization_info);
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    MarginalizationInfo* marginalization_info;
};
