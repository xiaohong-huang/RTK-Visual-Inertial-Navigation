#pragma once
#include<vector>
//this file is from Ceres-Solver https://github.com/ceres-solver/ceres-solver

class CostFunction {
 public:
  CostFunction() : num_residuals_(0) {}
  CostFunction(const CostFunction&) = delete;
  void operator=(const CostFunction&) = delete;

  virtual ~CostFunction() {}

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const = 0;

  virtual bool EvaluateLhsRhs(double const* const* parameters,
                        double* rhs,
                        double* lhs,double* gradient,double radius) {
                          return false;
                        }
  virtual bool SaveOrRestoreHidenStates(bool is_save){
                        return false;
  }

  virtual bool EvaluateCost(double const* const* parameters,double*cost){
                        return false;
  }


  virtual bool  UpdateInsideGaussStep(double const* const* parameters,double*gradient,double*gauss,
        double*gauss_newton_step_inside_square_norm,double*alpha1,double*gradient_squared_norm,double*gradient_dot_gauss_newton_inside){
                        return false;
        }
virtual bool  UpdateInsideStateAndEvaluateCost(double const* const* parameters,double*cost,double*model_cost_change,
        double gradient_scale,double gauss_scale){
                        return false;
        }

  const std::vector<int>& parameter_block_sizes() const {
    return parameter_block_sizes_;
  }

  int num_residuals() const { return num_residuals_; }

 protected:
  std::vector<int>* mutable_parameter_block_sizes() {
    return &parameter_block_sizes_;
  }

  void set_num_residuals(int num_residuals) { num_residuals_ = num_residuals; }

 private:
  // Cost function signature metadata: number of inputs & their sizes,
  // number of outputs (residuals).
  std::vector<int> parameter_block_sizes_;
  int num_residuals_;
};

enum DimensionType {
  DYNAMIC = -1,
};




template <int kNumResiduals, int... Ns>
class SizedCostFunction : public CostFunction {
 public:


  SizedCostFunction() {
    set_num_residuals(kNumResiduals);
    *mutable_parameter_block_sizes() = std::vector<int>{Ns...};
  }

  virtual ~SizedCostFunction() {}

};