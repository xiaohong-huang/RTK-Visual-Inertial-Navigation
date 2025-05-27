
#pragma once
#include <eigen3/Eigen/Dense>
//this file is from Ceres-Solver https://github.com/ceres-solver/ceres-solver

class  LossFunction {
  public:
    virtual ~LossFunction() {}
    virtual void Evaluate(double sq_norm, double out[3]) const = 0;
    virtual double Evaluate0(double sq_norm) const = 0;
};


class CauchyLoss : public LossFunction {
  public:
    explicit CauchyLoss(double a) : b_(a * a), c_(1 / b_) {}
    void Evaluate(double, double*) const override;
    double Evaluate0(double) const override;

  private:
    const double b_;
    const double c_;
};

class ArctanLoss : public LossFunction {
  public:
    explicit ArctanLoss(double a) : a_(a), b_(1 / (a * a)) {}
    void Evaluate(double, double*) const override;
    double Evaluate0(double) const override;

  private:
    const double a_;
    const double b_;
};


class HuberLoss : public LossFunction {
 public:
  explicit HuberLoss(double a) : a_(a), b_(a * a) {}
  void Evaluate(double, double*) const override;
  double Evaluate0(double) const override;

 private:
  const double a_;
  const double b_;
};


class SoftLOneLoss : public LossFunction {
 public:
  explicit SoftLOneLoss(double a) : b_(a * a), c_(1 / b_) {}
  void Evaluate(double, double*) const override;
  double Evaluate0(double) const override;

 private:
  const double b_;
  const double c_;
};