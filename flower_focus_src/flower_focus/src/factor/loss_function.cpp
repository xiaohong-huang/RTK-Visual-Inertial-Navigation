#include "loss_function.h"
//this file is from Ceres-Solver https://github.com/ceres-solver/ceres-solver
void CauchyLoss::Evaluate(double s, double rho[3]) const {
    const double sum = 1.0 + s * c_;
    const double inv = 1.0 / sum;
    rho[0] = b_ * log(sum);
    rho[1] = std::max(std::numeric_limits<double>::min(), inv);
    rho[2] = -c_ * (inv * inv);
}

void ArctanLoss::Evaluate(double s, double rho[3]) const {
    const double sum = 1 + s * s * b_;
    const double inv = 1 / sum;
    rho[0] = a_ * atan2(s, a_);
    rho[1] = std::max(std::numeric_limits<double>::min(), inv);
    rho[2] = -2.0 * s * b_ * (inv * inv);
}

void HuberLoss::Evaluate(double s, double rho[3]) const {
  if (s > b_) {
    const double r = sqrt(s);
    rho[0] = 2.0 * a_ * r - b_;
    rho[1] = std::max(std::numeric_limits<double>::min(), a_ / r);
    rho[2] = -rho[1] / (2.0 * s);
  } else {
    rho[0] = s;
    rho[1] = 1.0;
    rho[2] = 0.0;
  }
}

void SoftLOneLoss::Evaluate(double s, double rho[3]) const {
  const double sum = 1.0 + s * c_;
  const double tmp = sqrt(sum);
  
  rho[0] = 2.0 * b_ * (tmp - 1.0);
  rho[1] = std::max(std::numeric_limits<double>::min(), 1.0 / tmp);
  rho[2] = -(c_ * rho[1]) / (2.0 * sum);
}


double SoftLOneLoss::Evaluate0(double s) const {
  const double sum = 1.0 + s * c_;
  const double tmp = sqrt(sum);
  return 2.0 * b_ * (tmp - 1.0);
}

double CauchyLoss::Evaluate0(double s) const {
    const double sum = 1.0 + s * c_;
    return b_ * log(sum);
}

double ArctanLoss::Evaluate0(double s) const {
    return a_ * atan2(s, a_);
}

double HuberLoss::Evaluate0(double s) const {
  if (s > b_) {
    const double r = sqrt(s);
    return 2.0 * a_ * r - b_;
  } else {
    return s;
  }
}
