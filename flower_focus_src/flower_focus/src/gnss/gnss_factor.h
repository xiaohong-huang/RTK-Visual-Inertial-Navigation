#pragma once
#include "../solver/cost_function.h"
#include <eigen3/Eigen/Dense>
#include "../parameter/parameters.h"


// RB-SD carr-er-phase factor
class RTKCarrierPhaseFactor2 : public SizedCostFunction<1, 7, 1, 1> {
  public:
    RTKCarrierPhaseFactor2(double* satelite1_,
                          double L1_lam_,
                          double lam_,
                          double el_,
                          double base_rover_time_diff_,
                          double mea_var_,
                          double* base_pos_,
                          bool use_istd_, int sys_, int f_):
        satelite_pos(satelite1_),
        L1_lam(L1_lam_),
        lam(lam_),
        el(el_),
        base_rover_time_diff(base_rover_time_diff_),
        mea_var(mea_var_),
        base_pos(base_pos_),
        use_istd(use_istd_), sys(sys_), f(f_) {
        ASSERT(el != 0);
    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    double* satelite_pos; //satellite position
    double L1_lam;//RB-SD carrier-phase multiply wave length [m]
    double lam;//wave length [m]
    double el;//rover to satellite elevation
    double base_rover_time_diff;//time difference between rover frame and base frame
    double mea_var;//measurement covariance
    double* base_pos; //base position
    bool use_istd;//sqrt info is 1 when use_istd is 0
    int sys; //satellite system;0 for GPS;1 for Beidou;2 For Gaileo
    int f;//frequency of measurement
};

// RB-SD pseudorange factor
class RTKPseudorangeFactor2: public SizedCostFunction<1, 7, 1> {
  public:
    RTKPseudorangeFactor2(double* satelite1_,
                         double P1_,
                         double el_,
                         double base_rover_time_diff_,
                         double mea_var_,
                         double* base_pos_):
        satelite_pos(satelite1_),
        P1(P1_),
        el(el_),
        base_rover_time_diff(base_rover_time_diff_),
        mea_var(mea_var_),
        base_pos(base_pos_) {
        ASSERT(el != 0);
    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    double* satelite_pos; //satellite position
    double P1;//RB-SD pseudorange measurement
    double el;//rover to satellite elevation
    double base_rover_time_diff;//time difference between rover frame and base frame
    double mea_var;//measurement covariance
    double* base_pos; //base position

};

// rover-only pseudorange factor
class SppPseudorangeFactor2 : public SizedCostFunction<1, 7, 1> {
  public:
    SppPseudorangeFactor2(double* satelite1_, double P1_, double istd_, double* base_pos_):
        satelite_pos(satelite1_),
        P1(P1_),
        istd(istd_),
        base_pos(base_pos_) {}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    double* satelite_pos; //satellite position
    double P1;//rover-only pseudorange
    double istd;//fusion weight
    double* base_pos; //base position

};



// rover-only carrier-phase factor
class SppCarrierPhaseFactor2 : public SizedCostFunction<1, 7, 1, 1> {
  public:
    SppCarrierPhaseFactor2(double* satelite1_, double L1_lam_, double istd_, double* base_pos_, double lam_):
        satelite_pos(satelite1_),
        L1_lam(L1_lam_),
        istd(istd_),
        base_pos(base_pos_),
        lam(lam_) {}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    double* satelite_pos; //satellite position
    double L1_lam;//rover-only carrier-phase multiply wave length [m]
    double istd;//weight
    double* base_pos; //base position
    double lam;//wave length

};

// Doppler factor
class SppDopplerFactor2 : public SizedCostFunction<1, 9, 1, 7> {
  public:
    SppDopplerFactor2(double* satelitev1_, double* satelite1_, double* xyzt_, double D1_lam_, double istd_, double* base_pos_):
        satelitev1(satelitev1_),
        satelite_pos(satelite1_),
        D1_lam(D1_lam_),
        istd(istd_),
        base_pos(base_pos_) {

    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;


    double* satelitev1; //satellite velocity
    double* satelite_pos; //satellite position
    double D1_lam;//Doppler multiply wave length [m/s]
    double istd;  //weight
    double xyzglobal[3];//rover position under ECEF frame
    double* base_pos;//base position


};



// fixed integer factor
class FixedIntegerFactor2 : public SizedCostFunction<1, 1, 1> {
  public:
    FixedIntegerFactor2(double N21_, double istd_): N21(N21_), istd(istd_) {}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;
    double N21;//fixed integer
    double istd;//weight


};

class InitialBlackFactor2 : public SizedCostFunction<1, 1> {
  public:
    InitialBlackFactor2(double istd_): istd(istd_) {
    }
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const{
    residuals[0] = parameters[0][0] * istd;
    if (jacobians) {
        if (jacobians[0])jacobians[0][0] = 1 * istd;
    }
    return true;
}
    double istd;
};
