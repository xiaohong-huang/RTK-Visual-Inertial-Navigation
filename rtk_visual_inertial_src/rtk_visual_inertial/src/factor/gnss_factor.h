#pragma once
#include <ceres/ceres.h>
#include <eigen3/Eigen/Dense>



// RB-SD carr-er-phase factor
class RTKCarrierPhaseFactor : public ceres::SizedCostFunction<1, 7, 1, 1> {
  public:
    RTKCarrierPhaseFactor(double* satelite1_,
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
        assert(el != 0);
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
class RTKPseudorangeFactor: public ceres::SizedCostFunction<1, 7, 1> {
  public:
    RTKPseudorangeFactor(double* satelite1_,
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
        assert(el != 0);
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
class SppPseudorangeFactor : public ceres::SizedCostFunction<1, 7, 1> {
  public:
    SppPseudorangeFactor(double* satelite1_, double P1_, double istd_, double* base_pos_):
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
class SppCarrierPhaseFactor : public ceres::SizedCostFunction<1, 7, 1, 1> {
  public:
    SppCarrierPhaseFactor(double* satelite1_, double L1_lam_, double istd_, double* base_pos_, double lam_):
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
class SppDopplerFactor : public ceres::SizedCostFunction<1, 9, 1, 7> {
  public:
    SppDopplerFactor(double* satelitev1_, double* satelite1_, double* xyzt_, double D1_lam_, double istd_, double* base_pos_):
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
class FixedIntegerFactor : public ceres::SizedCostFunction<1, 1, 1> {
  public:
    FixedIntegerFactor(double N21_, double istd_): N21(N21_), istd(istd_) {}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;
    double N21;//fixed integer
    double istd;//weight


};

