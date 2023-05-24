#pragma once
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <vector>
#include <list>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include<vector>
#include "marginalization_factor.h"
#include "ceres/problem.h"
//common_function.c


#define LTHRE 15

#define clight  299792458.0
#define AZELMIN (25.0/180*PI)

#define NFREQ 2

#define MAXOBS  64


#define SYS_GPSS 0
#define SYS_CMP 1
#define SYS_GALL 2


#define MAXSATNUMGPS 39
#define MAXSATNUMCMP (39 + 37)
#define MAXSATNUMGAL (39 + 37 + 31)
#define MAXSATNUM  (39 + 37 + 31)

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */

#define SQR(x)   ((x)*(x))
#define ABS(x) ((x)>0?(x):(-(x)))
#define PI          3.1415926535897932  /* pi */






typedef struct {
    int time;
    double sec;
} gtime_tt;


class PBtype {
  public:
    double value;//phase bias estimation
    std::map<PBtype*, int>bestNtype; //not use
    bool use;//being observed by the remaind states
    bool isbest;//not use
    uint8_t LLI;//cycle slip flag
    uint8_t half_flag;//half flag
    uint8_t sys;//satellite system
    uint8_t f;//measured freequency
    int continue_count;//epoch for continuing tracked
    double last_update_time;//time stamp for the last observation
};


typedef struct {
    uint8_t sat;//satellite number
    uint8_t sys;//satellite system
    uint8_t SVH;//satellite health

    uint8_t	RTK_LLI[NFREQ];//RB-SD carrier-phase cycle slip flag
    uint8_t	SPP_LLI[NFREQ];//Rover-only carrier-phase cycle slip flag
    uint8_t half_flag[NFREQ];//RB-SD carrier-phase half flag

    double SPP_P[NFREQ];//Rover-only pseudorange measurement(with correction) [m]
    double SPP_L[NFREQ];//Rover-only carrier-phase measurement(with correction) [cycle]
    double SPP_D[NFREQ];//Rover-only doppler measurement(with correction) [cycle/s]
    double SPP_Lstd[NFREQ];//Rover-only pseudorange measurement std [m]
    double SPP_Pstd[NFREQ];//Rover-only carrier-phase measurement std [cycle]
    double SPP_Dstd[NFREQ];//Rover-only doppler measurement std [cycle/s]


    double RTK_P[NFREQ];//RB-SD pseudorange measurement[m]
    double RTK_L[NFREQ];//RB-SD carrier-phase measurement[cycle]
    double RTK_Pstd[NFREQ];//RB-SD pseudorange measurement std[m]
    double RTK_Lstd[NFREQ];//RB-SD carrier-phase measurement std[cycle]


    double SPP_P0[NFREQ];	//m

    double satellite_pos[3];// satellite position [m]
    double satellite_vel[3];//satellite velocity [m/s]

    double el;// elevation [rad]


    double sat_var; //satellite noise covariance [m^2]
    double ion_var; //iono noise covariance [m^2]
    double trop_var; //trop noise covariance [m^2]

    PBtype* RTK_Npoint[NFREQ];//RB-SD  Ambiguity pointer
    PBtype* SPP_Npoint[NFREQ];//Rover-only  Ambiguity pointer
    PBtype* SPP_Npoint_PCottections[NFREQ];//Pseudorange correction  Ambiguity pointer


} ObsMea;



typedef struct mea_t {
    int obs_count;//satellites observation number
    double ros_time;//fpga recording time stamp
    ObsMea obs_data[MAXOBS];//satellites observation measurements
    double base_xyz[3];//base position
    MarginalizationInfo* marg_info_gnss;// GNSS preprocess  information
    ceres::ResidualBlockId residualBlockId;// residualblockID in Ceres Solver
    double br_time_diff;//time difference between rover frame and base frame.
} mea_t;




void ecef2pos(const double* r, double* pos);
double satazel(const double* pos, const double* e, double* azel);
double distance(const double* rr, const double* rs, double* e);
void xyz2enu(const double* pos, double* E);
void ecef2enu(const double* pos, const double* r, double* e);
void matmul(const char* tr, int n, int k, int m, double alpha, const double* A, const double* B, double beta, double* C);
void matcpy(double* A, const double* B, int n, int m);
double* mat(int n, int m);
int* imat(int n, int m);
uint8_t matinv(double* A, int n);
uint8_t matinv33(double* a, double* b);
uint8_t matinv33f(float* a, float* b);
uint8_t lsq6( double* A,  double* y, int n, int m, double* x, double* Q);
uint8_t lsq3(double* A, double* y, int n, int m, double* x, double* Q);
int lsq(double* A, double* y, int n, int m, double* x, double* Q);
double dot(const double* a, const double* b, int n);
double norm(const double* a, int n);
double velecitydistance(const double* rr, const double* rs, const double* vr, const double* vs, double* e);
void update_azel(double globalxyz[3], mea_t* rover);

extern double lams[3][2];


