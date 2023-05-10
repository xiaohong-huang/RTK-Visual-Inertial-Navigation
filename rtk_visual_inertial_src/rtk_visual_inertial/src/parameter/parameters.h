

#pragma once


#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>


using namespace std;
//maximum duration for a avaliable imu pre-integration measurement
#define IMUINTEGRATIONTIMEMAX 10.0
//reduce the visual weight
#define FEATUREWEIGHTINVERSE 1.5
//focal length
#define FOCAL_LENGTH  1000.0
//maximum number of visual frames
#define FEATURE_WINDOW_SIZE 10
//maximum number of gnss frames
#define GNSS_WINDOW_SIZE 500
//feature point is avaliable when tracking frames are larger than FEATURE_CONTINUE
#define FEATURE_CONTINUE 2
//using the inverse depth as residual rather than the fearture position
#define USE_INVERSE_DEPTH 0
//number of epoch for estimaing the pseudorange correction [0.1s]
#define EstimatePcorrectionPerio 500
//inial depth values for feature points
#define INIT_DEPTH 5.0












//use visual measurements,must be 1 in current version
extern bool USE_IMAGE;
//use gnss measurements
extern bool USE_GNSS;
//use imu measuremtns
extern bool USE_IMU;
//use RB-SD pseudorange and carrier-phase measureemtns
extern bool USE_RTK;
//use RB-SD pseudorange measurements
extern bool USE_RTD;
//use doppler measurements
extern bool USE_DOPPLER;
//use rover-only carrier-phase measurements
extern bool USE_SPP_PHASE;
//use magnetormetor for yaw initialization
extern bool USE_MAG_INIT_YAW;
//correcting yaw using a magnetormetor when optimizing the states of the SWF
extern bool USE_MAG_CORRECT_YAW;
//obtain the covariance of phase biases when solving the Gauss-Newton problem
extern bool USE_DIRECT_N_RESOLVE;
//resolve the integer ambiguity
extern bool USE_N_RESOLVE;
//estimate the spp correction when first observing the pseudorange measurements. usefull for the satellites unstable-tracking environment
extern bool USE_SPP_CORRECTION;
//use a global ceres solver for solving states
extern bool USE_GLOBAL_OPTIMIZATION;
//use stereo visual measurements.
extern bool USE_STEREO;
//path for saving the extrinsic estimation.
extern std::string EX_CALIB_RESULT_PATH;
//path for saving the estimation
extern std::string RESULT_PATH;
//imu ros topic
extern std::string IMU_TOPIC;
//bag path
extern std::string ROS_PATH;
//other topic
extern std::string IMAGE0_TOPIC, IMAGE1_TOPIC, RTK_TOPIC, MAG_TOPIC, PPS_LOCAL_TOPIC, PPS_GPS_TOPIC, FEATURE_TOPIC;
//path of intrinc calibration
extern std::vector<std::string> CAM_NAMES;

//magnetometer installation rotation matrix
extern Eigen::Matrix3d MagMatrix;
//imu installation rotation matrix
extern Eigen::Matrix3d IMUMatrix;
//extrinsic calibration from imu to antenna
extern Eigen::Vector3d Pbg;
//anchor point for calculating the position under ENU frame.
extern Eigen::Vector3d ANCHOR_POINT;
//magnetometer installation translation matrix
extern Eigen::Vector3d MagVector;
//rotation from ECEF frame to ENU frame of the body frame
extern Eigen::Matrix3d Rwgw;
//rotation from camera to imu
extern std::vector<Eigen::Matrix3d> RIC;
//translation from camera to imu
extern std::vector<Eigen::Vector3d> TIC;
//gravity vector
extern Eigen::Vector3d G;
//use the first AVERAGE_IMU numbert of IMU measurements to initialize the states of the system.
extern double AVERAGE_IMU;
//skip the first SKIP_IMU number of IMU measurements.
extern double SKIP_IMU;
//skip the first AVERAGE_IMAGE number of visual measurements.
extern double AVERAGE_IMAGE;
//Levenberg trust region radius
extern double MAX_TRUST_REGION_RADIUS;
//ansac threshold (pixel)
extern double F_THRESHOLD;
//use pre-compute feature points. Set this to true can make the results more repeatable.
extern double USE_FEATURE;
// imu noise standard deviation
extern double ACC_N, ACC_W, GYR_N, GYR_W;
//max solver time for solving the Gauss-Newton problem
extern double MAX_SOLVER_TIME;
//keyframe selection threshold (pixel)
extern double MIN_PARALLAX;
//use the phase bias if the observing epochs is larger than CARRIER_PHASE_CONTINUE_THRESHOLD.
extern int CARRIER_PHASE_CONTINUE_THRESHOLD;
//fix the integer if the ratio-test sucess for FIX_CONTINUE_THRESHOLD epochs.
extern int FIX_CONTINUE_THRESHOLD;
//reset the phase biases if the fixed solutions are unavaliable for Phase_ALL_RESET_COUNT epochs.
extern int Phase_ALL_RESET_COUNT;
//publish tracking image as topic
extern int SHOW_TRACK;
//perform forward and backward optical flow to improve feature tracking accuracy
extern int FLOW_BACK;
//max feature number in feature tracking
extern int MAX_CNT;
//min distance between two features
extern int MIN_DIST;
//number of camera
extern int NUM_OF_CAM;
//estimate camera-imu calibration
extern int ESTIMATE_EXTRINSIC;
//max number iteration when solving the Gauss-Newton problem
extern int MAX_NUM_ITERATIONS;


void readParameters(std::string config_file);

enum SIZE_PARAMETERIZATION {
    SIZE_POSE = 7,
    SIZE_SPEEDBIAS = 9,
    SIZE_FEATURE = 1
};

enum StateOrder {
    O_P = 0,
    O_R = 3,
    O_V = 6,
    O_BA = 9,
    O_BG = 12
};

enum NoiseOrder {
    O_AN = 0,
    O_GN = 3,
    O_AW = 6,
    O_GW = 9
};
