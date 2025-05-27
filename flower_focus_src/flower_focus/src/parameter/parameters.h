

#pragma once


#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <omp.h>
#include "../utility/tic_toc.h"
#include<fstream>
#include <mutex>

using namespace std;

namespace flower_focus {
typedef Eigen::Matrix<double, 1, 1, Eigen::ColMajor> Vector1d;
typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor> Vector6d;
typedef Eigen::Matrix<double, 7, 1, Eigen::ColMajor> Vector7d;
typedef Eigen::Matrix<double, 9, 1, Eigen::ColMajor> Vector9d;
typedef Eigen::Matrix<double, 15, 1, Eigen::ColMajor> Vector15d;
typedef Eigen::Matrix<double, 3, 1, Eigen::ColMajor> Vector3d;
typedef Eigen::Matrix<double, 2, 1, Eigen::ColMajor> Vector2d;
typedef Eigen::VectorXd VectorXd;


typedef  Eigen::Matrix<double, 2, 2, Eigen::RowMajor> Matrix2d;
typedef  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matrix3d;
typedef  Eigen::Matrix<double, 9, 9, Eigen::RowMajor> Matrix9d;
typedef  Eigen::Matrix<double, 6, 6, Eigen::RowMajor> Matrix6d;
typedef  Eigen::Matrix<double, 6, 3, Eigen::RowMajor> Matrix6_3d;
typedef  Eigen::Matrix<double, 2, 6, Eigen::RowMajor> Matrix2_6d;
typedef  Eigen::Matrix<double, 2, 3, Eigen::RowMajor> Matrix2_3d;
typedef   Eigen::Matrix<double, 15, 6, Eigen::RowMajor> Matrix15_6d;
typedef   Eigen::Matrix<double, 15, 3, Eigen::RowMajor> Matrix15_3d;
typedef   Eigen::Matrix<double, 15, 9, Eigen::RowMajor> Matrix15_9d;
typedef  Eigen::Matrix<double, Eigen::Dynamic, 9, Eigen::RowMajor> MatrixX9d;
typedef  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXd;
typedef  Eigen::Matrix<double, 1, 3, Eigen::RowMajor> Matrix1_3d;
typedef  Eigen::Matrix<double, 1, 6, Eigen::RowMajor> Matrix1_6d;
typedef  Eigen::Matrix<double, 15, 15, Eigen::RowMajor> Matrix15d;

}


#define ENABLE_STEREO_MAP 0
#define REAL_TIME 0
#define USE_CUDA 0
#define USE_STEREO 0


#if USE_STEREO
#define STEREO_SCALE 1
#define ESTIMATE_EXTRINSIC2 0
#else
#define ESTIMATE_EXTRINSIC2 0
#endif


#define ESTIMATE_SCALE 0

extern double spp_psr_weight;
extern double spp_phase_weight;
#if REAL_TIME
extern std::mutex mutex_pose;
#endif
extern bool strategy2;
//switch for enabling debug mode. More information provide.
#define DEBUG 1
#if ENABLE_STEREO_MAP
extern std::map<double, Eigen::Matrix3f>PubRs;
extern std::map<double, Eigen::Vector3f>Pubts;
#endif

#define NLQ_THRESHOLD 0.01
//eps for peudo inverse
#define EPSS 1e-8
//feature weight
#define SQRT_INFO_FEATURE  FOCAL_LENGTH / F_NOISE_SIGMA
//default feature depth
#define INIT_DEPTH 2.0
//Dogleg parameter
#define initial_radius 10
//Dogleg parameter
#define IDEPTH_REGION 1
//amplifing the imu noise in the start phase
#define AMP 2
//amplifing the imu noise in the start phase
#define AMP_NUM 4
//switching for compensating the gray changing
#define USE_GRAY_SCALE 1
//minimum feature observation
#define FEATURE_CONTINUE 2
//time duration for computing the initial roll and pitch angles
#define AVERAGE_TIME 0.1
//outlier rejection threshold :[OUTLIER_SHRESHOLD1*SQRT_INFO_FEATURE] pixel
#define OUTLIER_SHRESHOLD1 4
//number of layer for LK optical flow algorithm.
#define OPTICAL_LAYER_NUM 4
//ransac threshold
#define F_THRESHOLD 1.0
//maximum number of feature points
#define MAX_CNT 250
// size of the block
#define SWF_SIZE_IN 10

#define MIN_PARALLAX 10
//size of the SWF
#define SWF_WINDOW_SIZE (SWF_SIZE_OUT * SWF_SIZE_IN)

#define LEFT 1

#define RIGHT 0
// number of block, the capacity of the SWF is: SWF_SIZE_IN*SWF_SIZE_OUT
enum SolverFlag {
    //not initialized
    Initial,
    //have initialized
    NonLinear
};

extern Eigen::Matrix3d RC2C1;
extern Eigen::Vector3d TC2C1;
extern bool needNreaolve;
extern bool USE_RTK;
extern bool USE_RTD;
extern bool USE_GNSS;
extern bool USE_DOPPLER;
extern bool USE_SPP_PHASE;
extern bool USE_SPP_CORRECTION;
extern double last_fix_time;
extern bool USE_REPARAM;
extern std::string IMAGE1_TOPIC;
extern int SWF_SIZE_OUT;
extern Eigen::Matrix3d Rwgw;
extern Eigen::Vector3d ANCHOR_POINT;
extern Eigen::Matrix3d InitRwgw;
extern std::string RTK_TOPIC;
extern std::string MAG_TOPIC;
extern Eigen::Vector3d Pbg;
extern bool init_gnss;
extern bool rtk_fix;
extern Eigen::Vector3d base_pos;
extern double blackvalue;
extern SolverFlag solver_flag;
extern Eigen::Matrix3d R_WI_WC;
extern Eigen::Vector3d P_WI_WC;
extern int ESTIMATE_R_WI_WC;
extern bool ESTIMATE_P_WI_WC;
extern bool ESTIMATE_PBG;
extern bool ESTIMATE_TD;
extern bool ESTIMATE_ACC_SCALE;
extern double kMinMu;
extern bool boolsmallmuwhenmarge;
extern cv::Mat mask_global;
//maximum number of keyframe for outlier rejection
extern int FRAME_NUM_FOR_REJECTION;
//enabling output
extern bool enable_output;
//inverse depth pointer
extern std::set<double*>idepth_pointers;
//file name of the exposure time
extern std::string exposure_time_file;
//maximum number of visual frames
extern int ESTIMATE_EXTRINSIC;
//have performing marginalization before.
extern bool have_hist;
//inverse depth pointer, divided by blocks
extern std::vector<std::unordered_map<int, double>> idepths_all;
//
extern int LEAK_NUM;
//file storage for logging information
extern std::ofstream LOG_OUT;
//optical flow searching window size
extern int OPTICAL_WINSIZE;
//path for saving the estimation
extern std::string RESULT_PATH;
//imu ros topic
extern std::string IMU_TOPIC;
//bag path
extern std::string ROS_PATH;
//other topic
extern std::string IMAGE0_TOPIC;
//path of intrinc calibration
extern std::vector<std::string> CAM_NAMES;
//anchor point for calculating the position under ENU frame.
extern Eigen::Vector3d Pgb;
extern std::vector<Eigen::Quaterniond> QIC;
//rotation from camera to imu
extern std::vector<Eigen::Matrix3d> RIC;
//translation from camera to imu
extern std::vector<Eigen::Vector3d> TIC;
//gravity vector
extern Eigen::Vector3d G;
// imu noise standard deviation
extern double ACC_N, INT_N, GYR_NONLINEARITY, ACC_W, GYR_N, GYR_W;
//publish tracking image as topic
extern int SHOW_TRACK;
//perform forward and backward optical flow to improve feature tracking accuracy
extern int FLOW_BACK;
//min distance between two features
extern int MIN_DIST;
//focal length
extern double FOCAL_LENGTH;
//standard deviation of feature noise.
extern double F_NOISE_SIGMA;
//time duration for skipping the first SKIP_TIME [s]
extern double SKIP_TIME;
//the system could fix the scale now
extern bool fix_scale;
//the equipment has excitation before
extern bool has_excitation;
extern double max_continue_time;
extern double phase_bias_weight;
extern omp_lock_t omp_lock;


void readParameters(std::string config_file);
int localSize(int size);
int globalSize(int size);

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



enum SaveRestoreFlag {
    Restore = 0,
    Save

};



#if DEBUG
#define USE_ASSERT 0
//number of threads
#define NUM_THREADS 2
#else
#define USE_ASSERT 0
//number of threads
#define NUM_THREADS 2
#endif

#if USE_ASSERT
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#endif


class PosInfo {
  public:
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    double time_stamp;
};

