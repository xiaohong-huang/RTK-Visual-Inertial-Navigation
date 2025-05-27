

#pragma once


#include <std_msgs/Header.h>
#include <eigen3/Eigen/Dense>
#include "../parameter/parameters.h"
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../feature/feature_tracker.h"
#include "../feature/feature_manager.h"
#include "../factor/gnss_imu_info.h"
#include "../factor/integration_base.h"
#include "../factor/visual_inerial_factor.h"
#include "../gnss/common_function.h"

class SWFOptimization {
  public:

    enum Optimization_mode {
        NormalMode = 0,
        MargeMode,
    };

    enum MarginalizationFlag {
        //marginalizing the oldest visual frame and the gnss frame between the oldest visual and second-oldest visual frame.
        MargImagOld = 0,
        //removing the second newest visual frame
        MargImagSecondNew = 1,
        //error mode
        MargNull = 4
    };


    SWFOptimization();
    //seting the global setting
    void SetParameter();
    //adding imu measurements to buffer
    void InputIMU(double t, const Eigen::Vector3d& linearAcceleration, const Eigen::Vector3d& angularVelocity);
    //adding visual measurements to buffer
    void InputImage(double t, const cv::Mat& _img, const cv::Mat& _img1 = cv::Mat());
    //measurement process, major loop
    void MeasurementProcess();
    //imu process; imu pre-integration
    void IMUProcess( double dt, const Eigen::Vector3d& linear_acceleration, const Eigen::Vector3d& angular_velocity, IntegrationBase*& pre_integration);
    //image postprocess; outlier rejection
    void ImagePostprocess();
    //image preprocess; adding feature points
    void ImagePreprocess(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image);
    //vector to optimization pointer
    void Vector2Double();
    //optimization pointer to vector
    void Double2Vector();
    //sliding out second-newest visual frame
    void SlideWindowNew(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0);
    //sliding out oldest visual frame
    void SlideWindowOld (Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0);
    //sliding out visual frames
    void SlideWindowFrame(int frameindex, int, bool updateIMU);
    //optimizing states
    void MyOptimization( );
    //initialize the constraint of the first frame
    void InitializeSqrtInfo();
    //initialize the orintation,position bias.
    void InitializePose();
    //reset the states//not complete;need to fix
    void ClearState();
    //marginalizing the selected frames
    void MargFrames();
    //obtaing imu measurements
    bool GetImuInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>>& accVector,
                        vector<pair<double, Eigen::Vector3d>>& gyrVector);
    //imu measurements avaliable
    bool ImuAvailable(double t);
    //imu pre-integration
    void ImuIntegrate();
    //reject the visual outliers
    void OutliersRejection( int image_count);
    //calculating the projection error
    double ReprojectionError(Eigen::Matrix3d& Rj, Eigen::Vector3d& Pj, Eigen::Matrix3d& ricj, Eigen::Vector3d& ticj, Eigen::Vector3d pts_w, Eigen::Vector3d& uvj, double weight);
    //pub estimation for visualization and storing
    void PubData();
    //slide out the visual frame
    void SlideWindow();
    //optimination or marginalization function
    void OptimizationOrMarginalization( int mode);
    //Adding factors to the SWF
    void AddFactors(std::vector<VisualInertialBase*>& visual_inertial_bases, std::set<int>factor_index_set);
    //Adding features to the SWF.
    void AddFeatures(FeaturePerId& it_per_id, std::vector<VisualInertialBase*>& visual_inertial_bases, std::set<int>& factor_index_set);
    //reject the tail of feature measurements and regenerate a new feature
    void RegenerateFeatures(FeaturePerId& it_per_id, int remain_size, std::vector<FeaturePerId>& new_features);
    //
    void UpdataGlobalFactors();
    //generating the position of all the visual frames.
    std::vector<PosInfo> RetriveAllPose();
    //saving position of the visual keyframes
    void SaveKefPos(int);
    //saving position of the visual non-keyframes.
    void SaveLocalPos(int, int);

    //position of the frames
    std::vector<Vector3d>        Ps;//ECEF if GNSS,IMU frame
    //velocity of the frames
    std::vector<Vector3d>        Vs;//ECEF if GNSS,IMU frame
    //rotation of the frames
    std::vector<Matrix3d>        Rs;//ECEF if GNSS,IMU frame
    //acc bias of the frames
    std::vector<Vector3d>        Bas;
    //gyr bias of the frames
    std::vector<Vector3d>        Bgs;
    //time stamp of the frames
    std::vector<double> headers;
    //position and orientation of the frames
    std::vector<double*> para_pose;
    //velocity, gyr bias and acc bias of the frames
    std::vector<double*> para_speed_bias;
    //imu-camera time shift.
    double td = 0;
    //imu-camera extrinsic calibration
    double para_extrinsic[7] = {0};
    //global orientation, only the tail 4 parameters are used
    double para_global_pos[7] = {0};
    //accerometer observation(not being added yet)
    queue<pair<double, Eigen::Vector3d>> acc_buf;
    //gyr observation(not being added yet)
    queue<pair<double, Eigen::Vector3d>> gyr_buf;
    //magnetometer observation(not being added yet)
    std::list<pair<double, Eigen::Vector3d>> mag_buf;
    //feature point observation(not being added yet)
    queue<pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>>> feature_buf;
    std::vector<GNSSIMUInfo*> gnss_imu_infos;

    //average acc measurements, used for initialize the roll and pitch angles
    Vector3d acc_mean;
    //number of image frame in SWF
    int image_count;
    //optical feature tracker
    FeatureTracker feature_tracker;
    //feature manager for hadling the feature points
    FeatureManager f_manager;
    //prior factor
    MarginalizationInfo* last_marg_info = 0;
    //which types of frames should be marginalized. only should be MargImagOld and MargImagSecondNew
    MarginalizationFlag  marg_flag;
    //solving frag. Initial when the system need a initialization. Only the visual-inertial only sysytem required a initialization process.

    //last accerameter and gyr measureements
    Vector3d acc_0, gyr_0;
    //last update time stamp and current time stamp
    double prev_time, cur_time;
    //0 when need a initialization for visualization
    bool pub_init = false;
    //visual frame index for marginalization
    int imag_marg_index;
    //number of imu measurements (for imu initialization)
    int acc_count = 0;
    //visual scale factor
    double scale_factor = 1;
    //time stamp for the first measurement
    double first_observe_time = 0;
    //for global optimization.
    std::vector<VisualInertialBase*> visual_inertial_bases_global;
    //the last time shift
    double old_time_shift = 0;
    //position of all the visual frames.
    std::map<double*, std::vector<PosInfo>> pos_save;
    //position of all the visual frames.
    std::vector<PosInfo> pos_save_all;
    //accerometor scale
    Eigen::Vector3d acc_scale;
    //gyroscope scale
    Eigen::Vector3d gyr_scale;


    void InputGnss(mea_t* rovernew2);
    std::list<mea_t*>rover_buf;
    void InputMag(double t, const Vector3d& Magdata);
    Vector3d mag_mean;
    int getVariableUseSingleNum(mea_t* obs_data);
    void set_last_marg_info_idx(int);
    Eigen::Vector3d old_P;
#if REAL_TIME
    std::mutex mutex_imu;
    std::mutex mutex_imu_feature;
    std::mutex mutex_feature;
#endif
#if USE_STEREO
    double para_extrinsic2[7] = {0};
#endif
};





