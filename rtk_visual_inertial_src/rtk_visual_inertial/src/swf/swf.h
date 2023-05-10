

#pragma once


#include <std_msgs/Header.h>
#include "ceres/ceres.h"
#include "ceres/schur_complement_solver.h"
#include <eigen3/Eigen/Dense>
#include "../parameter/parameters.h"
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../feature/feature_tracker.h"
#include "../feature/feature_manager.h"
#include "../factor/gnss_imu_factor.h"
#include "common_function.h"



class SWFOptimization {
  public:

    enum Optimization_mode {
        //normal optimization using ceres
        NormalMode = 0,
        //marginalizing the selected states using marginalization_factor(not efficient for large number set of states)
        MargeIncludeMode,
        //marginalizing states except the selected states
        MargeExcludeMode,
        //solving the Gauss-Newton problem and obtaining the hessian matrix in one-iteration (barely used)
        NormalMode2,
        //marginalizing the selected states using ceres(efficient for large number set of states)
        MargeIncludeMode2,
        //used for marginalizing the gnss frames
        GNSSMargIncludeMode
    };
    enum SolverFlag {
        //not initialized
        Initial,
        //have initialized
        NonLinear
    };
    enum MarginalizationFlag {
        //marginalizing the oldest visual frame and the gnss frame between the oldest visual and second-oldest visual frame.
        MargImagOld = 0,
        //removing the second newest visual frame
        MargImagSecondNew = 1,
        //not supported in current version
        MargRoverOld = 2,
        //not supported in current version
        MargRoverSecondNew = 3,
        //error mode
        MargNull = 4
    };
    enum FrameType {
        //the frame is gnss frame
        GnssFrame = 5,
        // the frame is visual frame
        ImagFrame = 6,
        //error occur
        ErroFrame = 7
    };
    enum Direction {
        //transform visual frame index and gnss frame index to frame index
        IR2Frame = 0,
        //transform frame index to visual frame index and gnss frame index
        Frame2IR = 1,
    };

    SWFOptimization();
    //seting the global setting
    void SetParameter();
    //adding imu measurements to buffer
    void InputIMU(double t, const Eigen::Vector3d& linearAcceleration, const Eigen::Vector3d& angularVelocity);
    //adding visual measurements to buffer
    void InputImage(double t, const cv::Mat& _img, const cv::Mat& _img1 = cv::Mat());
    //adding magnetometer measurements to buffer
    void InputMag(double t, const Eigen::Vector3d& Magdata);
    //adding gnss measurements to buffer
    void InputGnss(mea_t* rovernew2);
    //measurement process, major loop
    void MeasurementProcess();
    //mag process
    void MagProcess(double curtime);
    //imu process; imu pre-integration
    void IMUProcess( double dt, const Eigen::Vector3d& linear_acceleration, const Eigen::Vector3d& angular_velocity);
    //image postprocess; outlier rejection
    void ImagePostprocess();
    //image preprocess; adding feature points
    void ImagePreprocess(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image);
    //gnss processing
    void GnssProcess(mea_t* rover);
    //gnss preprocessing
    void GnssPreprocess(mea_t* data);
    //vector to optimization pointer
    void Vector2Double();
    //optimization pointer to vector
    void Double2Vector();
    //sliding out second-newest visual frame
    void SlideWindowNew();
    //sliding out oldest visual frame
    void SlideWindowOld (Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0, Eigen::Vector3d Pbg);
    //sliding out visual frames
    void SlideWindowFrame(int frameindex, int, bool updateIMU);
    //optimizing states
    void MyOptimization( );
    //my predefined elimination order.
    void MyOrdering(ceres::Problem& problem, ceres::Solver::Options& options);
    //adding all the residual to the problem
    void AddAllResidual(int mode, std::set<double*>MargePoint, MarginalizationInfo* marginalization_info, ceres::Problem& problem,
                        ceres::Solver::Options& options, bool usefeature, bool useimu, bool usegnss);
    //addding the gnss residual to the problem
    void AddGnssResidual(int mode, std::set<double*>MargePoint, MarginalizationInfo* marginalization_info, ceres::Problem& problem,
                         ceres::Solver::Options& options, int ir);
    //initialize the constraint of the first frame
    void InitializeSqrtInfo();
    //initialize the orintation,position bias.
    void InitializePos(Eigen::Matrix3d& Rwgw_);
    //reset the states//not complete;need to fix
    void ClearState();
    //marginzalized the middle gnss frame
    void MiddleMargGnssFrame();
    //marginalizing the selected frames
    void MargFrames(std::set<int> margeindex);
    //performing marginalization using the global ceres problem
    void GlobalMarge(std::set<double*>MargePoint);
    //marginalized the selected gnss frames
    MarginalizationInfo* MargGNSSFrames(std::set<int> margeindex, IMUGNSSBase* IMUGNSSmeasurement);
    //obtaing states for marginalization
    std::set <double*> FindMargSet(std::set<int> margeindex);
    //obtaing imu measurements
    bool GetImuInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>>& accVector,
                        vector<pair<double, Eigen::Vector3d>>& gyrVector);
    //imu measurements avaliable
    bool ImuAvailable(double t);
    //imu pre-integration
    void ImuIntegrate();
    //lambda searching for integer resolution
    void LambdaSearch();
    //prepare for integer resolution
    void IntegerSolve();
    //obtaining the heesian matrix and right hand side vector(from ceres solver)
    void UpdateSchur(ceres::Problem& problem);
    //obtaining the heesian matrix only (from ceres solver)
    void UpdateSchurHessianOnly(ceres::Problem& problem);
    //reject the visual outliers
    void OutliersRejection(set<int>& removeIndex, int image_count);
    //calculating the projection error
    double ReprojectionError(Eigen::Matrix3d& Rj, Eigen::Vector3d& Pj, Eigen::Matrix3d& ricj, Eigen::Vector3d& ticj, Eigen::Vector3d pts_w, Eigen::Vector3d& uvj);
    //adding parameter function to the ceres problem
    void AddParameter2Problem(ceres::Problem& problem, bool disable_constant);
    //adding feature points to the global optimization
    void AddFeature2Problem();
    //checking whether the status of ceres problem is as expection
    int CheckProblem();
    //saving the phase biases and seting them to zero
    void PhaseBiasSaveAndReset();
    //restore the values of phase biases
    void PhaseBiasRestore();
    //adding the gnss and imu measurements to the gnss-imu factor
    void UpdateImuGnssFactor();
    //reseting the imu-gnss factor//not efficiency; need to fixed
    void ResetImuGnssFactor(int i, MarginalizationInfo* marginfo);
    //pub estimation for visualization and storing
    void PubData();
    //updating the visual and gnss index
    int ImageRoverId2FrameId(int image_rover_count, int mode);
    //updating the visual and gnss index (i2f,g2f,f2g,f2i)
    void UpdateVisualGnssIndex();
    //copy the states of visual frame to the states of frame; or reverse
    void transformdataimage(Direction direction, int mode);
    //slide out the visual frame
    void SlideWindow();
    //prepare pointer for obtaining their hessian matrix and right hand side vector.
    //the hessian matrix and right hand side vector is obtained when solving the Gauss-Newton problem
    void UpdateNParameterHead();


    //extrinsic rotation
    Matrix3d ric[2];
    //extrinsic translation
    Vector3d tic[2];
    //extrinsic parameter
    double para_ex_Pose[2][SIZE_POSE];
    //magnetometer bias
    double para_bmg[3];
    //position of the frames
    Vector3d        Ps[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //velocity of the frames
    Vector3d        Vs[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //rotation of the frames
    Matrix3d        Rs[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //acc bias of the frames
    Vector3d        Bas[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //gyr bias of the frames
    Vector3d        Bgs[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //position of the visual frames
    Vector3d        Ps_image[(FEATURE_WINDOW_SIZE  + 1)];
    //velocity of the visual frames
    Vector3d        Vs_image[(FEATURE_WINDOW_SIZE  + 1)];
    //rotation of the visual frames
    Matrix3d        Rs_image[(FEATURE_WINDOW_SIZE  + 1)];
    //time stamp of the frames
    double headers[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //position and orientation of the frames
    double* para_pose[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)] = {0}; //[SIZE_POSE]
    //velocity, gyr bias and acc bias of the frames
    double* para_speed_bias[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)] = {0}; //[SIZE_SPEEDBIAS]
    //receiver clock parameter of the frames
    double* para_gnss_dt[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)] = {0}; //[6]
    //visual frame index to frame index
    int i2f[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //gnss frame index to frame index
    int g2f[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //frame index to visual frame index
    int f2i[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //frame index to gnss frame index
    int f2g[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //categories of the frames.
    FrameType frame_types[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //accerometer observation(not being added yet)
    queue<pair<double, Eigen::Vector3d>> acc_buf;
    //gyr observation(not being added yet)
    queue<pair<double, Eigen::Vector3d>> gyr_buf;
    //magnetometer observation(not being added yet)
    std::list<pair<double, Eigen::Vector3d>> mag_buf;
    //feature point observation(not being added yet)
    queue<pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1> > > > > > feature_buf;
    //gnss observation(not being added yet)
    std::list<mea_t*>rover_buf;
    //IMU delta time(being added)
    vector<double> dt_buf[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //accerometer observation(being added)
    vector<Vector3d> linear_acceleration_buf[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //gyr observation(being added)
    vector<Vector3d> angular_velocity_buf[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];
    //magnetometer observation(being added)
    Vector3d mags[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //gnss observation(being added)
    mea_t* rovers[FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1];
    //IMU pre-integration measurements
    IntegrationBase* pre_integrations[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)] = {0};
    //IMU-GNSS factors
    IMUGNSSBase* imu_gnss_factor[FEATURE_WINDOW_SIZE + 2];
    //global ceres problem;avaliable only when USE_GLOBAL_OPTIMIZATION is 1
    ceres::Problem my_problem;
    //global ceres option;avaliable only when USE_GLOBAL_OPTIMIZATION is 1
    ceres::Solver::Options my_options;
    //prior residual id in ceres ;avaliable only when USE_GLOBAL_OPTIMIZATION is 1
    ceres::ResidualBlockId marg_residual_block_id = 0;
    //average acc measurements, used for initialize the roll and pitch angles
    Vector3d acc_mean;
    //average gyr measurements, used for initialize the gyr bias
    Vector3d gyr_mean;
    //average magnetnometer measurements, used for initialize the yaw angle
    Vector3d mag_mean;
    //hessian matrix from ceres solver
    Eigen::MatrixXd A;
    //right hand side from ceres solver
    Eigen::VectorXd b;
    //parameter for obtained the corresponding  hessian matrix when solving the Gauss-Newton problem
    std::vector<const double*>parameter_block_addr;
    //parameter size of parameter_block_addr
    std::vector<int>parameter_block_global_size;
    //number of image frame in SWF
    int image_count;
    //number of gnss frame in SWF
    int rover_count;
    //optical feature tracker
    FeatureTracker feature_tracker;
    //feature manager for hadling the feature points
    FeatureManager f_manager;
    //storing the RB-SD phase biases
    std::list<PBtype> rtk_phase_bias_variables[MAXSATNUM * 2];
    //storing the rover-only phase biases
    std::list<PBtype> spp_phase_bias_variables[MAXSATNUM * 2];
    //storing the pseudorange correction
    std::list<PBtype> pseudorange_correction_variables[MAXSATNUM * 2];
    //storing the phase bias values
    std::vector<double>phase_bias_save;
    //no use variables;use to avoid bug
    double blackvalue = 0, blackvalue2 = 0;
    //prior factor
    MarginalizationInfo* last_marg_info = 0;
    //which types of frames should be marginalized. only should be MargImagOld and MargImagSecondNew
    MarginalizationFlag  marg_flag;
    //solving frag. Initial when the system need a initialization. Only the visual-inertial only sysytem required a initialization process.
    SolverFlag solver_flag;
    //base position
    Eigen::Vector3d base_pos;
    //rotation from ECEF frame to ENU frame at the anchor point.
    Eigen::Matrix3d InitRwgw;
    //last accerameter and gyr measureements
    Vector3d acc_0, gyr_0;
    //last update time stamp and current time stamp
    double prev_time, prev_time2, cur_time;
    //0 when need a gnss initialization process
    bool init_gnss;
    //0 when need a imu initialization process
    bool imu_initialize;
    //0 when need a initialization for visualization
    bool pub_init;
    //extimate camera-imu extrinsic calibration.
    bool open_ex_estimation;
    //need the first imu measurement for initialize
    bool first_imu;
    //visual frame index for marginalization
    int imag_marg_index;
    //first fix solution avaliable
    bool fix;
    //number of imu measurements (for imu initialization)
    int gyr_count = 0;
    //number of fixed solution
    int gnss_fix_solution_count;
    //current fix solution avaliable
    bool rtk_fix;
    //time stamp when the last fixed solution is avaliable
    double last_fix_time = 0;
    //yaw estimation from magnetometer
    double mag_yaw;
    //number of epoch that fixed solutions are unavaliable
    int not_fix_count;
    //last epoch is fixed
    bool last_fix = false;
    //number of rover frames in the history
    int rover_count_accumulate = 0;
    //corrent frame required a integer resolution
    bool need_Nresolve = false;


};





