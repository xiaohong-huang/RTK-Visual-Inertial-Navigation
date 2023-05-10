
#include "parameters.h"
#include <ros/ros.h>
#include "../utility/utility.h"
#include <fstream>
#include <map>

double MIN_PARALLAX;
double ACC_N, ACC_W;
double GYR_N, GYR_W;
double MAX_SOLVER_TIME;
double F_THRESHOLD;
double USE_FEATURE;
double AVERAGE_IMU;
double SKIP_IMU;
double AVERAGE_IMAGE;
double MAX_TRUST_REGION_RADIUS;


std::vector<Eigen::Matrix3d> RIC;
std::vector<Eigen::Vector3d> TIC;
Eigen::Vector3d G{0.0, 0.0, 9.8};
Eigen::Matrix3d MagMatrix;
Eigen::Matrix3d IMUMatrix;
Eigen::Vector3d Pbg;
Eigen::Vector3d ANCHOR_POINT;
Eigen::Vector3d MagVector;
Eigen::Matrix3d Rwgw;

int MAX_NUM_ITERATIONS;
int ESTIMATE_EXTRINSIC;
int NUM_OF_CAM;
int MAX_CNT;
int MIN_DIST;
int SHOW_TRACK;
int FLOW_BACK;
int CARRIER_PHASE_CONTINUE_THRESHOLD;
int FIX_CONTINUE_THRESHOLD;
int Phase_ALL_RESET_COUNT;

std::string EX_CALIB_RESULT_PATH;
std::string RESULT_PATH;
std::string IMU_TOPIC;
std::string RTK_TOPIC;
std::string MAG_TOPIC;
std::string PPS_LOCAL_TOPIC;
std::string PPS_GPS_TOPIC;
std::string FEATURE_TOPIC;
std::string IMAGE0_TOPIC, IMAGE1_TOPIC;
std::vector<std::string> CAM_NAMES;
std::string ROS_PATH;


bool USE_STEREO;
bool USE_IMAGE;
bool USE_GNSS;
bool USE_IMU;
bool USE_RTK;
bool USE_RTD;
bool USE_MAG_INIT_YAW;
bool USE_MAG_CORRECT_YAW;
bool USE_DOPPLER;
bool USE_SPP_PHASE;
bool USE_SPP_CORRECTION;
bool USE_GLOBAL_OPTIMIZATION;
bool USE_DIRECT_N_RESOLVE;
bool USE_N_RESOLVE;












void readParameters(std::string config_file) {

    FILE* fh = fopen(config_file.c_str(), "r");
    if (fh == NULL) {
        printf("config_file dosen't exist; wrong config_file path");
        // abort();
        abort();
        return;
    }
    fclose(fh);

    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    if (!fsSettings.isOpened()) {
        std::cerr << "ERROR: Wrong path to settings" << std::endl;
    }

    fsSettings["USE_IMAGE"] >> USE_IMAGE;
    fsSettings["USE_GNSS"] >> USE_GNSS;
    fsSettings["USE_IMU"] >> USE_IMU;
    fsSettings["USE_RTK"] >> USE_RTK;
    fsSettings["USE_RTD"] >> USE_RTD;
    fsSettings["USE_MAG_INIT_YAW"] >> USE_MAG_INIT_YAW;
    fsSettings["USE_MAG_CORRECT_YAW"] >> USE_MAG_CORRECT_YAW;
    fsSettings["USE_SPP_PHASE"] >> USE_SPP_PHASE;
    fsSettings["USE_SPP_CORRECTION"] >> USE_SPP_CORRECTION;
    fsSettings["USE_DOPPLER"] >> USE_DOPPLER;
    fsSettings["USE_GLOBAL_OPTIMIZATION"] >> USE_GLOBAL_OPTIMIZATION;
    fsSettings["AVERAGE_IMU"] >> AVERAGE_IMU;
    fsSettings["SKIP_IMU"] >> SKIP_IMU;
    fsSettings["AVERAGE_IMAGE"] >> AVERAGE_IMAGE;
    fsSettings["CARRIER_PHASE_CONTINUE_THRESHOLD"] >> CARRIER_PHASE_CONTINUE_THRESHOLD;
    fsSettings["FIX_CONTINUE_THRESHOLD"] >> FIX_CONTINUE_THRESHOLD;
    fsSettings["Phase_ALL_RESET_COUNT"] >> Phase_ALL_RESET_COUNT;
    fsSettings["USE_DIRECT_N_RESOLVE"] >> USE_DIRECT_N_RESOLVE;
    fsSettings["USE_N_RESOLVE"] >> USE_N_RESOLVE;
    fsSettings["MAX_TRUST_REGION_RADIUS"] >> MAX_TRUST_REGION_RADIUS;
    fsSettings["USE_STEREO"] >> USE_STEREO;
    fsSettings["image0_topic"] >> IMAGE0_TOPIC;
    fsSettings["image1_topic"] >> IMAGE1_TOPIC;
    fsSettings["rtk_potic"] >> RTK_TOPIC;
    fsSettings["mag_potic"] >> MAG_TOPIC;
    fsSettings["pps_local_topic"] >> PPS_LOCAL_TOPIC;
    fsSettings["pps_gps_topic"] >> PPS_GPS_TOPIC;
    fsSettings["feature_topic"] >> FEATURE_TOPIC;
    fsSettings["imu_topic"] >> IMU_TOPIC;

    MAX_SOLVER_TIME = fsSettings["MAX_SOLVER_TIME"];
    MAX_NUM_ITERATIONS = fsSettings["MAX_NUM_ITERATIONS"];
    MIN_PARALLAX = fsSettings["keyframe_parallax"];
    MIN_PARALLAX = MIN_PARALLAX / FOCAL_LENGTH;
    MAX_CNT = fsSettings["max_cnt"];
    MIN_DIST = fsSettings["min_dist"];
    F_THRESHOLD = fsSettings["F_THRESHOLD"];
    USE_FEATURE = fsSettings["USE_FEATURE"];
    SHOW_TRACK = fsSettings["SHOW_TRACK"];
    FLOW_BACK = fsSettings["FLOW_BACK"];
    ACC_N = fsSettings["acc_n"];
    ACC_W = fsSettings["acc_w"];
    GYR_N = fsSettings["gyr_n"];
    GYR_W = fsSettings["gyr_w"];
    G.z() = fsSettings["g_norm"];
    ESTIMATE_EXTRINSIC = fsSettings["ESTIMATE_EXTRINSIC"];

    cv::Mat cv_M;
    fsSettings["Mag_Matrix"] >> cv_M;
    Eigen::Matrix4d M;
    cv::cv2eigen(cv_M, M);
    MagMatrix = M.block<3, 3>(0, 0);
    MagVector = M.block<3, 1>(0, 3);

    fsSettings["IMU_Matrix"] >> cv_M;
    cv::cv2eigen(cv_M, M);
    IMUMatrix = M.block<3, 3>(0, 0);

    cv::Mat cv_V;
    Eigen::Vector3d V;

    fsSettings["Pbg"] >> cv_V;
    cv::cv2eigen(cv_V, V);
    Pbg = V;

    fsSettings["ANCHOR_POINT"] >> cv_V;
    cv::cv2eigen(cv_V, V);
    ANCHOR_POINT = V;


    std::cout << "result path " << RESULT_PATH << std::endl;




    if (ESTIMATE_EXTRINSIC == 2) {
        printf("have no prior about extrinsic param, calibrate extrinsic param");
        RIC.push_back(Eigen::Matrix3d::Identity());
        TIC.push_back(Eigen::Vector3d::Zero());
        EX_CALIB_RESULT_PATH = "/media/huang/A/extrinsic_parameter.csv";
    } else {
        if ( ESTIMATE_EXTRINSIC == 1) {
            printf(" Optimize extrinsic param around initial guess!");
            EX_CALIB_RESULT_PATH =  "/media/huang/A/extrinsic_parameter.csv";
        }
        if (ESTIMATE_EXTRINSIC == 0)
            printf(" fix extrinsic param ");

        cv::Mat cv_T;
        fsSettings["body_T_cam0"] >> cv_T;
        Eigen::Matrix4d T;
        cv::cv2eigen(cv_T, T);
        RIC.push_back(T.block<3, 3>(0, 0));
        TIC.push_back(T.block<3, 1>(0, 3));
    }
    NUM_OF_CAM = fsSettings["num_of_cam"];
    printf("camera number %d\n", NUM_OF_CAM);

    if (NUM_OF_CAM != 1 && NUM_OF_CAM != 2) {
        printf("num_of_cam should be 1 or 2\n");
        assert(0);
    }


    int pn = config_file.find_last_of('/');
    std::string configPath = config_file.substr(0, pn);

    std::string cam0Calib;
    fsSettings["cam0_calib"] >> cam0Calib;
    std::string cam0Path = configPath + "/" + cam0Calib;
    CAM_NAMES.push_back(cam0Path);

    if (NUM_OF_CAM == 2) {
        std::string cam1Calib;
        fsSettings["cam1_calib"] >> cam1Calib;
        std::string cam1Path = configPath + "/" + cam1Calib;
        //printf("%s cam1 path\n", cam1Path.c_str() );
        CAM_NAMES.push_back(cam1Path);

        cv::Mat cv_T;
        fsSettings["body_T_cam1"] >> cv_T;
        Eigen::Matrix4d T;
        cv::cv2eigen(cv_T, T);
        RIC.push_back(T.block<3, 3>(0, 0));
        TIC.push_back(T.block<3, 1>(0, 3));
    }







    fsSettings.release();
}
