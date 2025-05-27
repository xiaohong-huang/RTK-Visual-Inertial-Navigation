#include <csignal>
#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "swf/swf.h"
#include "parameter/parameters.h"
#include "utility/visualization.h"
#include "utility/utility.h"
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <std_msgs/ByteMultiArray.h>
#include <sensor_msgs/MagneticField.h>
#include <sensor_msgs/NavSatFix.h>
#include<random>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>
#include "utility/depth_map.h"

SWFOptimization* swf_optimization;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::ImageConstPtr> img0_buf;
queue<sensor_msgs::ImageConstPtr> img1_buf;
double start_timestamp = 0;//1896.28
TicToc system_time;
double last_system_time;
#if REAL_TIME
std::mutex m_buf;
#endif
#include "gnss/common_function.h"
void mGNSS_callback(const std_msgs::ByteMultiArrayPtr& gnss_msg) {


    mea_t gnss_obs;
    char* rover_c = (char*)(&gnss_obs);
    memset(rover_c, 0, sizeof(mea_t));
    for (int i = 0; i < (int)(gnss_msg->data.size()); i++)
        rover_c[i] = gnss_msg->data[i];
    for (int i = 0; i < gnss_obs.obs_count; i++) {
        ObsMea* d = gnss_obs.obs_data + i;
        d->RTK_L[1] = d->SPP_L[1] = 0;
        d->RTK_P[1] = d->SPP_P[1] = 0;

    }

    mea_t* rovernew = new (mea_t);
    memcpy(rovernew, &gnss_obs, sizeof(mea_t));
    if (rovernew->ros_time < start_timestamp)return;

    swf_optimization->InputGnss(rovernew);


}


void mag_callback(const sensor_msgs::MagneticFieldPtr& mag_msg) {
    if (mag_msg->header.stamp.toSec() < start_timestamp)return;
    double t = mag_msg->header.stamp.toSec();
    double mx = mag_msg->magnetic_field.x * 1e6;
    double my = mag_msg->magnetic_field.y * 1e6;
    double mz = mag_msg->magnetic_field.z * 1e6;
    Vector3d mag(mx, my, mz);
    swf_optimization->InputMag(t, mag);

    return;
}


cv::Mat getImageFromMsg(const sensor_msgs::ImageConstPtr& img_msg) {
    cv_bridge::CvImageConstPtr ptr;
    if (img_msg->encoding == "8UC1") {
        sensor_msgs::Image img;
        img.header = img_msg->header;
        img.height = img_msg->height;
        img.width = img_msg->width;
        img.is_bigendian = img_msg->is_bigendian;
        img.step = img_msg->step;
        img.data = img_msg->data;
        img.encoding = "mono8";
        ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
    } else
        ptr = cv_bridge::toCvCopy(img_msg, sensor_msgs::image_encodings::MONO8);
    cv::Mat img = ptr->image.clone();
    return img;
}


void processoneimage() {
    cv::Mat image0, image1;
    double time = 0;

#if USE_STEREO
    {
        if (!img0_buf.empty() && !img1_buf.empty()) {
            double time0 = img0_buf.front()->header.stamp.toSec();
            double time1 = img1_buf.front()->header.stamp.toSec();
            if (fabs(time0 - time1) < 1e-3)
                time0 = time1;
            if (time0 < time1) {
                time = img0_buf.front()->header.stamp.toSec();
                image0 = getImageFromMsg(img0_buf.front());
                img0_buf.pop();
            } else if (time0 > time1) {
                img1_buf.pop();
                printf("throw img1\n");
            } else {
                time = img0_buf.front()->header.stamp.toSec();
                image0 = getImageFromMsg(img0_buf.front());
                img0_buf.pop();
                image1 = getImageFromMsg(img1_buf.front());
                img1_buf.pop();
            }
        }
    }
#else
    {
        if (!img0_buf.empty() ) {
            time = img0_buf.front()->header.stamp.toSec();
            image0 = getImageFromMsg(img0_buf.front());
            img0_buf.pop();
        }
    }
#endif

    if (!image0.empty())
        swf_optimization->InputImage(time, image0, image1);
}


void img0_callback(const sensor_msgs::ImageConstPtr& img_msg) {
    if (img_msg->header.stamp.toSec() < start_timestamp)
        return;
    img0_buf.push(img_msg);
#if ENABLE_STEREO_MAP
    stereo_img0_callback(img_msg);
#endif
    processoneimage();
}
void img1_callback(const sensor_msgs::ImageConstPtr& img_msg) {
    if (img_msg->header.stamp.toSec() < start_timestamp)
        return;
#if ENABLE_STEREO_MAP
    stereo_img1_callback(img_msg);
#endif
#if USE_STEREO
    img1_buf.push(img_msg);
    processoneimage();
#endif
}


void imu_callback(const sensor_msgs::ImuConstPtr& imu_msg) {

    if (imu_msg->header.stamp.toSec() < start_timestamp - 1)
        return;
    double t = imu_msg->header.stamp.toSec();
    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Vector3d acc(dx, dy, dz);
    Vector3d gyr(rx, ry, rz);
    acc = acc;
    gyr = gyr;
    static double old_dt;
    if (old_dt != t)
        swf_optimization->InputIMU(t, acc, gyr);
    old_dt = t;
    return;
}





void bind_cpu(std::vector<int>cpu_set) {

}



void sig_handler( int sig ) {
    std::cout << "\tabort_\r\n";
    exit( 0 );
}




#if REAL_TIME

void MeasurementProcess() {
    while (1) {
        swf_optimization->MeasurementProcess();
        if (swf_optimization->feature_buf.empty()) {
            std::chrono::milliseconds dura(1);
            std::this_thread::sleep_for(dura);
        }
    }

}

int main(int argc, char** argv) {


    string config_file = argv[1];
    ROS_PATH = argv[2];
    RESULT_PATH = argv[3];
    printf("config_file: %s\n", argv[1]);
    readParameters(config_file);
#if 1
    ros::init(argc, argv, "flower_focus");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    registerPub(n);
#endif
    LEAK_NUM = SWF_SIZE_IN - 1;


#if DEBUG
    time_t t = time(nullptr);
    struct tm* now = localtime(&t);
    std::stringstream timeStr;
    timeStr << now->tm_year + 1900 << "-";
    timeStr << now->tm_mon + 1 << "-";
    timeStr << now->tm_mday << " ";
    timeStr << now->tm_hour << ":";
    timeStr << now->tm_min << ":";
    timeStr << now->tm_sec << ".log";
    std::string LOG_PATH = "log/" + RESULT_PATH + "--" + timeStr.str();
    LOG_OUT = std::ofstream(LOG_PATH, std::ios::out);
    LOG_OUT.precision(10);
#endif


    swf_optimization = new SWFOptimization();
    swf_optimization->SetParameter();


    printf("waiting for image and imu...");
    signal( SIGINT, sig_handler );
    std::cout << IMU_TOPIC << std::endl;
    std::cout << IMAGE0_TOPIC << std::endl;
#if ENABLE_STEREO_MAP
    stereo_init(config_file, n);
#endif

    ros::Subscriber sub_imu = n.subscribe("/" + IMU_TOPIC, 5000, imu_callback);
    ros::Subscriber sub_image = n.subscribe("/" + IMAGE0_TOPIC, 5000, img0_callback);
#if ENABLE_STEREO_MAP||USE_STEREO
    ros::Subscriber sub_image2 = n.subscribe("/" + IMAGE1_TOPIC, 5000, img1_callback);
#endif
    ros::Subscriber sub_mag;
    ros::Subscriber sub_gnss;

    if (RTK_TOPIC != "")
        sub_gnss = n.subscribe("/" + RTK_TOPIC, 5000, mGNSS_callback);
    if (MAG_TOPIC != "")
        sub_mag = n.subscribe( "/" + MAG_TOPIC, 5000, mag_callback);

    std::thread measurement_process {MeasurementProcess};
#if ENABLE_STEREO_MAP
    std::thread stereo_process {process_stereo};
#endif
#if 1
    ros::spin();
#else
    //As same as ORB-SLAM3, we generate the final resutls for evaluation.
    std::vector<PosInfo> camera_poses = swf_optimization->RetriveAllPose();
    ofstream foutC(RESULT_PATH, ios::out);
    foutC.setf(ios::fixed, ios::floatfield);

    double scale_factor = swf_optimization->scale_factor;
    static double old_time = 0;
    for (int i = 0; i < (int)camera_poses.size(); i++) {
        Eigen::Matrix3d Ri = R_WI_WC * camera_poses[i].R * RIC[0].transpose();
        Eigen::Vector3d ti = R_WI_WC * camera_poses[i].t * scale_factor - Ri * TIC[0] - Ri * Pgb + P_WI_WC;
        Eigen::Quaterniond Qi = Eigen::Quaterniond(Ri);
        double time_stamp = camera_poses[i].time_stamp;
        if (old_time)assert(time_stamp > old_time);

        foutC.precision(9);
        foutC << time_stamp << " ";
        foutC.precision(5);
        foutC << ti.x() << " "
              << ti.y() << " "
              << ti.z() << " "
              << Qi.x() << " "
              << Qi.y() << " "
              << Qi.z() << " "
              << Qi.w() << endl;
        old_time = time_stamp;

    }
    foutC.close();

    std::cout << "finish\n";
    std::cout << "finish\n";
#endif
    return 0;
}
#else
int main(int argc, char** argv) {


    string config_file = argv[1];
    ROS_PATH = argv[2];
    RESULT_PATH = argv[3];
    printf("config_file: %s\n", argv[1]);
    readParameters(config_file);
#if ENABLE_STEREO_MAP||DEBUG
    ros::init(argc, argv, "flower_focus");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    registerPub(n);
#endif
    LEAK_NUM = SWF_SIZE_IN - 1;

#if DEBUG
    time_t t = time(nullptr);
    struct tm* now = localtime(&t);
    std::stringstream timeStr;
    timeStr << now->tm_year + 1900 << "-";
    timeStr << now->tm_mon + 1 << "-";
    timeStr << now->tm_mday << " ";
    timeStr << now->tm_hour << ":";
    timeStr << now->tm_min << ":";
    timeStr << now->tm_sec << ".log";
    std::string LOG_PATH = "log/" + RESULT_PATH + "--" + timeStr.str();
    LOG_OUT = std::ofstream(LOG_PATH, std::ios::out);
    LOG_OUT.precision(10);
#endif


    swf_optimization = new SWFOptimization();
    swf_optimization->SetParameter();
    printf("waiting for image and imu...");
    rosbag::Bag bag;
    bag.open(ROS_PATH, rosbag::bagmode::Read);
    rosbag::View view(bag);
#if ENABLE_STEREO_MAP
    stereo_init(config_file, n);
#endif
    signal( SIGINT, sig_handler );
    for ( rosbag::View::iterator it = view.begin(); it != view.end(); ++it) {
        auto m = *it;
        if (m.getTopic() == IMU_TOPIC)
            imu_callback(m.instantiate<sensor_msgs::Imu>());
        else if (m.getTopic() == IMAGE0_TOPIC)
            img0_callback(m.instantiate<sensor_msgs::Image>());
#if ENABLE_STEREO_MAP||USE_STEREO
        else if (m.getTopic() == IMAGE1_TOPIC)
            img1_callback(m.instantiate<sensor_msgs::Image>());
#endif
        else if (m.getTopic() == RTK_TOPIC)
            mGNSS_callback(m.instantiate<std_msgs::ByteMultiArray>());
        else if (m.getTopic() == MAG_TOPIC)
            mag_callback(m.instantiate<sensor_msgs::MagneticField>());
        processoneimage();
#if !REAL_TIME
        process_stereo();
#endif
    }
    while (swf_optimization->feature_buf.size() > 10);


    //As same as ORB-SLAM3, we generate the final resutls for evaluation.
    std::vector<PosInfo> camera_poses = swf_optimization->RetriveAllPose();
    ofstream foutC(RESULT_PATH, ios::out);
    foutC.setf(ios::fixed, ios::floatfield);

    double scale_factor = swf_optimization->scale_factor;
    static double old_time = 0;
    for (int i = 0; i < (int)camera_poses.size(); i++) {
        Eigen::Matrix3d Ri = R_WI_WC * camera_poses[i].R * RIC[0].transpose();
        Eigen::Vector3d ti = R_WI_WC * camera_poses[i].t * scale_factor - Ri * TIC[0] - Ri * Pgb + P_WI_WC;
        Ri = Rwgw.transpose() * Ri;
        ti = Rwgw.transpose() * ti;

        Eigen::Quaterniond Qi = Eigen::Quaterniond(Ri);
        double time_stamp = camera_poses[i].time_stamp;
        if (old_time)assert(time_stamp > old_time);

        foutC.precision(9);
        foutC << time_stamp << " ";
        foutC.precision(5);
        foutC << ti.x() << " "
              << ti.y() << " "
              << ti.z() << " "
              << Qi.x() << " "
              << Qi.y() << " "
              << Qi.z() << " "
              << Qi.w() << endl;
        old_time = time_stamp;

    }
    foutC.close();

    std::cout << "finish\n";
    std::cout << "finish\n";
    return 0;
}
#endif



