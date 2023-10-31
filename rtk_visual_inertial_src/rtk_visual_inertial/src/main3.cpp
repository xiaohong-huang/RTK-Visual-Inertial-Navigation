

#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "swf/swf.h"
#include "parameter/parameters.h"
#include "utility/visualization.h"
#include "utility/utility.h"
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include "common_function.h"
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


SWFOptimization* swf_optimization;

queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::ImageConstPtr> img0_buf;
queue<sensor_msgs::ImageConstPtr> img1_buf;
std::mutex m_buf;

double start_timestamp = 0;
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
    if (USE_STEREO) {
        if (!img0_buf.empty() && !img1_buf.empty()) {
            double time0 = img0_buf.front()->header.stamp.toSec();
            double time1 = img1_buf.front()->header.stamp.toSec();
            if (fabs(time0 - time1) < 1e-3)time0 = time1;
            // std::cout<<time0<<","<<time1<<std::endl;
            if (time0 < time1) {
                img0_buf.pop();
                printf("throw img0\n");
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
    } else {
        if (!img0_buf.empty() ) {
            time = img0_buf.front()->header.stamp.toSec();
            image0 = getImageFromMsg(img0_buf.front());
            img0_buf.pop();

        }
    }

    if (!image0.empty())
        swf_optimization->InputImage(time, image0, image1);
}

void img0_callback(const sensor_msgs::ImageConstPtr& img_msg) {
    if (img_msg->header.stamp.toSec() < start_timestamp)return;
    img0_buf.push(img_msg);
    processoneimage();
}

void img1_callback(const sensor_msgs::ImageConstPtr& img_msg) {
    if (img_msg->header.stamp.toSec() < start_timestamp)return;
    img1_buf.push(img_msg);
    processoneimage();
}




void imu_callback(const sensor_msgs::ImuConstPtr& imu_msg) {
    if (imu_msg->header.stamp.toSec() < start_timestamp - 1)return;
    static double oldt;
    double t = imu_msg->header.stamp.toSec();

    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Vector3d acc(dx, dy, dz);
    Vector3d gyr(rx, ry, rz);
    acc = IMUMatrix * acc;
    gyr = IMUMatrix * gyr;
    swf_optimization->InputIMU(t, acc, gyr);
    oldt = t;

    return;
}

void mag_callback(const sensor_msgs::MagneticFieldPtr& mag_msg) {
    if (mag_msg->header.stamp.toSec() < start_timestamp)return;
    double t = mag_msg->header.stamp.toSec();
    double mx = mag_msg->magnetic_field.x * 1e6;
    double my = mag_msg->magnetic_field.y * 1e6;
    double mz = mag_msg->magnetic_field.z * 1e6;
    Vector3d mag(mx, my, mz);
    mag = MagMatrix * mag + MagVector;
    swf_optimization->InputMag(t, mag);

    return;
}






void mGNSS_callback(const std_msgs::ByteMultiArrayPtr& gnss_msg) {


    mea_t gnss_obs;
    char* rover_c = (char*)(&gnss_obs);
    memset(rover_c, 0, sizeof(mea_t));
    for (int i = 0; i < (int)(gnss_msg->data.size()); i++) {
        rover_c[i] = gnss_msg->data[i];
    }
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

void feature_callback(const sensor_msgs::PointCloudConstPtr& feature_msg) {
    if (feature_msg->header.stamp.toSec() < start_timestamp)return;
    static int inputImageCnt;

    if (inputImageCnt < AVERAGE_IMAGE) {
        inputImageCnt++;
        return;
    }
    if (USE_IMAGE) {

        double timestamp = feature_msg->header.stamp.toSec();
        std::map<int, std::vector<std::pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
        for (unsigned int i = 0; i < feature_msg->points.size(); i++) {
            int feature_id = feature_msg->channels[0].values[i];
            int camera_id = feature_msg->channels[1].values[i];
            double x = feature_msg->points[i].x;
            double y = feature_msg->points[i].y;
            double z = feature_msg->points[i].z;
            double p_u = feature_msg->channels[2].values[i];
            double p_v = feature_msg->channels[3].values[i];
            double velocity_x = feature_msg->channels[4].values[i];
            double velocity_y = feature_msg->channels[5].values[i];
            assert(z == 1);
            Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
            xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
            featureFrame[feature_id].emplace_back(camera_id,  xyz_uv_velocity);
        }
        swf_optimization->feature_buf.push(make_pair(timestamp, featureFrame));

        swf_optimization->MeasurementProcess();
    }
}

int main(int argc, char** argv) {
    ros::init(argc, argv, "rtk_visual_inertial");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);



    string config_file = argv[1];
    ROS_PATH = argv[2];
    RESULT_PATH = argv[3];

    printf("config_file: %s\n", argv[1]);


    readParameters(config_file);
    swf_optimization = new SWFOptimization();
    swf_optimization->SetParameter();


    printf("waiting for image and imu...");

    registerPub(n);

    rosbag::Bag bag;
    bag.open(ROS_PATH, rosbag::bagmode::Read);
    rosbag::View view(bag);

    for ( rosbag::View::iterator it = view.begin(); it != view.end(); ++it) {
        auto m = *it;
        if (m.getTopic() == IMU_TOPIC) {
            imu_callback(m.instantiate<sensor_msgs::Imu>());
        } else if (m.getTopic() == IMAGE0_TOPIC && !USE_FEATURE) {
            img0_callback(m.instantiate<sensor_msgs::Image>());
        } else if (m.getTopic() == IMAGE1_TOPIC && !USE_FEATURE) {
            img1_callback(m.instantiate<sensor_msgs::Image>());
        } else if (m.getTopic() == RTK_TOPIC) {
            mGNSS_callback(m.instantiate<std_msgs::ByteMultiArray>());
        } else if (m.getTopic() == MAG_TOPIC) {
            mag_callback(m.instantiate<sensor_msgs::MagneticField>());
        } else if (m.getTopic() == FEATURE_TOPIC && USE_FEATURE) {
            feature_callback(m.instantiate<sensor_msgs::PointCloud>());
        }
        processoneimage();
    }

    return 0;
}




