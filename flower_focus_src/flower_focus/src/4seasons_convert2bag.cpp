#include <csignal>
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
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <std_msgs/ByteMultiArray.h>
#include <sensor_msgs/MagneticField.h>
#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>

rosbag::Bag bag_write;

void imu_pub(double* acc, double* gyr, std::string imu_frame_id, double timestamp) {

    sensor_msgs::Imu msg;
    msg.header.stamp =  ros::Time(timestamp);
    msg.header.frame_id = imu_frame_id;

    msg.linear_acceleration.x = acc[0];
    msg.linear_acceleration.y = acc[1];
    msg.linear_acceleration.z = acc[2];

    msg.angular_velocity.x = gyr[0];
    msg.angular_velocity.y = gyr[1];
    msg.angular_velocity.z = gyr[2];

    bag_write.write(IMU_TOPIC, msg.header.stamp, msg);
}


void image_pub(cv::Mat img, std::string img_frame_id, double timestamp) {


    sensor_msgs::ImagePtr ptr = boost::make_shared<sensor_msgs::Image>();
    sensor_msgs::Image& imgMessage = *ptr;

    imgMessage.header.stamp = ros::Time(timestamp);
    imgMessage.header.frame_id = img_frame_id;
    imgMessage.height = 400;
    imgMessage.width = 800;


    imgMessage.step = imgMessage.width ;

    imgMessage.data.resize(imgMessage.width * imgMessage.height);
    imgMessage.encoding = sensor_msgs::image_encodings::MONO8;
    memcpy((char*)(&imgMessage.data[0]), img.data, imgMessage.height * imgMessage.width);

    bag_write.write(img_frame_id, ros::Time(timestamp), ptr);

}

void sig_handler( int sig ) {
    if ( sig == SIGINT) {
        exit( 3 );
        LOG_OUT << "abort_\r\n";
    }
}



int main(int argc, char** argv) {

    string config_file = argv[1];
    readParameters(config_file);
    std::string ROS_PATH = argv[2];
    std::string imu_file = argv[3];
    std::string picture_path = argv[4];


    bag_write.open(ROS_PATH, rosbag::bagmode::Write); //打开一个bag文件

//imu
    {
        static int count = 0;
        std::string line;
        std::ifstream in = std::ifstream(imu_file);
        std::cout << imu_file << std::endl;
        if (in.fail())
            cout << "File not found" << endl;
        getline(in, line);
        while (getline(in, line)  && in.good() ) {

            vector<double>b;
            std::stringstream sin(line);
            string field;
            string::size_type size; //13

            while (getline(sin, field, ' ')) {
                double d = stod(field, &size);  //1234.56789
                b.push_back(d);
            }

            double time = b[0] / 1e9;
            double acc[3], gyr[3];
            for (int i = 0; i < 3; i++)gyr[i] = b[1 + i];
            for (int i = 0; i < 3; i++)acc[i] = b[4 + i];
            count++;
            if (count % 2000 == 0)printf("%.9f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\r\n", time, gyr[0], gyr[1], gyr[2], acc[0], acc[1], acc[2]);

            imu_pub(acc, gyr, IMU_TOPIC, time);
        }
    }



//camera
    {
        static int count = 0;
        std::string line;
        std::ifstream in = std::ifstream(picture_path + "/times.txt");
        if (in.fail())cout << "File not found" << endl;
        getline(in, line);
        while (getline(in, line)  && in.good() ) {

            vector<double>b;
            std::stringstream sin(line);
            string field;
            string::size_type size; //13

            while (getline(sin, field, ' ')) {
                double d = stod(field, &size);  //1234.56789
                b.push_back(d);
            }

            double time = b[1] ;
            double exposure_time = b[2];
            std::string name = std::to_string((long)b[0]) + ".png";
            cv::Mat image1 = cv::imread(picture_path + "/cam0/" + name,cv::IMREAD_GRAYSCALE);
            // cv::cvtColor(image1, image1, cv::COLOR_BGR2GRAY);
            cv::imshow("left camera", image1);
            cv::waitKey(1);
            image_pub(image1.clone(), IMAGE0_TOPIC, time);
            count++;
            if (count % 30 == 0)printf("%.9f\t%.5f\r\n", time,exposure_time);
        }
    }

    bag_write.close();



    return 0;
}




