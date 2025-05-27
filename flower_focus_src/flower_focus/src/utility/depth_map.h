#pragma once
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


void stereo_img0_callback(const sensor_msgs::ImageConstPtr& img_msg);
void stereo_img1_callback(const sensor_msgs::ImageConstPtr& img_msg);
void stereo_init(string config_file,ros::NodeHandle& n);
void process_stereo() ;