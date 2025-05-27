

#pragma once

#include <ros/ros.h>

#include <eigen3/Eigen/Dense>

#include <sensor_msgs/PointCloud.h>
#include "../swf/swf.h"
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>



void registerPub(ros::NodeHandle& n);

void printStatistics(const SWFOptimization& swf_optimization, double t);

void pubOdometry(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubCameraPose(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubPointCloud(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void resetpot(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void save_result(const SWFOptimization& swf_optimization,double time);