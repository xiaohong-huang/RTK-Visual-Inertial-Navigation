

#pragma once

#include <ros/ros.h>

#include <eigen3/Eigen/Dense>

#include "../swf/swf.h"



void registerPub(ros::NodeHandle& n);

void printStatistics(const SWFOptimization& swf_optimization, double t);

void pubOdometry(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubInitialGuess(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubKeyPoses(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubCameraPose(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubPointCloud(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubTF(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void pubKeyframe(const SWFOptimization& swf_optimization);

void pubRelocalization(const SWFOptimization& swf_optimization);

void pubCar(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void resetpot(const SWFOptimization& swf_optimization, const std_msgs::Header& header);

void save_result(const SWFOptimization& swf_optimization);