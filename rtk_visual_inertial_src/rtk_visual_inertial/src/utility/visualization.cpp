

#include "visualization.h"
#include "../parameter/parameters.h"
#include <fstream>
#include <std_msgs/Header.h>
#include <sensor_msgs/PointCloud.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>
#include "camera_pose_visualization.h"
#include "common_function.h"
using namespace ros;
using namespace Eigen;
ros::Publisher pub_odometry, pub_latest_odometry;
ros::Publisher pub_path;
ros::Publisher pub_point_cloud, pub_margin_cloud;
ros::Publisher pub_key_poses;
ros::Publisher pub_camera_pose;
ros::Publisher pub_camera_pose_right;
ros::Publisher pub_rectify_pose_left;
ros::Publisher pub_rectify_pose_right;
ros::Publisher pub_camera_pose_visual;
nav_msgs::Path path;

ros::Publisher pub_keyframe_pose;
ros::Publisher pub_keyframe_point;
ros::Publisher pub_extrinsic;

camera_pose_visualization cameraposevisual(1, 0, 0, 1);
static double sum_of_path = 0;
static Vector3d last_path(0.0, 0.0, 0.0);
#define scalefactor 1
size_t pub_counter = 0;
Vector3d        Ps0;
void registerPub(ros::NodeHandle& n) {
    pub_latest_odometry = n.advertise<nav_msgs::Odometry>("imu_propagate", 1000);
    pub_path = n.advertise<nav_msgs::Path>("path", 1000);
    pub_odometry = n.advertise<nav_msgs::Odometry>("odometry", 1000);
    pub_point_cloud = n.advertise<sensor_msgs::PointCloud>("point_cloud", 1000);
    pub_margin_cloud = n.advertise<sensor_msgs::PointCloud>("margin_cloud", 1000);
    pub_camera_pose = n.advertise<nav_msgs::Odometry>("camera_pose", 1000);
    pub_camera_pose_right = n.advertise<nav_msgs::Odometry>("camera_pose_right", 1000);
    pub_rectify_pose_left = n.advertise<geometry_msgs::PoseStamped>("rectify_pose_left", 1000);
    pub_rectify_pose_right = n.advertise<geometry_msgs::PoseStamped>("rectify_pose_right", 1000);
    pub_camera_pose_visual = n.advertise<visualization_msgs::MarkerArray>("camera_pose_visual", 1000);
    pub_keyframe_pose = n.advertise<nav_msgs::Odometry>("keyframe_pose", 1000);
    pub_keyframe_point = n.advertise<sensor_msgs::PointCloud>("keyframe_point", 1000);
    pub_extrinsic = n.advertise<nav_msgs::Odometry>("extrinsic", 1000);

    // cameraposevisual.setScale(10);
    // cameraposevisual.setLineWidth(10);
}

Eigen::Vector3d ecef2geo_google_map(const Eigen::Vector3d& xyz) {

#define EARTH_ECCE_2            6.69437999014e-3    // WGS 84 (Earth eccentricity)^2 (m^2)
#define EARTH_MEAN_RADIUS       6371009             // Mean R of ellipsoid(m) IU Gedosey& Geophysics
#define EARTH_SEMI_MAJOR        6378137             // WGS 84 Earth semi-major axis (m)
#define R2D                     (180.0/M_PI)        // radius to degree


    Eigen::Vector3d lla = Eigen::Vector3d::Zero();
    if (xyz.x() == 0 && xyz.y() == 0) {
        LOG(ERROR) << "LLA coordinate is not defined if x = 0 and y = 0";
        return lla;
    }

    double e2 = EARTH_ECCE_2;
    double a = EARTH_SEMI_MAJOR;
    double a2 = a * a;
    double b2 = a2 * (1 - e2);
    double b = sqrt(b2);
    double ep2 = (a2 - b2) / b2;
    double p = xyz.head<2>().norm();

    // two sides and hypotenuse of right angle triangle with one angle = theta:
    double s1 = xyz.z() * a;
    double s2 = p * b;
    double h = sqrt(s1 * s1 + s2 * s2);
    double sin_theta = s1 / h;
    double cos_theta = s2 / h;

    // two sides and hypotenuse of right angle triangle with one angle = lat:
    s1 = xyz.z() + ep2 * b * pow(sin_theta, 3);
    s2 = p - a * e2 * pow(cos_theta, 3);
    h = sqrt(s1 * s1 + s2 * s2);
    double tan_lat = s1 / s2;
    double sin_lat = s1 / h;
    double cos_lat = s2 / h;
    double lat = atan(tan_lat);
    double lat_deg = lat * R2D;

    double N = a2 * pow((a2 * cos_lat * cos_lat + b2 * sin_lat * sin_lat), -0.5);
    double altM = p / cos_lat - N;

    double lon = atan2(xyz.y(), xyz.x());
    double lon_deg = lon * R2D;
    lla << lat_deg, lon_deg, altM;
    return lla;
}


void printStatistics(const SWFOptimization& swf_optimization, double t) {

    if (swf_optimization.solver_flag != SWFOptimization::SolverFlag::NonLinear)
        return;

    int index = swf_optimization.rover_count + swf_optimization.image_count - 1;
    printf("time:%.3f\r\n", swf_optimization.headers[index]);
    // assert(swf_optimization.headers[index]<2024.681);
    std::cout << "pos: " << (swf_optimization.InitRwgw.transpose()*swf_optimization.Ps[index]).transpose() << " " << (swf_optimization.Ps[index]).transpose() << "\t";
    std::cout << "vel: " << (swf_optimization.InitRwgw.transpose()*swf_optimization.Vs[index]).transpose() << "\t";
    std::cout << "orientation: " << Utility::R2ypr(swf_optimization.InitRwgw.transpose()*swf_optimization.Rs[index]).transpose() << "\t";
    std::cout << "gyro bias: " << swf_optimization.Bgs[index].transpose() << "\t";
    std::cout << "acc bias: " << swf_optimization.Bas[index].transpose() << "\t";
    std::cout << "rover_count:" << swf_optimization.rover_count << "\t image_count:" << swf_optimization.image_count << "\t";

    static int count = 0;
    static double count2 = 0;
    count += swf_optimization.rover_count;
    count2 += 1;
    std::cout << "rover ratio:" << count / count2 << "\t";


    for (int i = 0; i < swf_optimization.image_count + 1; i++) {
        if (swf_optimization.imu_gnss_factor[i])
            std::cout << swf_optimization.imu_gnss_factor[i]->update_round * 1.0 / swf_optimization.imu_gnss_factor[i]->evaluate_round << ",";
    }
    std::cout << std::endl;



    if (ESTIMATE_EXTRINSIC) {
        cv::FileStorage fs(EX_CALIB_RESULT_PATH, cv::FileStorage::WRITE);
        for (int i = 0; i < NUM_OF_CAM; i++) {
            Eigen::Matrix4d eigen_T = Eigen::Matrix4d::Identity();
            eigen_T.block<3, 3>(0, 0) = swf_optimization.ric[i];
            eigen_T.block<3, 1>(0, 3) = swf_optimization.tic[i];
            cv::Mat cv_T;
            cv::eigen2cv(eigen_T, cv_T);
            if (i == 0)
                fs << "body_T_cam0" << cv_T ;
            else
                fs << "body_T_cam1" << cv_T ;
        }
        fs.release();
    }



    static double sum_of_time = 0;
    static int sum_of_calculation = 0;
    sum_of_time += t;
    sum_of_calculation++;

    sum_of_path += (swf_optimization.Ps[swf_optimization.rover_count + swf_optimization.image_count - 1] - last_path).norm();
    last_path = swf_optimization.Ps[swf_optimization.rover_count + swf_optimization.image_count - 1];

}
void resetpot(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {
    Ps0 = swf_optimization.Ps[0];
}
void pubOdometry(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {
    Vector3d        Ps[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];


    for (int i = 0; i < swf_optimization.rover_count + swf_optimization.image_count; i++) {
        Ps[i] = swf_optimization.Ps[i] - Ps0;
    }
    if (swf_optimization.solver_flag == SWFOptimization::SolverFlag::NonLinear) {
        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(swf_optimization.InitRwgw.transpose() * swf_optimization.Rs[swf_optimization.rover_count + swf_optimization.image_count - 1]);
        Vector3d tmpP = swf_optimization.InitRwgw.transpose() * Ps[swf_optimization.rover_count + swf_optimization.image_count - 1] * scalefactor;

        Vector3d tmpv = swf_optimization.InitRwgw.transpose() * swf_optimization.Vs[swf_optimization.rover_count + swf_optimization.image_count - 1] * scalefactor;

        odometry.pose.pose.position.x = tmpP.x();
        odometry.pose.pose.position.y = tmpP.y();
        odometry.pose.pose.position.z = tmpP.z();
        odometry.pose.pose.orientation.x = tmp_Q.x();
        odometry.pose.pose.orientation.y = tmp_Q.y();
        odometry.pose.pose.orientation.z = tmp_Q.z();
        odometry.pose.pose.orientation.w = tmp_Q.w();
        odometry.twist.twist.linear.x = tmpv.x();
        odometry.twist.twist.linear.y = tmpv.y();
        odometry.twist.twist.linear.z = tmpv.z();
        pub_odometry.publish(odometry);

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose = odometry.pose.pose;
        path.header = header;
        path.header.frame_id = "world";
        path.poses.push_back(pose_stamped);
        pub_path.publish(path);
        if (swf_optimization.rover_count + swf_optimization.image_count - 1 < 0)return;
        // write result to file

        {
            nav_msgs::Odometry odometry;
            odometry.header = header;
            odometry.header.frame_id = "world";
            odometry.pose.pose.position.x = swf_optimization.tic[0].x();
            odometry.pose.pose.position.y = swf_optimization.tic[0].y();
            odometry.pose.pose.position.z = swf_optimization.tic[0].z();
            Quaterniond tmp_q{swf_optimization.ric[0]};
            odometry.pose.pose.orientation.x = tmp_q.x();
            odometry.pose.pose.orientation.y = tmp_q.y();
            odometry.pose.pose.orientation.z = tmp_q.z();
            odometry.pose.pose.orientation.w = tmp_q.w();
            pub_extrinsic.publish(odometry);
        }

    }
}

void save_result(const SWFOptimization& swf_optimization) {


    static bool f;
    if (!f) {
        f = true;
        std::ofstream fout(RESULT_PATH, std::ios::out);
        fout << "time,px,py,pz,vx,vy,vz,yaw,pitch,roll,ax,ay,az,bax,bay,baz,bgx,bgy,bgz,rtk_fix,mag_yaw,havegps,trajectory_sum,lat,lon,pbgx,pbgy,pbgz";
        fout << std::endl;
        fout.close();
    }


    ofstream foutC(RESULT_PATH, ios::app);
    foutC.setf(ios::fixed, ios::floatfield);
    foutC.precision(0);
    foutC << swf_optimization.headers[swf_optimization.rover_count + swf_optimization.image_count - 1] * 1e9 << ",";
    foutC.precision(10);
    Eigen::Vector3d xyz = swf_optimization.InitRwgw.transpose() * (
                              swf_optimization.base_pos + swf_optimization.Ps[swf_optimization.rover_count + swf_optimization.image_count - 1]
                          );
    Eigen::Vector3d vxyz = swf_optimization.InitRwgw.transpose() * swf_optimization.Vs[swf_optimization.rover_count + swf_optimization.image_count - 1];
    Eigen::Vector3d ypr = Utility::R2ypr(swf_optimization.InitRwgw.transpose() * swf_optimization.Rs[swf_optimization.rover_count + swf_optimization.image_count - 1]);
    Eigen::Vector3d acc_w = swf_optimization.InitRwgw.transpose() * (swf_optimization.Rs[swf_optimization.rover_count + swf_optimization.image_count - 1] * swf_optimization.acc_0) - G;

    Eigen::Vector3d ba = swf_optimization.Bas[swf_optimization.rover_count + swf_optimization.image_count - 1];
    Eigen::Vector3d bg = swf_optimization.Bgs[swf_optimization.rover_count + swf_optimization.image_count - 1];
    Eigen::Vector3d Pbgw = swf_optimization.InitRwgw.transpose() * (swf_optimization.Rs[swf_optimization.rover_count + swf_optimization.image_count - 1] * Pbg);
    int havegps = (swf_optimization.frame_types[swf_optimization.rover_count + swf_optimization.image_count - 1] == swf_optimization.GnssFrame ||
                   swf_optimization.frame_types[swf_optimization.rover_count + swf_optimization.image_count - 2] == swf_optimization.GnssFrame ||
                   swf_optimization.frame_types[swf_optimization.rover_count + swf_optimization.image_count - 3] == swf_optimization.GnssFrame) ? 1 : 0;

    int rtk_fix = fabs(swf_optimization.last_fix_time - swf_optimization.headers[swf_optimization.rover_count + swf_optimization.image_count - 1]) < 0.2 ? 1 : 0;
    double mag_yaw = swf_optimization.mag_yaw;
    double distance = 0;

    Eigen::Vector3d xyz2 = swf_optimization.base_pos + swf_optimization.Ps[swf_optimization.rover_count + swf_optimization.image_count - 1];

    Eigen::Vector3d lla = ecef2geo_google_map( xyz2);

    foutC << xyz.x() << "," << xyz.y() << "," << xyz.z() << ","
          << vxyz.x() << "," << vxyz.y() << "," << vxyz.z() << ","
          << ypr.x() << "," << ypr.y() << "," << ypr.z() << ","
          << acc_w.x() << "," << acc_w.y() << "," << acc_w.z() << ","
          << ba.x() << "," << ba.y() << "," << ba.z() << ","
          << bg.x() << "," << bg.y() << "," << bg.z() << ","
          << rtk_fix << ","
          << mag_yaw << ","
          << havegps << ","
          << distance << ","
          << lla.x() << "," << lla.y() << ","
          << Pbgw.x() << "," << Pbgw.y() << "," << Pbgw.z();

    foutC << std::endl;
    foutC.close();
}


void pubCameraPose(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {

    Vector3d        Ps[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];


    for (int i = 0; i < swf_optimization.rover_count + swf_optimization.image_count; i++) {
        Ps[i] = swf_optimization.Ps[i] - Ps0;
    }

    int idx2 = swf_optimization.rover_count + swf_optimization.image_count - 1;
    if (swf_optimization.solver_flag == SWFOptimization::SolverFlag::NonLinear) {
        int i = idx2;
        Vector3d P = swf_optimization.InitRwgw.transpose() * (Ps[i] + swf_optimization.Rs[i] * swf_optimization.tic[0]) * scalefactor;
        Quaterniond R = Quaterniond(swf_optimization.InitRwgw.transpose() * (swf_optimization.Rs[i] * swf_optimization.ric[0]));

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();

        if (USE_STEREO) {
            Vector3d P_r = swf_optimization.InitRwgw.transpose() * (Ps[i] + swf_optimization.Rs[i] * swf_optimization.tic[1]) * scalefactor;
            Quaterniond R_r = Quaterniond(swf_optimization.InitRwgw.transpose() * (swf_optimization.Rs[i] * swf_optimization.ric[1]));

            nav_msgs::Odometry odometry_r;
            odometry_r.header = header;
            odometry_r.header.frame_id = "world";
            odometry_r.pose.pose.position.x = P_r.x();
            odometry_r.pose.pose.position.y = P_r.y();
            odometry_r.pose.pose.position.z = P_r.z();
            odometry_r.pose.pose.orientation.x = R_r.x();
            odometry_r.pose.pose.orientation.y = R_r.y();
            odometry_r.pose.pose.orientation.z = R_r.z();
            odometry_r.pose.pose.orientation.w = R_r.w();
            pub_camera_pose_right.publish(odometry_r);
        }

        pub_camera_pose.publish(odometry);

        cameraposevisual.reset();
        cameraposevisual.add_pose(P, R);
        if (USE_STEREO) {
            Vector3d P = swf_optimization.InitRwgw.transpose() * (Ps[i] + swf_optimization.Rs[i] * swf_optimization.tic[1]) * scalefactor;
            Quaterniond R = Quaterniond(swf_optimization.InitRwgw.transpose() * (swf_optimization.Rs[i] * swf_optimization.ric[1]));
            cameraposevisual.add_pose(P, R);
        }
        cameraposevisual.publish_by(pub_camera_pose_visual, odometry.header);
    }
}


void pubPointCloud(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {
    Vector3d        Ps[(FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1)];


    for (int i = 0; i < swf_optimization.rover_count + swf_optimization.image_count; i++) {
        Ps[i] = swf_optimization.Ps[i] - Ps0;
    }

    sensor_msgs::PointCloud point_cloud, loop_point_cloud;
    point_cloud.header = header;
    loop_point_cloud.header = header;

    for (auto& it_per_id : swf_optimization.f_manager.feature) {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < swf_optimization.rover_count + swf_optimization.image_count - 1 - 2))
            continue;
        if (it_per_id.start_frame > (swf_optimization.rover_count + swf_optimization.image_count - 1) * 3.0 / 4.0 || it_per_id.solve_flag != 1)
            continue;

        Vector3d w_pts_i = swf_optimization.InitRwgw.transpose() * (it_per_id.ptsInWorld - Ps0);

        geometry_msgs::Point32 p;
        p.x = w_pts_i(0);
        p.y = w_pts_i(1);
        p.z = w_pts_i(2);
        point_cloud.points.push_back(p);
    }
    pub_point_cloud.publish(point_cloud);


    // pub margined potin
    sensor_msgs::PointCloud margin_cloud;
    margin_cloud.header = header;

    for (auto& it_per_id : swf_optimization.f_manager.feature) {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < swf_optimization.rover_count + swf_optimization.image_count - 1 - 2))
            continue;

        if (it_per_id.start_frame == 0 && it_per_id.feature_per_frame.size() <= 2
                && it_per_id.solve_flag == 1 ) {
            Vector3d w_pts_i = swf_optimization.InitRwgw.transpose() * (it_per_id.ptsInWorld - Ps0);

            geometry_msgs::Point32 p;
            p.x = w_pts_i(0);
            p.y = w_pts_i(1);
            p.z = w_pts_i(2);
            margin_cloud.points.push_back(p);
        }
    }
    pub_margin_cloud.publish(margin_cloud);
}



void pubKeyframe(const SWFOptimization& swf_optimization) {
    // pub camera pose, 2D-3D points of keyframe
    if (swf_optimization.solver_flag == SWFOptimization::SolverFlag::NonLinear && swf_optimization.marg_flag == 0) {
        int i = swf_optimization.i2f[FEATURE_WINDOW_SIZE - 2];
        Vector3d P = swf_optimization.Ps[i] - swf_optimization.Rs[i] * Pbg;
        Quaterniond R = Quaterniond(swf_optimization.Rs[i]);

        nav_msgs::Odometry odometry;
        odometry.header.stamp = ros::Time(swf_optimization.headers[i]);
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();
        //printf("time: %f t: %f %f %f r: %f %f %f %f\n", odometry.header.stamp.toSec(), P.x(), P.y(), P.z(), R.w(), R.x(), R.y(), R.z());

        pub_keyframe_pose.publish(odometry);


        sensor_msgs::PointCloud point_cloud;
        point_cloud.header.stamp = ros::Time(swf_optimization.headers[i]);
        point_cloud.header.frame_id = "world";
        for (auto& it_per_id : swf_optimization.f_manager.feature) {
            int frame_size = it_per_id.feature_per_frame.size();
            if (it_per_id.start_frame < FEATURE_WINDOW_SIZE - 2
                    && it_per_id.start_frame + frame_size - 1 >= FEATURE_WINDOW_SIZE - 2
                    && it_per_id.solve_flag == 1) {

                Vector3d w_pts_i = it_per_id.ptsInWorld;
                geometry_msgs::Point32 p;
                p.x = w_pts_i(0);
                p.y = w_pts_i(1);
                p.z = w_pts_i(2);
                point_cloud.points.push_back(p);

                int imu_j = FEATURE_WINDOW_SIZE - 2 - it_per_id.start_frame;
                sensor_msgs::ChannelFloat32 p_2d;
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.y());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.y());
                p_2d.values.push_back(it_per_id.feature_id);
                point_cloud.channels.push_back(p_2d);
            }

        }
        pub_keyframe_point.publish(point_cloud);
    }
}