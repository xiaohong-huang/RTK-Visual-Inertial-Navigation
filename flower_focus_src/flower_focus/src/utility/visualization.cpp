
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

#include "visualization.h"
#include "../parameter/parameters.h"
#include <fstream>
#include <std_msgs/Header.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>
#include "camera_pose_visualization.h"



using namespace ros;
using namespace Eigen;

nav_msgs::Path path;

ros::Publisher pub_odometry, pub_latest_odometry;
ros::Publisher pub_path;
ros::Publisher pub_point_cloud_short, pub_margin_cloud, pub_point_cloud_long;
ros::Publisher pub_key_poses;
ros::Publisher pub_camera_pose;
ros::Publisher pub_camera_pose_right;
ros::Publisher pub_rectify_pose_left;
ros::Publisher pub_rectify_pose_right;
ros::Publisher pub_camera_pose_visual;
ros::Publisher pub_keyframe_pose;
ros::Publisher pub_keyframe_point;
ros::Publisher pub_extrinsic;

camera_pose_visualization cameraposevisual(1, 0, 0, 1);
Vector3d        Ps0;

#define PUB_INDEX swf_optimization.image_count - 1


void registerPub(ros::NodeHandle& n) {
    pub_latest_odometry = n.advertise<nav_msgs::Odometry>("imu_propagate", 1000);
    pub_path = n.advertise<nav_msgs::Path>("path", 1000);
    pub_odometry = n.advertise<nav_msgs::Odometry>("odometry", 1000);
    pub_point_cloud_short = n.advertise<sensor_msgs::PointCloud>("point_cloud", 1000);

    pub_margin_cloud = n.advertise<sensor_msgs::PointCloud>("margin_cloud", 1000);
    pub_camera_pose = n.advertise<nav_msgs::Odometry>("camera_pose", 1000);
    pub_camera_pose_right = n.advertise<nav_msgs::Odometry>("camera_pose_right", 1000);
    pub_rectify_pose_left = n.advertise<geometry_msgs::PoseStamped>("rectify_pose_left", 1000);
    pub_rectify_pose_right = n.advertise<geometry_msgs::PoseStamped>("rectify_pose_right", 1000);
    pub_camera_pose_visual = n.advertise<visualization_msgs::MarkerArray>("camera_pose_visual", 1000);
    pub_keyframe_pose = n.advertise<nav_msgs::Odometry>("keyframe_pose", 1000);
    pub_keyframe_point = n.advertise<sensor_msgs::PointCloud>("keyframe_point", 1000);
    pub_extrinsic = n.advertise<nav_msgs::Odometry>("extrinsic", 1000);
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

void save_result(const SWFOptimization& swf_optimization, double time) {
    if (!enable_output)return;

    static bool f;
    if (!f) {
        f = true;
        std::ofstream fout(RESULT_PATH + "_real", std::ios::out);
        fout.close();
    }


    ofstream foutC(RESULT_PATH + "_real", ios::app);
    foutC.setf(ios::fixed, ios::floatfield);
    foutC.precision(5);
    foutC << time << " ";

// Rwgw.transpose()*R_WI_WC.transpose() * (swf_optimization.Ps[i] - P_WI_WC)
    // if(Pbg.norm())Pgb=-Pbg;
    Eigen::Vector3d xyz = InitRwgw.transpose() * R_WI_WC.transpose() * (
                              swf_optimization.Ps[ swf_optimization.image_count - 1] - swf_optimization.Rs[ swf_optimization.image_count - 1] * Pgb - P_WI_WC
                          );

    foutC << xyz.x() << " " << xyz.y() << " " << xyz.z() << " "
          << 0 << " " << 0 << " " << 0 << " "
          << 0 ;

    foutC << std::endl;
    foutC.close();


}


void saveMultipleExtrinsicsToYaml(const std::vector<Eigen::Vector3d>& TICs, const std::vector<Eigen::Quaterniond>& QICs, Eigen::Vector3d& Pbg, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "无法创建文件: " << filename << std::endl;
        return;
    }

    if (TICs.size() != QICs.size()) {
        std::cerr << "TICs 和 QICs 的大小不一致！" << std::endl;
        return;
    }

    for (size_t i = 0; i < TICs.size(); ++i) {
        Eigen::Vector3d TIC = TICs[i];
        Eigen::Quaterniond QIC = QICs[i];

        // 创建变换矩阵
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        T.block<3, 3>(0, 0) = QIC.toRotationMatrix();  // 旋转部分
        T.block<3, 1>(0, 3) = TIC;                    // 平移部分

        // 创建 YAML 文件内容
        outFile << "body_T_cam" << i << ": !!opencv-matrix\n";
        outFile << "   rows: 4\n";
        outFile << "   cols: 4\n";
        outFile << "   dt: d\n";
        outFile << "   data: [ ";

        // 输出矩阵数据
        for (int j = 0; j < 4; ++j) {
            outFile << "\n      "; // 换行并缩进
            for (int k = 0; k < 4; ++k) {
                outFile << T(j, k);
                if (k != 3) outFile << ", ";
            }
            if (j != 3) outFile << ",";
        }

        outFile << "]\n\n"; // 添加换行以分隔不同的相机
    }

    outFile << "Pbg: !!opencv-matrix\n";
    outFile << "   rows: 3\n";
    outFile << "   cols: 1\n";
    outFile << "   dt: d\n";
    outFile << "   data: [ ";

    // 输出平移向量数据
    for (int i = 0; i < 3; ++i) {
        outFile << Pbg(i);
        if (i != 2) outFile << ", ";
    }

    outFile << " ]\n";

    outFile.close();
}

void printStatistics(const SWFOptimization& swf_optimization, double t) {

    if (solver_flag != NonLinear)
        return;

    int index = PUB_INDEX;
    LOG_OUT <<  "time:" << swf_optimization.headers[index] << std::endl;
    cout << "pos: " << (Rwgw.transpose() * swf_optimization.Ps[index]).transpose() << " " << (Rwgw.transpose() * (R_WI_WC.transpose() * (swf_optimization.Ps[index] - P_WI_WC))).transpose() << std::endl;
    LOG_OUT << "vel: " << (swf_optimization.Vs[index]).transpose() << std::endl;
    LOG_OUT << "orientation: " << Utility::R2ypr(swf_optimization.Rs[index]).transpose() << std::endl;
    LOG_OUT << "gyro bias: " << swf_optimization.Bgs[index].transpose() << std::endl;
    LOG_OUT << "acc bias: " << swf_optimization.Bas[index].transpose() << std::endl;
    if (ESTIMATE_TD) LOG_OUT << "td: " << swf_optimization.td << std::endl;
    if (ESTIMATE_EXTRINSIC)LOG_OUT << "extrinsic: " << TIC[0].transpose() << ",\t" << Utility::R2ypr(RIC[0]).transpose() << std::endl;
#if ESTIMATE_EXTRINSIC2&&USE_STEREO
    LOG_OUT << "extrinsic2: " << TIC[1].transpose() << ",\t" << Utility::R2ypr(RIC[1]).transpose() << std::endl;
#endif
    if ( ESTIMATE_ACC_SCALE) LOG_OUT << "imu scale: " << swf_optimization.acc_scale.transpose() << std::endl;
    if ( ESTIMATE_PBG)LOG_OUT << "Pbg: " << Pbg.transpose() << std::endl;
    LOG_OUT << "RWGW: " << Utility::R2ypr(Rwgw).transpose() << std::endl;
    LOG_OUT << "R_WI_WC: " << Utility::R2ypr(R_WI_WC).transpose() << std::endl;
    LOG_OUT << "P_WI_WC:" << P_WI_WC.transpose() << std::endl;
    LOG_OUT << "image_count:" << swf_optimization.image_count << std::endl;

    if (ESTIMATE_PBG || (ESTIMATE_EXTRINSIC2 && USE_STEREO) || ESTIMATE_EXTRINSIC)
        saveMultipleExtrinsicsToYaml(TIC, QIC, Pbg, RESULT_PATH + "extrinsic.yaml");

}

void resetpot(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {
    // Ps0 = swf_optimization.Ps[0];
}

void pubOdometry(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {
    std::vector<Vector3d>        Ps(swf_optimization.image_count);


    for (int i = 0; i < swf_optimization.image_count; i++)
        Ps[i] = R_WI_WC.transpose() * (swf_optimization.Ps[i] - P_WI_WC);
    if (solver_flag == NonLinear) {
        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(Rwgw.transpose() * R_WI_WC.transpose() * swf_optimization.Rs[PUB_INDEX]);
        Vector3d tmpP = Rwgw.transpose() * Ps[PUB_INDEX];

        Vector3d tmpv = swf_optimization.Vs[PUB_INDEX];



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
        if (PUB_INDEX < 0)
            return;


        {
            nav_msgs::Odometry odometry;
            odometry.header = header;
            odometry.header.frame_id = "world";
            odometry.pose.pose.position.x = TIC[0].x();
            odometry.pose.pose.position.y = TIC[0].y();
            odometry.pose.pose.position.z = TIC[0].z();
            Quaterniond tmp_q{RIC[0]};
            odometry.pose.pose.orientation.x = tmp_q.x();
            odometry.pose.pose.orientation.y = tmp_q.y();
            odometry.pose.pose.orientation.z = tmp_q.z();
            odometry.pose.pose.orientation.w = tmp_q.w();
            pub_extrinsic.publish(odometry);
        }
    }
}


void pubCameraPose(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {

    std::vector<Vector3d>        Ps(swf_optimization.image_count);


    for (int i = 0; i < swf_optimization.image_count; i++)
        Ps[i] = R_WI_WC.transpose() * (swf_optimization.Ps[i] - P_WI_WC);

    int idx2 = PUB_INDEX;
    if (solver_flag == NonLinear) {
        int i = idx2;
        Vector3d P = Rwgw.transpose() * (Ps[i] + swf_optimization.Rs[i] * TIC[0]);
        Quaterniond R = Quaterniond((Rwgw.transpose() * R_WI_WC.transpose() * swf_optimization.Rs[i] * RIC[0]));
#if ENABLE_STEREO_MAP
#if REAL_TIME
        mutex_pose.lock();
#endif
        PubRs[header.stamp.toSec()] = R.toRotationMatrix().cast<float>();
        Pubts[header.stamp.toSec()] = P.cast<float>();
#if REAL_TIME
        mutex_pose.unlock();
#endif
#endif
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
#if USE_STEREO
        {
            Vector3d P_r =  Ps[i] + swf_optimization.Rs[i] * TIC[1];
            Quaterniond R_r = Quaterniond(swf_optimization.Rs[i] * RIC[1]);

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
#endif
        pub_camera_pose.publish(odometry);

        cameraposevisual.reset();
        cameraposevisual.add_pose(P, R);
#if USE_STEREO
        {
            Vector3d P = Rwgw.transpose() * (Ps[i] + swf_optimization.Rs[i] * TIC[1]);
            Quaterniond R = Quaterniond((Rwgw.transpose() * R_WI_WC.transpose() * swf_optimization.Rs[i] * RIC[1]));

            cameraposevisual.add_pose(P, R);
        }
#endif
        cameraposevisual.publish_by(pub_camera_pose_visual, odometry.header);
    }
}


void pubPointCloud(const SWFOptimization& swf_optimization, const std_msgs::Header& header) {
    std::vector<Vector3d>        Ps(swf_optimization.image_count);


    for (int i = 0; i < swf_optimization.image_count; i++)
        Ps[i] = R_WI_WC.transpose() * (swf_optimization.Ps[i] - P_WI_WC);

    sensor_msgs::PointCloud point_cloud;
    point_cloud.header = header;

    for (auto& it_per_id : swf_optimization.f_manager.feature) {
        if (!it_per_id.valid)continue;
        if (it_per_id.start_frame > SWF_SIZE_IN * (SWF_SIZE_OUT - 1))
            continue;
        Eigen::Vector3d ptsInW;

        bool is_cross = (it_per_id.start_frame / SWF_SIZE_IN != (/*it_per_id.endFrame()*/it_per_id.start_frame + (int)it_per_id.feature_per_frame.size() - 1/*it_per_id.endFrame()*/ - 1) / SWF_SIZE_IN) || it_per_id.start_frame == 0;
        if (is_cross) {

            int imu_i = it_per_id.start_frame;

            if  (it_per_id.start_frame % SWF_SIZE_IN != 0)
                imu_i = (it_per_id.start_frame / SWF_SIZE_IN + 1) * SWF_SIZE_IN;

            if (idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id] < 0 || 1 / idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id] > 20)continue;
            ASSERT(idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id] != 0);
            ptsInW = swf_optimization.Rs[imu_i] * (
                         RIC[0] * (it_per_id.feature_per_frame[imu_i - it_per_id.start_frame].point /
                                   idepths_all[imu_i / SWF_SIZE_IN][it_per_id.feature_id]) + TIC[0]
                     ) + swf_optimization.Ps[imu_i];
        } else {
            if (idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] < 0 || 1 / idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] > 20)continue;
            ptsInW = swf_optimization.Rs[it_per_id.start_frame] * (
                         RIC[0] * (it_per_id.feature_per_frame[0].point /
                                   idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id]) + TIC[0]
                     ) + swf_optimization.Ps[it_per_id.start_frame];
        }

        Vector3d w_pts_i = Rwgw.transpose() * (R_WI_WC.transpose() * (ptsInW - P_WI_WC));

        geometry_msgs::Point32 p;
        p.x = w_pts_i(0);
        p.y = w_pts_i(1);
        p.z = w_pts_i(2);

        point_cloud.points.push_back(p);


    }
    pub_point_cloud_short.publish(point_cloud);

}