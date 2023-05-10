

#ifndef FEATURE_MANAGER_H
#define FEATURE_MANAGER_H


using namespace std;
#include "ceres/problem.h"
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>



class FeaturePerFrame {
  public:
    FeaturePerFrame(const Eigen::Matrix<double, 7, 1>& _point) {
        point.x() = _point(0);
        point.y() = _point(1);
        point.z() = _point(2);
        uv.x() = _point(3);
        uv.y() = _point(4);
        velocity.x() = _point(5);
        velocity.y() = _point(6);
        is_stereo = false;
        is_in_optimize = false;
    }
    void rightObservation(const Eigen::Matrix<double, 7, 1>& _point) {
        pointRight.x() = _point(0);
        pointRight.y() = _point(1);
        pointRight.z() = _point(2);
        uvRight.x() = _point(3);
        uvRight.y() = _point(4);
        velocityRight.x() = _point(5);
        velocityRight.y() = _point(6);
        is_stereo = true;
    }
    Eigen::Vector3d point, pointRight;
    Eigen::Vector2d uv, uvRight;
    Eigen::Vector2d velocity, velocityRight;
    bool is_stereo;
    bool is_in_optimize;
};

class FeaturePerId {
  public:
    int feature_id;
    int start_frame;
    vector<FeaturePerFrame> feature_per_frame;
    int used_num;
    bool valid;
    Eigen::Vector3d ptsInWorld;
    double idepth_;

    int solve_flag; // 0 haven't solve yet; 1 solve succ; 2 solve fail;

    FeaturePerId(int _feature_id, int _start_frame)
        : feature_id(_feature_id), start_frame(_start_frame),
          used_num(0), valid(false), solve_flag(0) {
    }

    int endFrame();
};

class FeatureManager {
  public:
    FeatureManager(Eigen::Matrix3d _Rs[]);

    void setRic(Eigen::Matrix3d _ric[]);
    void ClearState();
    int getFeatureCount();
    bool addFeatureCheckParallax(int image_index, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image);
    vector<pair<Eigen::Vector3d, Eigen::Vector3d>> getCorresponding(int frame_count_l, int frame_count_r);
    //void updateDepth(const VectorXd &x);
    void removeFailures(ceres::Problem& my_problem);
    void clearDepth();
    void triangulate(Eigen::Vector3d Ps[], Eigen::Matrix3d Rs[], Eigen::Vector3d tic[], Eigen::Matrix3d ric[], Eigen::Vector3d Pbg_);
    void triangulatePoint(Eigen::Matrix<double, 3, 4>& Pose0, Eigen::Matrix<double, 3, 4>& Pose1,
                          Eigen::Vector2d& point0, Eigen::Vector2d& point1, Eigen::Vector3d& point_3d);
    void initFramePoseByPnP(int frameCnt, Eigen::Vector3d Ps[], Eigen::Matrix3d Rs[], Eigen::Vector3d tic[], Eigen::Matrix3d ric[], Eigen::Vector3d Pbg_);
    bool solvePoseByPnP(Eigen::Matrix3d& R_initial, Eigen::Vector3d& P_initial,
                        vector<cv::Point2f>& pts2D, vector<cv::Point3f>& pts3D);
    void removeBack(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0, Eigen::Vector3d Pbg_, ceres::Problem& my_problem);
    void removeFront(int image_index, ceres::Problem& my_problem);
    void removeOut(int windowsize, ceres::Problem& my_problem);
    void removeOut2(int windowsize, ceres::Problem& my_problem);
    bool CheckParallax(int);
    list<FeaturePerId> feature;
    int last_track_num;
    double last_average_parallax;
    int new_feature_num;
    int long_track_num;
    int back_up_feature_index = -1;

  private:
    double compensatedParallax2(const FeaturePerId& it_per_id, int image_index);
    const Eigen::Matrix3d* Rs;
    Eigen::Matrix3d ric[2];
};

#endif
