

#ifndef FEATURE_MANAGER_H
#define FEATURE_MANAGER_H


using namespace std;
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>

#include "../parameter/parameters.h"
#include "../factor/visual_inerial_factor.h"

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
};

class FeaturePerId {
  public:
    FeaturePerId(int _feature_id, int _start_frame)
        : feature_id(_feature_id), start_frame(_start_frame), valid(false), solve_flag(0) {
    }
    void RemoveIdepth() {
        for (int i = start_frame / SWF_SIZE_IN; i < endFrame() / SWF_SIZE_IN + 1; i++) {
            if (idepths_all[i].find(feature_id) != idepths_all[i].end())idepths_all[i].erase(feature_id);
        }
    }
    int endFrame() {
        return start_frame + feature_per_frame.size() - 1;
    }

    int feature_id;
    int start_frame;
    bool valid;
    int solve_flag; // 0 haven't solve yet; 1 solve succ; 2 solve fail;

    vector<FeaturePerFrame> feature_per_frame;

};

class FeatureManager {
  public:
    FeatureManager();

    void ClearState();
    void addFeature(int image_index, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image);
    vector<pair<Eigen::Vector3d, Eigen::Vector3d>> getCorresponding(int frame_count_l, int frame_count_r);
    void removeFailures(std::vector<VisualInertialBase*>& visual_inertial_bases);
    void clearDepth();
    void triangulate(std::vector<Eigen::Vector3d>& Ps, std::vector<Eigen::Matrix3d>& Rs);
    void triangulatePoint(Eigen::Matrix<double, 3, 4>& Pose0, Eigen::Matrix<double, 3, 4>& Pose1,
                          Eigen::Vector2d& point0, Eigen::Vector2d& point1, Eigen::Vector3d& point_3d);
    void removeBack(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0);
    void removeFront(int image_index, Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0);
    void removeOut(int windowsize);
    bool CheckParallax(int);
    double compensatedParallax(const FeaturePerId& it_per_id, int image_index);
    void GetFailurePointers(std::set<double*>& feature_marge_pointer, std::set<double*>& prior_pointer);
    void removeFeatures(std::set<int>feature_ids);

    list<FeaturePerId> feature;
    std::unordered_map<int, FeaturePerId*>feature_map;


};

#endif
