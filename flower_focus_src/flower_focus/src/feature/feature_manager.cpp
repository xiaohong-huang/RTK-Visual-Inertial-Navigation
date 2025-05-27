

#include "feature_manager.h"
#include "../parameter/parameters.h"





FeatureManager::FeatureManager() {
}



void FeatureManager::ClearState() {
    feature.clear();
}




void FeatureManager::addFeature(int image_index, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& image) {


    for (auto& id_pts : image) {
        FeaturePerFrame f_per_fra(id_pts.second[0].second);
        ASSERT(id_pts.second[0].first == 0);
        if (id_pts.second.size() == 2) {
            f_per_fra.rightObservation(id_pts.second[1].second);
            ASSERT(id_pts.second[1].first == 1);
        }

        int feature_id = id_pts.first;
        ASSERT(feature_id >= 0);
        FeaturePerId* it = 0;
        if (feature_map.find(feature_id) != feature_map.end()) {
            it = feature_map[feature_id];
            ASSERT(it->feature_id == feature_id);
        }

        if (it == 0) {
            feature.push_back(FeaturePerId(feature_id, image_index));
            feature.back().feature_per_frame.push_back(f_per_fra);
            feature_map[feature_id] = &feature.back();
        } else if (it->feature_id == feature_id)
            it->feature_per_frame.push_back(f_per_fra);

    }



}



bool FeatureManager::CheckParallax(int image_index) {
    double parallax_sum = 0;
    int parallax_num = 0;

    int last_track_num = 0;
    int new_feature_num = 0;
    int long_track_num = 0;

    for (auto& it_per_id : feature) {
        if (it_per_id.start_frame <= image_index - 2 &&
                it_per_id.start_frame + int(it_per_id.feature_per_frame.size()) - 1 >= image_index - 1) {
            parallax_sum += compensatedParallax(it_per_id, image_index);
            parallax_num++;
        }
        if (it_per_id.endFrame() == image_index) {
            if (it_per_id.feature_per_frame.size() == 1)new_feature_num++;
            if (it_per_id.feature_per_frame.size() >= 2)last_track_num++;
            if (it_per_id.feature_per_frame.size() >= 4)long_track_num++;
        }
    }

    if (image_index < 2 || last_track_num < 20 || long_track_num < 40 || new_feature_num > 0.5 * last_track_num)
        return true;

    if (parallax_num == 0)
        return true;
    else {
        if (have_hist)assert(image_index > 10);
        return parallax_sum / parallax_num >= (image_index > 10 ? MIN_PARALLAX/FOCAL_LENGTH : 30/FOCAL_LENGTH);
    }
}

vector<pair<Eigen::Vector3d, Eigen::Vector3d>> FeatureManager::getCorresponding(int frame_count_l, int frame_count_r) {
    vector<pair<Eigen::Vector3d, Eigen::Vector3d>> corres;
    for (auto& it : feature) {
        if (it.start_frame <= frame_count_l && it.endFrame() >= frame_count_r) {
            Eigen::Vector3d a = Eigen::Vector3d::Zero(), b = Eigen::Vector3d::Zero();
            int idx_l = frame_count_l - it.start_frame;
            int idx_r = frame_count_r - it.start_frame;

            a = it.feature_per_frame[idx_l].point;

            b = it.feature_per_frame[idx_r].point;

            corres.push_back(make_pair(a, b));
        }
    }
    return corres;
}
void FeatureManager::GetFailurePointers(std::set<double*>& feature_marge_pointer, std::set<double*>& prior_pointer) {
    for (auto it = feature.begin(), it_next = feature.begin(); it != feature.end(); it = it_next) {
        it_next++;
        if (it->solve_flag == 2) {
            if  (idepths_all[0].find(it->feature_id) != idepths_all[0].end() && prior_pointer.find(&idepths_all[0][it->feature_id]) != prior_pointer.end()) {
                feature_marge_pointer.insert(&idepths_all[0][it->feature_id]);
                ASSERT(it->start_frame == 0);
            }
        }
    }
}

void FeatureManager::removeFailures( std::vector<VisualInertialBase*>& visual_inertial_bases_global) {
    for (auto it = feature.begin(), it_next = feature.begin(); it != feature.end(); it = it_next) {
        it_next++;
        if (it->solve_flag == 2) {

            for (int i = 0; i < (int)visual_inertial_bases_global.size(); i++) {
                if (visual_inertial_bases_global[i])visual_inertial_bases_global[i]->RemoveFeature(it->feature_id);
            }
            it->RemoveIdepth();
            feature_map.erase(it->feature_id);
            feature.erase(it);

        }
    }
}

void FeatureManager::removeFeatures(std::set<int>feature_ids) {
    for (auto it = feature.begin(), it_next = feature.begin();
            it != feature.end(); it = it_next) {
        it_next++;
        if (feature_ids.find(it->feature_id) != feature_ids.end()) {
            it->RemoveIdepth();
            feature_map.erase(it->feature_id);
            feature.erase(it);
        }

    }
}


void FeatureManager::clearDepth() {
    for (auto& it_per_id : feature)
        it_per_id.valid = false;
}


void FeatureManager::triangulatePoint(Eigen::Matrix<double, 3, 4>& Pose0, Eigen::Matrix<double, 3, 4>& Pose1,
                                      Eigen::Vector2d& point0, Eigen::Vector2d& point1, Eigen::Vector3d& point_3d) {
    Eigen::Matrix4d design_matrix = Eigen::Matrix4d::Zero();
    design_matrix.row(0) = point0[0] * Pose0.row(2) - Pose0.row(0);
    design_matrix.row(1) = point0[1] * Pose0.row(2) - Pose0.row(1);
    design_matrix.row(2) = point1[0] * Pose1.row(2) - Pose1.row(0);
    design_matrix.row(3) = point1[1] * Pose1.row(2) - Pose1.row(1);
    Eigen::Vector4d triangulated_point;
    triangulated_point =
        design_matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();
    point_3d(0) = triangulated_point(0) / triangulated_point(3);
    point_3d(1) = triangulated_point(1) / triangulated_point(3);
    point_3d(2) = triangulated_point(2) / triangulated_point(3);
}


void FeatureManager::triangulate(std::vector<Eigen::Vector3d>& Ps, std::vector<Eigen::Matrix3d>& Rs) {
    for (auto& it_per_id : feature) {
        if (it_per_id.feature_per_frame.size() < FEATURE_CONTINUE)continue;
        if (it_per_id.valid) continue;
        #if USE_STEREO
        if ( it_per_id.feature_per_frame[0].is_stereo&&STEREO_SCALE) {
            int imu_i = it_per_id.start_frame;
            Eigen::Matrix<double, 3, 4> leftPose;
            Eigen::Vector3d t0 = Ps[imu_i] + Rs[imu_i] * TIC[0];
            Eigen::Matrix3d R0 = Rs[imu_i] * RIC[0];
            leftPose.leftCols<3>() = R0.transpose();
            leftPose.rightCols<1>() = -R0.transpose() * t0;

            Eigen::Matrix<double, 3, 4> rightPose;
            Eigen::Vector3d t1 = Ps[imu_i] + Rs[imu_i] * TIC[1];
            Eigen::Matrix3d R1 = Rs[imu_i] * RIC[1];
            rightPose.leftCols<3>() = R1.transpose();
            rightPose.rightCols<1>() = -R1.transpose() * t1;

            Eigen::Vector2d point0, point1;
            Eigen::Vector3d point3d;
            point0 = it_per_id.feature_per_frame[0].point.head(2);
            point1 = it_per_id.feature_per_frame[0].pointRight.head(2);

            triangulatePoint(leftPose, rightPose, point0, point1, point3d);
            Eigen::Vector3d localPoint;
            localPoint = leftPose.leftCols<3>() * point3d + leftPose.rightCols<1>();
            double depth = localPoint.z();
            if (depth <= 0)
                depth = INIT_DEPTH;
            idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] = 1. / depth;
            it_per_id.valid = true;
            continue;
        } else
#endif
        if (it_per_id.feature_per_frame.size() > 1) {
            int imu_i = it_per_id.start_frame;
            Eigen::Matrix<double, 3, 4> leftPose;
            Eigen::Vector2d point0, point1;
            Eigen::Vector3d t0 = Ps[imu_i] + Rs[imu_i] * TIC[0];
            Eigen::Matrix3d R0 = Rs[imu_i] * RIC[0];
            leftPose.leftCols<3>() = R0.transpose();
            leftPose.rightCols<1>() = -(R0.transpose() * t0);
            point0 = it_per_id.feature_per_frame[imu_i - it_per_id.start_frame].point.head(2);

            imu_i = it_per_id.endFrame();
            Eigen::Matrix<double, 3, 4> rightPose;
            Eigen::Vector3d t1 = Ps[imu_i] + Rs[imu_i] * TIC[0];
            Eigen::Matrix3d R1 = Rs[imu_i] * RIC[0];
            rightPose.leftCols<3>() = R1.transpose();
            rightPose.rightCols<1>() = -(R1.transpose() * t1);
            point1 = it_per_id.feature_per_frame[imu_i - it_per_id.start_frame].point.head(2);


            Eigen::Vector3d point3d;

            triangulatePoint(leftPose, rightPose, point0, point1, point3d);
            Eigen::Vector3d localPoint;
            localPoint = leftPose.leftCols<3>() * point3d + leftPose.rightCols<1>();
            double depth = localPoint.z();
            if (depth <= 0)
                depth = INIT_DEPTH;
            idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] = 1. / depth;
            it_per_id.valid = true;
            continue;
        }
        ASSERT(0);
        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;

        Eigen::MatrixXd svd_A(2 * it_per_id.feature_per_frame.size(), 4);
        int svd_idx = 0;

        Eigen::Matrix<double, 3, 4> P0;
        Eigen::Vector3d t0 = Ps[imu_i] + Rs[imu_i] * TIC[0];
        Eigen::Matrix3d R0 = Rs[imu_i] * RIC[0];
        P0.leftCols<3>() = Eigen::Matrix3d::Identity();
        P0.rightCols<1>() = Eigen::Vector3d::Zero();

        for (auto& it_per_frame : it_per_id.feature_per_frame) {
            imu_j++;

            Eigen::Vector3d t1 = Ps[imu_j] + Rs[imu_j] * TIC[0];
            Eigen::Matrix3d R1 = Rs[imu_j] * RIC[0];
            Eigen::Vector3d t = R0.transpose() * (t1 - t0);
            Eigen::Matrix3d R = R0.transpose() * R1;
            Eigen::Matrix<double, 3, 4> P;
            P.leftCols<3>() = R.transpose();
            P.rightCols<1>() = -R.transpose() * t;
            Eigen::Vector3d f = it_per_frame.point.normalized();
            svd_A.row(svd_idx++) = f[0] * P.row(2) - f[2] * P.row(0);
            svd_A.row(svd_idx++) = f[1] * P.row(2) - f[2] * P.row(1);

            if (imu_i == imu_j)
                continue;
        }
        ASSERT(svd_idx == svd_A.rows());
        Eigen::Vector4d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(svd_A, Eigen::ComputeThinV).matrixV().rightCols<1>();
        double depth = svd_V[2] / svd_V[3];
        if (depth < 0.1)
            depth = INIT_DEPTH;
        idepths_all[it_per_id.start_frame / SWF_SIZE_IN][it_per_id.feature_id] = 1. / depth;
        it_per_id.valid = true;
    }
}


void FeatureManager::removeBack(Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0) {

    for (auto it = feature.begin(), it_next = feature.begin();
            it != feature.end(); it = it_next) {

        it_next++;
        bool drop_flag = false;
        if (it->start_frame >= SWF_SIZE_IN)
            it->start_frame -= SWF_SIZE_IN;
        else {
            for (int i = 0; i < SWF_SIZE_IN - it->start_frame; i++) {

                it->feature_per_frame.erase(it->feature_per_frame.begin());
                if (it->feature_per_frame.size() <= 1) {
                    it->RemoveIdepth();
                    feature_map.erase(it->feature_id);
                    feature.erase(it);
                    drop_flag = true;
                    break;
                }
            }
            if (!drop_flag)
                it->start_frame = 0;
        }
    }

    for (int i = 0; i < (int)idepths_all.size() - 1; i++)
        idepths_all[i].swap(idepths_all[i + 1]);
    idepths_all[idepths_all.size() - 1].clear();



}


void FeatureManager::removeFront(int image_index, Eigen::Vector3d P0, Eigen::Matrix3d R0, Eigen::Vector3d P1, Eigen::Matrix3d R1, Eigen::Vector3d tic0, Eigen::Matrix3d ric0) {

    int marg_index = image_index - 1;
    for (auto it = feature.begin(), it_next = feature.begin(); it != feature.end(); it = it_next) {
        it_next++;

        if (marg_index % SWF_SIZE_IN == 0
                && idepths_all[marg_index / SWF_SIZE_IN].find(it->feature_id) != idepths_all[marg_index / SWF_SIZE_IN].end()
                && idepths_all[marg_index / SWF_SIZE_IN][it->feature_id]
                && it->valid) {
            Eigen::Vector3d ptsInW = R0 * (  ric0 * (it->feature_per_frame[marg_index - it->start_frame].point / idepths_all[marg_index / SWF_SIZE_IN][it->feature_id]) + tic0 ) + P0;
            Eigen::Vector3d pts_cj = ric0.transpose() * ( R1.transpose() * (ptsInW - P1) - tic0);
            idepths_all[marg_index / SWF_SIZE_IN][it->feature_id] = 1.0 / pts_cj.z();
        }

        if (it->start_frame == image_index)
            it->start_frame--;
        else {
            int j = marg_index - it->start_frame;
            if (it->endFrame() < marg_index)
                continue;
            it->feature_per_frame.erase(it->feature_per_frame.begin() + j);

            if (it->feature_per_frame.size() == 0) {
                it->RemoveIdepth();
                feature_map.erase(it->feature_id);
                feature.erase(it);
            }

        }

    }
}


void FeatureManager::removeOut(int windowsize) {
    for (auto it = feature.begin(), it_next = feature.begin(); it != feature.end(); it = it_next) {
        it_next++;
        if (it->endFrame() != windowsize - 1 && !it->valid) {
            it->RemoveIdepth();
            feature_map.erase(it->feature_id);
            feature.erase(it);
        }
    }
}


double FeatureManager::compensatedParallax(const FeaturePerId& it_per_id, int image_index) {
    //check the second last frame is keyframe or not
    //parallax betwwen seconde last frame and third last frame
    const FeaturePerFrame& frame_i = it_per_id.feature_per_frame[image_index - 2 - it_per_id.start_frame];
    const FeaturePerFrame& frame_j = it_per_id.feature_per_frame[image_index - 1 - it_per_id.start_frame];

    double ans = 0;
    Eigen::Vector3d p_j = frame_j.point;

    double u_j = p_j(0);
    double v_j = p_j(1);

    Eigen::Vector3d p_i = frame_i.point;
    Eigen::Vector3d p_i_comp;

    p_i_comp = p_i;
    double dep_i = p_i(2);
    double u_i = p_i(0) / dep_i;
    double v_i = p_i(1) / dep_i;
    double du = u_i - u_j, dv = v_i - v_j;

    double dep_i_comp = p_i_comp(2);
    double u_i_comp = p_i_comp(0) / dep_i_comp;
    double v_i_comp = p_i_comp(1) / dep_i_comp;
    double du_comp = u_i_comp - u_j, dv_comp = v_i_comp - v_j;

    ans = max(ans, sqrt(min(du * du + dv * dv, du_comp * du_comp + dv_comp * dv_comp)));

    return ans;
}
