#include "swf.h"
#include "../factor/mag_factor.h"
#include "../utility/visualization.h"

void SWFOptimization::InitializePos(Matrix3d& Rwgw_) {

    for (int i = 0; i < FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1; i++) {
        Bgs[i] = gyr_mean;
    }

    std::cout << "initial bgs:" << Bgs[0];
    printf("averge acc %f %f %f\n", acc_mean.x(), acc_mean.y(), acc_mean.z());

    if (mag_mean.norm() < 1e-2) {
        assert(0);
        mag_mean(0) = 0;
        mag_mean(1) = 1;
        mag_mean(2) = 0;
    }
    Matrix3d Rwb0;
    Eigen::Vector3d z0 = acc_mean.normalized();
    Eigen::Vector3d x0 = (Utility::skewSymmetric(mag_mean) * z0).normalized();
    Eigen::Vector3d y0 = (Utility::skewSymmetric(z0) * x0).normalized();
    Rwb0.block(0, 0, 1, 3) = x0.transpose();
    Rwb0.block(1, 0, 1, 3) = y0.transpose();
    Rwb0.block(2, 0, 1, 3) = z0.transpose();

    Rwgw = Rwgw_;

    for (int i = 0; i < FEATURE_WINDOW_SIZE + GNSS_WINDOW_SIZE + 1; i++) {
        Rs[i] = Rwgw * Rwb0;
    }
    cout << "init R0: " << endl << Utility::R2ypr(Rs[0]).transpose() << "," << Utility::R2ypr(Rwb0).transpose() << "," << Utility::R2ypr(Rwgw).transpose() << "  " << endl;

    Pbg = Pbg;

}


void SWFOptimization::InputIMU(double t, const Vector3d& linearAcceleration, const Vector3d& angularVelocity) {
    static int gyrcount0;
    if (gyrcount0 < SKIP_IMU) {
        gyrcount0++;
        return;
    }

    if (gyr_count < AVERAGE_IMU) {
        acc_mean += linearAcceleration;
        gyr_mean += angularVelocity;
        gyr_count++;
        return;
    }

    if (gyr_count >= AVERAGE_IMU && !imu_initialize && (!USE_MAG_INIT_YAW || !mag_buf.empty())) {
        imu_initialize = true;
        acc_mean /= gyr_count;
        gyr_mean /= gyr_count;
        InitializePos(Rwgw);
    }


    acc_buf.push(make_pair(t, linearAcceleration));
    gyr_buf.push(make_pair(t, angularVelocity));



}


void SWFOptimization::InputMag(double t, const Vector3d& Magdata) {

    if (USE_MAG_INIT_YAW) {
        mag_mean = Magdata;
        mag_buf.push_back(make_pair(t, Magdata));
    }



}


bool SWFOptimization::GetImuInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>>& accVector,
                                     vector<pair<double, Eigen::Vector3d>>& gyrVector) {
    if (acc_buf.empty()) {
        printf("not receive imu\n");
        return false;
    }
    if (t1 <= acc_buf.back().first) {
        while (acc_buf.front().first <= t0) {
            acc_buf.pop();
            gyr_buf.pop();
        }
        while (acc_buf.front().first < t1) {
            accVector.push_back(acc_buf.front());
            acc_buf.pop();
            gyrVector.push_back(gyr_buf.front());
            gyr_buf.pop();
        }
        accVector.push_back(acc_buf.front());
        gyrVector.push_back(gyr_buf.front());
    } else {
        printf("wait for imu\n");
        return false;
    }
    return true;
}


bool SWFOptimization::ImuAvailable(double t) {
    if (!acc_buf.empty() && t <= acc_buf.back().first)
        return true;
    else
        return false;
}


void SWFOptimization::ImuIntegrate() {

    vector<pair<double, Eigen::Vector3d>> accVector, gyrVector;
    if (USE_IMU)GetImuInterval(prev_time, cur_time, accVector, gyrVector);
    assert(accVector.size() > 1);
    if (USE_IMU) {
        Ps[rover_count + image_count - 1] = Ps[rover_count + image_count - 1] - Rs[rover_count + image_count - 1] * Pbg;
        if (rover_count + image_count - 1 > 5)
            Vs[rover_count + image_count - 1] = Vs[rover_count + image_count - 1] - Rs[rover_count + image_count - 1] * (Utility::skewSymmetric(pre_integrations[rover_count + image_count - 2]->gyrj - Bgs[rover_count + image_count - 1]) * Pbg);
        assert(prev_time == prev_time2);
        assert(cur_time >= prev_time);
        double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0;
        for (size_t i = 0; i < accVector.size(); i++) {
            static double current_time = -1;
            double t = accVector[i].first;
            double img_t = cur_time;

            if (t <= img_t) {
                if (current_time < 0)
                    current_time = t;
                double dt = t - current_time;
                assert(dt >= 0);
                current_time = t;

                dx = accVector[i].second.x();
                dy = accVector[i].second.y();
                dz = accVector[i].second.z();
                rx = gyrVector[i].second.x();
                ry = gyrVector[i].second.y();
                rz = gyrVector[i].second.z();
                headers[image_count + rover_count - 1] = t;
                IMUProcess(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
            } else {
                assert(dx + dy + dz + rx + ry + rz != 0);
                double dt_1 = img_t - current_time;
                double dt_2 = t - img_t;
                current_time = img_t;
                assert(dt_1 >= 0);
                assert(dt_2 >= 0);
                assert(dt_1 + dt_2 > 0);
                double w1 = dt_2 / (dt_1 + dt_2);
                double w2 = dt_1 / (dt_1 + dt_2);
                dx = w1 * dx + w2 * accVector[i].second.x();
                dy = w1 * dy + w2 * accVector[i].second.y();
                dz = w1 * dz + w2 * accVector[i].second.z();
                rx = w1 * rx + w2 * gyrVector[i].second.x();
                ry = w1 * ry + w2 * gyrVector[i].second.y();
                rz = w1 * rz + w2 * gyrVector[i].second.z();
                headers[image_count + rover_count - 1] = img_t;
                IMUProcess(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));

                // printf("dimu: dt:%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",dt_1,dt_2,img_t,w1,w2,dx,dy,dz,rx,ry,rz);
            }
        }

        Ps[rover_count + image_count - 1] = Ps[rover_count + image_count - 1] + Rs[rover_count + image_count - 1] * Pbg;
        if (rover_count + image_count - 1 > 5)
            Vs[rover_count + image_count - 1] = Vs[rover_count + image_count - 1] + Rs[rover_count + image_count - 1] * (Utility::skewSymmetric(pre_integrations[rover_count + image_count - 1]->gyrj - Bgs[rover_count + image_count - 1]) * Pbg);
        prev_time = cur_time;
    }
}

void SWFOptimization::IMUProcess( double dt, const Vector3d& linear_acceleration, const Vector3d& angular_velocity) {
    if (first_imu) {
        first_imu = false;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }
    if (!pre_integrations[image_count + rover_count - 1]) {
        pre_integrations[image_count + rover_count - 1] = new IntegrationBase{acc_0, gyr_0, Bas[image_count + rover_count - 1], Bgs[image_count + rover_count - 1]};
        if (image_count + rover_count - 2 >= 0 && USE_GLOBAL_OPTIMIZATION) {
            IMUFactor* factor = new IMUFactor(pre_integrations[image_count + rover_count - 1]);
            my_problem.AddResidualBlock(factor, 0, para_pose[image_count + rover_count - 2],
                                        para_speed_bias[image_count + rover_count - 2],
                                        para_pose[image_count + rover_count - 1],
                                        para_speed_bias[image_count + rover_count - 1]);
        }

    }
    if (image_count + rover_count - 1 != 0) {
        pre_integrations[image_count + rover_count - 1]->push_back(dt, linear_acceleration, angular_velocity);
        dt_buf[image_count + rover_count - 1].push_back(dt);
        linear_acceleration_buf[image_count + rover_count - 1].push_back(linear_acceleration);
        angular_velocity_buf[image_count + rover_count - 1].push_back(angular_velocity);
        int j = image_count + rover_count - 1;
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - Rwgw * G;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - Rwgw * G;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
    save_result(*this);
}





void SWFOptimization::MagProcess(double curtime) {
    mags[rover_count + image_count - 1].setZero();
    if (mag_buf.size() > 0) {
        auto selectmagit = mag_buf.begin();
        if (mag_buf.size() > 1) {
            auto it = mag_buf.begin();
            it++;
            for (; it != mag_buf.end(); it++) {
                if (abs(it->first - curtime) < abs(selectmagit->first - curtime)) {
                    selectmagit = it;
                }
            }
        }
        mags[rover_count + image_count - 1] = selectmagit->second;
        for (auto it = mag_buf.begin(), it_next = mag_buf.begin(); it != mag_buf.end(); it = it_next) {
            it_next++;
            if (curtime - it->first > 10)mag_buf.erase(it);
        }
    }

    if (USE_MAG_CORRECT_YAW && USE_IMU && USE_GLOBAL_OPTIMIZATION) {
        MAGFactor* factor = new MAGFactor(mags[rover_count + image_count - 1]);
        my_problem.AddResidualBlock(factor, 0, para_pose[rover_count + image_count - 1], para_bmg);
    }

    Eigen::Vector3d z0 = (Rwgw.transpose() * Rs[rover_count + image_count - 1]).transpose() * G;
    Eigen::Vector3d x0 = (Utility::skewSymmetric(mags[rover_count + image_count - 1] + Eigen::Vector3d(para_bmg)) * z0).normalized();
    Eigen::Vector3d y0 = (Utility::skewSymmetric(z0) * x0).normalized();
    Matrix3d Rwb0;
    Rwb0.block(0, 0, 1, 3) = x0.transpose();
    Rwb0.block(1, 0, 1, 3) = y0.transpose();
    Rwb0.block(2, 0, 1, 3) = z0.transpose();
    mag_yaw = Utility::R2ypr(Rwb0).x();
}
