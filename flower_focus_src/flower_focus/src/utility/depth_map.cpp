

#include <cv_bridge/cv_bridge.h>
#include "../parameter/parameters.h"
#include "../utility/visualization.h"
#include <sensor_msgs/PointCloud.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/convolution_3d.h>
#include <pcl/search/kdtree.h>
#include <pcl/common/gaussian.h>
#include <pcl/surface/mls.h>
#include <thread>

cv::Mat map1x, map1y, map2x, map2y;
int image_col;
int image_row;
double focal_length_px;  // 使用相机内参中的焦距
double baseline_m;
cv::Ptr<cv::StereoSGBM> stereo;
double new_fx, new_fy, new_cx, new_cy;
ros::Publisher pub_point_cloud_dense;
std::mutex mutex_feature;


void convertPointCloud2ToPointCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr& pcl_cloud, sensor_msgs::PointCloud& ros_cloud) {
    // 初始化 ROS PointCloud 消息
    // ros_cloud.header = pcl_cloud->header;


    // 填充点数据
    for (size_t i = 0; i < pcl_cloud->points.size(); ++i) {
        const pcl::PointXYZ& point = pcl_cloud->points[i];
        geometry_msgs::Point32 p;
        // std::cout<<point.x<<","<<point.y<<","<<point.z<<std::endl;
        p.x = point.x;
        p.y = point.y;
        p.z = point.z;
        // std::cout<<p.x<<","<<p.y<<","<<p.z<<std::endl;
        ros_cloud.points.push_back(p);
    }


}


void point_filter(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, ros::Publisher& pub_point_cloud_dense, sensor_msgs::PointCloud& point_cloud) {

// // 统计滤波器去除离群值
//     pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
//     sor.setInputCloud(cloud);
//     sor.setMeanK(50);
//     sor.setStddevMulThresh(1.0);
//     pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
//     sor.filter(*cloud_filtered);

// 半径滤波器去除孤立点
    pcl::RadiusOutlierRemoval<pcl::PointXYZ> ror;
    ror.setInputCloud(cloud);
    ror.setRadiusSearch(0.5);
    ror.setMinNeighborsInRadius(20);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered2(new pcl::PointCloud<pcl::PointXYZ>);
    ror.filter(*cloud_filtered2);

    convertPointCloud2ToPointCloud(cloud_filtered2, point_cloud);
    pub_point_cloud_dense.publish(point_cloud);
}

bool readFromYamlFile(const std::string& filename, cv::Mat& camera_matrix, cv::Mat& dist_coeffs) {
    cv::FileStorage fs(filename, cv::FileStorage::READ);

    if (!fs.isOpened())
        return false;

    if (!fs["model_type"].isNone()) {
        std::string sModelType;
        fs["model_type"] >> sModelType;

        if (sModelType.compare("PINHOLE") != 0 && sModelType.compare("PINHOLE_FULL") != 0)
            return false;
    }


    image_col = static_cast<int>(fs["image_width"]);
    image_row = static_cast<int>(fs["image_height"]);

    cv::FileNode n = fs["distortion_parameters"];
    double m_k1 = static_cast<double>(n["k1"]);
    double m_k2 = static_cast<double>(n["k2"]);
    double m_p1 = static_cast<double>(n["p1"]);
    double m_p2 = static_cast<double>(n["p2"]);

    n = fs["projection_parameters"];
    double m_fx = static_cast<double>(n["fx"]);
    double m_fy = static_cast<double>(n["fy"]);
    double m_cx = static_cast<double>(n["cx"]);
    double m_cy = static_cast<double>(n["cy"]);

    camera_matrix = (cv::Mat_<double>(3, 3) << m_fx, 0, m_cx, 0, m_fy, m_cy, 0, 0, 1);
    dist_coeffs = (cv::Mat_<double>(1, 4) << m_k1, m_k2, m_p1, m_p2);

    return true;
}


void stereo_init(string config_file, ros::NodeHandle& n) {

    cv::Mat camera_matrix1, dist_coeffs1, camera_matrix2, dist_coeffs2;

    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    int pn = config_file.find_last_of('/');
    std::string configPath = config_file.substr(0, pn);

    std::string cam0Calib;
    fsSettings["cam0_calib"] >> cam0Calib;
    std::string cam0Path = configPath + "/" + cam0Calib;

    std::string cam1Calib;
    fsSettings["cam1_calib"] >> cam1Calib;
    std::string cam1Path = configPath + "/" + cam1Calib;

    pub_point_cloud_dense = n.advertise<sensor_msgs::PointCloud>("point_cloud_dense", 1000);


    readFromYamlFile(cam0Path, camera_matrix1, dist_coeffs1);
    readFromYamlFile(cam1Path, camera_matrix2, dist_coeffs2);





    cv::Mat T0, T1;
    fsSettings["body_T_cam0"] >> T0;
    fsSettings["body_T_cam1"] >> T1;




    // 计算 T01 = inv(T1) * T0
    cv::Mat T01 = T1.inv() * T0;

    // 提取旋转矩阵 R 和平移向量 T
    cv::Mat R = T01(cv::Range(0, 3), cv::Range(0, 3));
    cv::Mat T = T01(cv::Range(0, 3), cv::Range(3, 4));
    std::cout << "T:\n" << T.t() << std::endl;

    // 立体校正
    cv::Mat R1, R2, P1, P2, Q;
    cv::stereoRectify(camera_matrix1, dist_coeffs1, camera_matrix2, dist_coeffs2,
                      cv::Size(image_col, image_row), R, T,
                      R1, R2, P1, P2, Q);
    focal_length_px = P1.at<double>(0, 0) ;
    new_fx = P1.at<double>(0, 0);
    new_fy = P1.at<double>(1, 1);
    new_cx = P1.at<double>(0, 2);
    new_cy = P1.at<double>(1, 2);

    // 计算校正映射

    cv::initUndistortRectifyMap(camera_matrix1, dist_coeffs1, R1, P1, cv::Size(image_col, image_row), CV_32FC1, map1x, map1y);
    cv::initUndistortRectifyMap(camera_matrix2, dist_coeffs2, R2, P2, cv::Size(image_col, image_row), CV_32FC1, map2x, map2y);
    baseline_m = cv::norm(T);   // 基线距离

    // 初始化 StereoSGBM 对象
    int window_size = 9;
    int min_disp = 0;
    int num_disp = 128;

    stereo = cv::StereoSGBM::create(
                 min_disp, num_disp, window_size,
                 8 * 3 * window_size * window_size,
                 32 * 3 * window_size * window_size,
                 1, 60, 15, 100, 2, cv::StereoSGBM::MODE_SGBM_3WAY
             );




}



queue<sensor_msgs::ImageConstPtr> img_left_buf;
queue<sensor_msgs::ImageConstPtr> img_right_buf;
std::map<double, cv::Mat, less<double>>depths;


cv::Mat getImageFromMsg(const sensor_msgs::ImageConstPtr& img_msg) {
    cv_bridge::CvImageConstPtr ptr;
    if (img_msg->encoding == "8UC1") {
        sensor_msgs::Image img;
        img.header = img_msg->header;
        img.height = img_msg->height;
        img.width = img_msg->width;
        img.is_bigendian = img_msg->is_bigendian;
        img.step = img_msg->step;
        img.data = img_msg->data;
        img.encoding = "mono8";
        ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
    } else
        ptr = cv_bridge::toCvCopy(img_msg, sensor_msgs::image_encodings::MONO8);
    cv::Mat img = ptr->image.clone();
    return img;
}
std::deque<std::vector<pcl::PointXYZ>> depth_map_queue;

void publish_point_cloud(double time_stamp, cv::Mat dense_map, Eigen::Matrix3f R, Eigen::Vector3f t) {




    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(time_stamp);
    sensor_msgs::PointCloud point_cloud;
    point_cloud.header = header;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    std::vector<pcl::PointXYZ> tmp;

    for (int i = 0; i < dense_map.rows; i += 4) {
        for (int j = 0; j < dense_map.cols; j += 4) {
            float dense = dense_map.at<float>(i, j);

            if (dense > 0) {
                Eigen::Vector3f point;
                point << (j - new_cx) / new_fx* dense, (i - new_cy) / new_fy* dense, dense;
                point = R * point + t;
                geometry_msgs::Point32 p;
                p.x = point(0);
                p.y = point(1);
                p.z = point(2);
                tmp.push_back(pcl::PointXYZ(p.x, p.y, p.z));

                // point_cloud.points.push_back(p);




            }
        }
    }

    depth_map_queue.push_back(tmp);
    if (depth_map_queue.size() > 1)depth_map_queue.pop_front();

    for (auto& it : depth_map_queue) {
        for (int i = 0; i < (int)it.size(); i++)
            cloud->push_back(it[i]);
    }


    point_filter(cloud, pub_point_cloud_dense, point_cloud);


}


void process_stereo() {
#if REAL_TIME
    while (1) {
        cv::Mat image0, image1;
        double time = 0;
        {
            mutex_feature.lock();
            if (!img_left_buf.empty() && !img_right_buf.empty()) {
                while (img_left_buf.size() > 5 && img_right_buf.size() > 5) {
                    img_left_buf.pop();
                    img_right_buf.pop();
                    std::cout << "pop stereo\r\n";
                }
                double time0 = img_left_buf.front()->header.stamp.toSec();
                double time1 = img_right_buf.front()->header.stamp.toSec();
                if (fabs(time0 - time1) < 1e-3)
                    time0 = time1;
                if (time0 < time1) {
                    time = img_left_buf.front()->header.stamp.toSec();
                    // image0 = getImageFromMsg(img_left_buf.front());
                    img_left_buf.pop();
                } else if (time0 > time1) {
                    img_right_buf.pop();
                    printf("throw img1\n");
                } else {
                    time = img_left_buf.front()->header.stamp.toSec();
                    image0 = getImageFromMsg(img_left_buf.front());
                    img_left_buf.pop();
                    image1 = getImageFromMsg(img_right_buf.front());
                    img_right_buf.pop();
                }
            }
            mutex_feature.unlock();
        }
        static int count = 0;

        if (!image0.empty() && count++ % 3 == 0) {

            // 读取左、右图像
            cv::Mat left_image = image0;
            cv::Mat right_image = image1;
            TicToc ta;


            // 应用映射进行校正
            cv::Mat left_gray, right_gray;
            assert(!image0.empty());
            cv::remap(left_image, left_gray, map1x, map1y, cv::INTER_LINEAR);
            assert(!image1.empty());
            cv::remap(right_image, right_gray, map2x, map2y, cv::INTER_LINEAR);

            // 显示校正后的图像
            // cv::imshow("Left Map", right_gray);
            // cv::imshow("Right Map", right_gray);
            // std::cout << "1:" << ta.toc() << std::endl;




            // 计算视差图
            cv::Mat disparity_map;
            stereo->compute(left_gray, right_gray, disparity_map);
            // std::cout << "2:" << ta.toc() << std::endl;
            disparity_map.convertTo(disparity_map, CV_32F, 1.0 / 16.0);
            // 获取相机焦距（像素）和基线距离（米）

            // cv::Mat disparity_map_filtered;
            // cv::medianBlur(disparity_map, disparity_map_filtered, 5);
            // 计算深度图
            cv::Mat depth_map = cv::Mat::zeros(disparity_map.size(), CV_32F);
            for (int i = 0; i < disparity_map.rows; ++i) {
                for (int j = 0; j < disparity_map.cols; ++j) {
                    float disparity = disparity_map.at<float>(i, j);
                    if (disparity > 0.1) {
                        double dist = (focal_length_px * baseline_m) / (disparity + 1e-6);
                        if (dist < 20)
                            depth_map.at<float>(i, j) = dist;
                        // if (depth_map.at<float>(i, j) > 20)depth_map.at<float>(i, j) = 20;
                    }
                }
            }
            // std::cout << "3:" << ta.toc() << std::endl;
            depths[time] = depth_map;
#if ENABLE_STEREO_MAP
            double select_time = -1;
            for (auto it = Pubts.begin(), it_next = Pubts.begin(); it != Pubts.end(); it = it_next) {
                it_next++;
                if (depths.find(it->first) != depths.end()) {
#if REAL_TIME
                    mutex_pose.lock();
#endif
                    Eigen::Matrix3f R = PubRs[it->first];
                    Eigen::Vector3f P = Pubts[it->first];
                    PubRs.erase(it->first);
                    Pubts.erase(it->first);
                    double time = it->first;
#if REAL_TIME
                    mutex_pose.unlock();
#endif
                    publish_point_cloud(time, depths[time], R, P);
                    depths.erase(time);
                    select_time = time;
                    break;
                }
            }


            for (auto it = depths.begin(), it_next = depths.begin(); it != depths.end(); it = it_next) {
                it_next++;
                if (it->first < select_time)depths.erase(it);
            }
#if REAL_TIME
            mutex_pose.lock();
#endif
            for (auto it = Pubts.begin(), it_next = Pubts.begin(); it != Pubts.end(); it = it_next) {
                it_next++;
                if (it->first < select_time)Pubts.erase(it);
            }
            for (auto it = PubRs.begin(), it_next = PubRs.begin(); it != PubRs.end(); it = it_next) {
                it_next++;
                if (it->first < select_time)PubRs.erase(it);
            }
#if REAL_TIME
            mutex_pose.unlock();
#endif
            std::cout << depths.size() << "," << Pubts.size() << "," << Pubts.size() << std::endl;

#endif
            // std::cout << "4:" << ta.toc() << std::endl;


            // depth_map.at<float>(0, 0) = 20;
            // std::cout << "3:" << ta.toc() << std::endl;

            // 归一化深度图以便于可视化
            // cv::Mat depth_map_normalized;
            // cv::normalize(depth_map, depth_map_normalized, 0, 255, cv::NORM_MINMAX, CV_8U);
            // std::cout << "4:" << ta.toc() << std::endl;


            // // 可视化深度图
            // cv::Mat depth_map_colored;
            // cv::applyColorMap(depth_map_normalized, depth_map_colored, cv::COLORMAP_BONE);
            // cv::imshow("Depth Map", depth_map_colored);
            // printf("%.3f,%.3f\r\n",time,ta.toc());
            // std::cout << "5:" <<time<<","<< ta.toc() << std::endl;
            // std::chrono::milliseconds dura(1000);
            // std::this_thread::sleep_for(dura);
            // cv::waitKey(1);
        }
    }
    if (img_right_buf.size() && img_left_buf.size()) {
        std::chrono::milliseconds dura(1);
        std::this_thread::sleep_for(dura);
    }
#else
{
        cv::Mat image0, image1;
        double time = 0;
        {
            mutex_feature.lock();
            if (!img_left_buf.empty() && !img_right_buf.empty()) {
                while (img_left_buf.size() > 5 && img_right_buf.size() > 5) {
                    img_left_buf.pop();
                    img_right_buf.pop();
                    std::cout << "pop stereo\r\n";
                }
                double time0 = img_left_buf.front()->header.stamp.toSec();
                double time1 = img_right_buf.front()->header.stamp.toSec();
                if (fabs(time0 - time1) < 1e-3)
                    time0 = time1;
                if (time0 < time1) {
                    time = img_left_buf.front()->header.stamp.toSec();
                    // image0 = getImageFromMsg(img_left_buf.front());
                    img_left_buf.pop();
                } else if (time0 > time1) {
                    img_right_buf.pop();
                    printf("throw img1\n");
                } else {
                    time = img_left_buf.front()->header.stamp.toSec();
                    image0 = getImageFromMsg(img_left_buf.front());
                    img_left_buf.pop();
                    image1 = getImageFromMsg(img_right_buf.front());
                    img_right_buf.pop();
                }
            }
            mutex_feature.unlock();
        }
        static int count = 0;

        if (!image0.empty() && count++ % 3 == 0) {

            // 读取左、右图像
            cv::Mat left_image = image0;
            cv::Mat right_image = image1;
            TicToc ta;


            // 应用映射进行校正
            cv::Mat left_gray, right_gray;
            assert(!image0.empty());
            cv::remap(left_image, left_gray, map1x, map1y, cv::INTER_LINEAR);
            assert(!image1.empty());
            cv::remap(right_image, right_gray, map2x, map2y, cv::INTER_LINEAR);

            // 显示校正后的图像
            // cv::imshow("Left Map", right_gray);
            // cv::imshow("Right Map", right_gray);
            // std::cout << "1:" << ta.toc() << std::endl;




            // 计算视差图
            cv::Mat disparity_map;
            stereo->compute(left_gray, right_gray, disparity_map);
            // std::cout << "2:" << ta.toc() << std::endl;
            disparity_map.convertTo(disparity_map, CV_32F, 1.0 / 16.0);
            // 获取相机焦距（像素）和基线距离（米）

            // cv::Mat disparity_map_filtered;
            // cv::medianBlur(disparity_map, disparity_map_filtered, 5);
            // 计算深度图
            cv::Mat depth_map = cv::Mat::zeros(disparity_map.size(), CV_32F);
            for (int i = 0; i < disparity_map.rows; ++i) {
                for (int j = 0; j < disparity_map.cols; ++j) {
                    float disparity = disparity_map.at<float>(i, j);
                    if (disparity > 0.1) {
                        double dist = (focal_length_px * baseline_m) / (disparity + 1e-6);
                        if (dist < 20)
                            depth_map.at<float>(i, j) = dist;
                        // if (depth_map.at<float>(i, j) > 20)depth_map.at<float>(i, j) = 20;
                    }
                }
            }
            // std::cout << "3:" << ta.toc() << std::endl;
            depths[time] = depth_map;
#if ENABLE_STEREO_MAP
            double select_time = -1;
            for (auto it = Pubts.begin(), it_next = Pubts.begin(); it != Pubts.end(); it = it_next) {
                it_next++;
                if (depths.find(it->first) != depths.end()) {
#if REAL_TIME
                    mutex_pose.lock();
#endif
                    Eigen::Matrix3f R = PubRs[it->first];
                    Eigen::Vector3f P = Pubts[it->first];
                    PubRs.erase(it->first);
                    Pubts.erase(it->first);
                    double time = it->first;
#if REAL_TIME
                    mutex_pose.unlock();
#endif
                    publish_point_cloud(time, depths[time], R, P);
                    depths.erase(time);
                    select_time = time;
                    break;
                }
            }


            for (auto it = depths.begin(), it_next = depths.begin(); it != depths.end(); it = it_next) {
                it_next++;
                if (it->first < select_time)depths.erase(it);
            }
#if REAL_TIME
            mutex_pose.lock();
#endif
            for (auto it = Pubts.begin(), it_next = Pubts.begin(); it != Pubts.end(); it = it_next) {
                it_next++;
                if (it->first < select_time)Pubts.erase(it);
            }
            for (auto it = PubRs.begin(), it_next = PubRs.begin(); it != PubRs.end(); it = it_next) {
                it_next++;
                if (it->first < select_time)PubRs.erase(it);
            }
#if REAL_TIME
            mutex_pose.unlock();
#endif
            std::cout << depths.size() << "," << Pubts.size() << "," << Pubts.size() << std::endl;

#endif
            // std::cout << "4:" << ta.toc() << std::endl;


            // depth_map.at<float>(0, 0) = 20;
            // std::cout << "3:" << ta.toc() << std::endl;

            // 归一化深度图以便于可视化
            // cv::Mat depth_map_normalized;
            // cv::normalize(depth_map, depth_map_normalized, 0, 255, cv::NORM_MINMAX, CV_8U);
            // std::cout << "4:" << ta.toc() << std::endl;


            // // 可视化深度图
            // cv::Mat depth_map_colored;
            // cv::applyColorMap(depth_map_normalized, depth_map_colored, cv::COLORMAP_BONE);
            // cv::imshow("Depth Map", depth_map_colored);
            // printf("%.3f,%.3f\r\n",time,ta.toc());
            // std::cout << "5:" <<time<<","<< ta.toc() << std::endl;
            // std::chrono::milliseconds dura(1000);
            // std::this_thread::sleep_for(dura);
            // cv::waitKey(1);
        }
    }
#endif



}

void stereo_img1_callback(const sensor_msgs::ImageConstPtr& img_msg) {
    mutex_feature.lock();
    img_right_buf.push(img_msg);
    mutex_feature.unlock();
}

void stereo_img0_callback(const sensor_msgs::ImageConstPtr& img_msg) {
    mutex_feature.lock();
    img_left_buf.push(img_msg);
    mutex_feature.unlock();

}



