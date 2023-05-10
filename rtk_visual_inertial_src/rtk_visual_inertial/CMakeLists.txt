cmake_minimum_required(VERSION 2.8.3)
project(rtk_visual_inertial)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    geometry_msgs
    nav_msgs
    tf
    cv_bridge
    camera_models
    image_transport
    rosbag)
    
    
find_package(OpenCV 3 REQUIRED)
include(../../ceres-solver-modified/build/lib/cmake/Ceres/CeresConfig.cmake)
find_package(Eigen3)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${OpenCV_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  src/gnss/include
  src/swf/
  src/factor/
  src/utility/
  src/initial/
  src/feature/
)

catkin_package()
message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")
message(STATUS "    catkin libraries: ${catkin_LIBRARIES}")

add_library(rtk_visual_inertial_lib
    src/parameter/parameters.cpp
    src/swf/swf.cpp
    src/swf/swf_imu.cpp
    src/swf/swf_image.cpp
    src/swf/swf_gnss.cpp
    src/swf/swf_core.cpp
    src/swf/swf_lambda.cpp
    src/feature/feature_manager.cpp
    src/feature/feature_tracker.cpp
    src/factor/gnss_imu_factor.cpp    
    src/factor/pose_local_parameterization.cpp
    src/factor/projection_factor.cpp
    src/factor/marginalization_factor.cpp
    src/factor/mag_factor.cpp
    src/factor/gnss_factor.cpp
    src/factor/imu_factor.cpp
    src/factor/integration_base.cpp
    src/factor/initial_factor.cpp
    src/factor/pose0_factor.cpp
    src/utility/visualization.cpp
    src/utility/camera_pose_visualization.cpp
    src/gnss/src/lambda.cpp
    src/gnss/src/common_function.cpp
    )
    
target_link_libraries(rtk_visual_inertial_lib ${OpenCV_LIBS} ${catkin_LIBRARIES}  ${CERES_LIBRARIES})


add_executable(rtk_visual_inertial_node src/main3.cpp)
target_link_libraries(rtk_visual_inertial_node rtk_visual_inertial_lib) 




