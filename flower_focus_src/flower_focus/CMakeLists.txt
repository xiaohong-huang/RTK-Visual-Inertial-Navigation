cmake_minimum_required(VERSION 2.8.3)
project(flower_focus)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(cv_bridge_DIR /usr/local/share/cv_bridge/cmake)
find_package(PCL REQUIRED)
find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    geometry_msgs
    nav_msgs
    tf
    cv_bridge
    camera_models
    image_transport
    rosbag
    PCL
    )
    


# find_package(OpenCV 4 REQUIRED)
#include(../../../../../opencv_2/build/OpenCVConfig.cmake)
include(../../opencv/build/OpenCVConfig.cmake)
include(../../ceres-solver-modified/build/lib/cmake/Ceres/CeresConfig.cmake)
find_package(Eigen3)
find_package ( OpenMP )


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${OpenMP_EXE_LINKER_FLAGS}")

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${OpenCV_INCLUDE_DIRS}
  src/swf/
  src/factor/
  src/utility/
  src/initial/
  src/feature/
  src/gnss/
  ${CERES_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
)


catkin_package()
message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")
message(STATUS "    catkin libraries: ${catkin_LIBRARIES}")

add_library(flower_focus_lib
    src/parameter/parameters.cpp
    src/swf/swf.cpp
    src/swf/swf_imu.cpp
    src/swf/swf_image.cpp
    src/swf/swf_core.cpp
    src/feature/feature_manager.cpp
    src/feature/feature_tracker.cpp
    src/factor/marginalization_factor.cpp
    src/factor/imu_factor.cpp
    src/factor/integration_base.cpp
    src/factor/initial_factor.cpp
    src/factor/visual_inerial_factor.cpp
    src/factor/loss_function.cpp
    src/factor/visual_inerial_base.cpp
    src/factor/visual_inerial_base_trans.cpp
    src/factor/visual_inerial_base_rhs.cpp
    src/factor/gnss_imu_info.cpp
    src/solver/solver.cpp
    src/solver/solver_residualblock.cpp
    src/utility/visualization.cpp
    src/utility/camera_pose_visualization.cpp
    src/utility/depth_map.cpp
    
    src/gnss/lambda.cpp
    src/gnss/common_function.cpp
    src/gnss/gnss_ceres_factor.cpp
    src/gnss/gnss_factor.cpp
    src/gnss/gnss_imu_factor.cpp
    src/swf/swf_gnss.cpp

    )
add_definitions(${PCL_DEFINITIONS})
target_link_libraries(flower_focus_lib ${OpenCV_LIBS} ${catkin_LIBRARIES} ${CERES_LIBRARIES} ${PCL_LIBRARY_DIRS})


add_executable(flower_focus_node src/main3.cpp)
target_link_libraries(flower_focus_node flower_focus_lib) 

add_executable(4seasons_convert2bag_node src/4seasons_convert2bag.cpp)
target_link_libraries(4seasons_convert2bag_node flower_focus_lib) 

# add_executable(depth_map_node src/depth_map.cpp)
# target_link_libraries(depth_map_node flower_focus_lib) 