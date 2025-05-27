# FLOWER-FOCUS




## 1. Prerequisites
### 1.1 C++14 Compiler
This package requires some features of C++14.

### 1.2 ROS
This package is developed and tested under [ROS Melodic (Ubuntu 18.04)](http://wiki.ros.org/melodic) or [ROS Noetic (Ubuntu 20.04)](http://wiki.ros.org/noetic) environment.
### 1.3 OCTOMAP
```
#ROS Noetic
sudo apt install ros-noetic-octomap-ros
sudo apt install ros-noetic-octomap-rviz-plugins
sudo apt-get install ros-noetic-pcl-ros

#ROS Melodic
#sudo apt install ros-melodic-octomap-ros
#sudo apt install ros-melodic-octomap-rviz-plugins
#sudo apt-get install ros-melodic-pcl-ros
```

### 1.4 PCL
```
sudo apt-get install libpcl-dev
```

## 2. Build FLOWER-FOCUS
Clone the repository to your catkin workspace (for example `~/catkin_ws/`):
```
cd ~/catkin_ws/src/
git clone https://github.com/xiaohong-huang/RTK-Visual-Inertial-Navigation
git checkout FLOWER-FOCUS
```
Build the OpenCV4 (>=4.3.0):
```
#Clone the Opencv to the folder
cd ~/catkin_ws/src/RTK-Visual-Inertial-Navigation
git clone https://github.com/opencv/opencv.git
cd opencv

#Switch the blanch to OpenCV 4.10.0. 
git checkout 4.10.0

#build opencv
mkdir build
cd build
cmake ..
make -j8
#Do not use "make install"

#build ceres-solver-modified
cd ~/catkin_ws/src/RTK-Visual-Inertial-Navigation
tar -xvf ceres-solver-modified.tar
cd ceres-solver-modified/
sh build.sh
```
Note, the OpenCV version must be larger than 4.3.0 (the newest version is 4.10.0, which is work well in our project). Otherwise, some of the function may not work well.

Then build the package with:
```
cd ~/catkin_ws/
catkin_make
```



## 3. Run RTK-Visual-Inertial-Navigation with our dataset
Our equipment is shown as follows: A grayscale camera (MT9V034 752x480@25HZ), a MEMS-grade IMU (BMI088 400HZ), a $360^o$ prism, and a GNSS receiver (ublox ZED-F9P 10HZ) are installed together with a small GNSS antenna (BT-560). A Trimble S9 total station is installed in a fixed position and observes the prism to generate the ground truth of the rover station every 0.1 seconds with mm-level accuracy. A GNSS receiver (ublox ZED-F9P 1HZ) with an experimental-level antenna is installed in a fixed position for the base station.
![image](https://github.com/xiaohong-huang/RTK-Visual-Inertial-Navigation/blob/main/fig/equipment.png)
The experiment environment is shown as follows.
![image](https://github.com/xiaohong-huang/RTK-Visual-Inertial-Navigation/blob/main/fig/experiment_sense.png)

Download our [Dataset](https://1drv.ms/f/s!ApdCy_pJvU0qyVsLB906CNjAEQiH) and launch the rviz via:
```
source ~/catkin_ws/devel/setup.bash
roslaunch flower_focus visual_inertial_rviz.launch
```
Open another terminal and run the project by:
```
source ~/catkin_ws/devel/setup.bash
rosrun flower_focus flower_focus_node src/RTK-Visual-Inertial-Navigation/yaml/SETTING.yaml YOUR_BAG_FOLDER/BAG_NAME.bag ourput.csv
```
YOUR_BAG_FOLDER is the folder where you save our dataset. BAG_NAME is the name of our dataset. SETTING.yaml is the setting for RTK-Visual-Inertial-Navigation. You could use the following settings to perform different types of navigation.
```
rtk_visual_inertial_config.yaml     #RTK-Visual-Inertial-Navigation
rtd_visual_inertial_config.yaml     #RTD-Visual-Inertial-Navigation
spp_visual_inertial_config.yaml     #SPP-Visual-Inertial-Navigation
spp_CP_visual_inertial_config.yaml  #SPP-Visual-Inertial-Navigation with carrier-phase fusion
visual_inertial_config.yaml         #Visual-Inertial-Navigation
```



## 4. Acknowledgements
The VIO framework is adapted from [FLOWER-VIO](https://github.com/xiaohong-huang/FLOWER-VIO). The Ceres-Solver-Modified is developed base on [Ceres-Solver](http://ceres-solver.org/).
## 5. Licence
The source code is released under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license.


