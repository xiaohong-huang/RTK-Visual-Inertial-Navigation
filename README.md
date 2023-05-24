# RTK-Visual-Inertial-Navigation

A Sliding Window Filter with GNSS-State Constraint for RTK-Visual-Inertial Navigation. [paper link]()

Authors: Xiaohong Huang, Cui Yang

**RTK-Visual-Inertial-Navigation** is a navigation system that tightly fuses GNSS, visual, and inertial measurements. It uses a sliding window filter (SWF) with GNSS-state constraints for sensor fusion. That is, the GNSS states (i.e., position, orientation, and velocity of the body and inertial biases at the time of capturing GNSS measurements) are retained in the SWF to construct more appropriate constraints between measurements and states. It also uses the parallel elimination strategy in a predefined elimination ordering, which can solve the Gauss-Newton problem and simultaneously obtain the covariance for ambiguity resolution. The system can perform the following types of navigation:

- RTK-Visual-Inertial Navigation;
- RTD-Visual-Inertial Navigation;
- SPP-Visual-Inertial Navigation;
- SPP-Visual-Inertial Navigation with Carrier-Phase Fusion
- Visual-Inertial navigation.


## 1. Prerequisites
### 1.1 C++11 Compiler
This package requires some features of C++11.

### 1.2 ROS
This package is developed under [ROS Kinetic](http://wiki.ros.org/kinetic) environment.

### 1.3 Opencv 3
Our code uses [Opencv 3](https://github.com/opencv/opencv/tree/3.4) for image process.

## 2. Build RTK-Visual-Inertial-Navigation
Clone the repository to your catkin workspace (for example `~/catkin_ws/`):
```
cd ~/catkin_ws/src/
git clone https://github.com/xiaohong-huang/RTK-Visual-Inertial-Navigation.git
```
In our source code, we have developed our solving strategy based on [Ceres-Solver](http://ceres-solver.org/). The original version of [Ceres-Solver](http://ceres-solver.org/) is not satisfied for our project. To build the project, you need to build our modified Ceres-Solver with:
```
# CMake
sudo apt-get install cmake
# Eigen3
sudo apt-get install libeigen3-dev
# Ceres-Solver-Modified
cd ~/catkin_ws/src/RTK-Visual-Inertial-Navigation
tar -xvf ceres-solver-modified.tar
cd ceres-solver-modified/
sh build.sh
```
The modified version will only be installed in the current folder. So you don't need to worry that the installation will change the settings of your computer.

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
roslaunch rtk_visual_inertial rtk_visual_inertial_rviz.launch
```
Open another terminal and run the project by:
```
source ~/catkin_ws/devel/setup.bash
rosrun rtk_visual_inertial rtk_visual_inertial_node src/RTK-Visual-Inertial-Navigation/yaml/SETTING.yaml YOUR_BAG_FOLDER/BAG_NAME.bag ourput.csv
```
YOUR_BAG_FOLDER is the folder where you save our dataset. BAG_NAME is the name of our dataset. SETTING.yaml is the setting for RTK-Visual-Inertial-Navigation. You could use the following settings to perform different types of navigation.
```
rtk_visual_inertial_config.yaml     #RTK-Visual-Inertial-Navigation
rtd_visual_inertial_config.yaml     #RTD-Visual-Inertial-Navigation
spp_visual_inertial_config.yaml     #SPP-Visual-Inertial-Navigation
spp_CP_visual_inertial_config.yaml  #SPP-Visual-Inertial-Navigation with carrier-phase fusion
visual_inertial_config.yaml         #Visual-Inertial-Navigation
```
We have also provide a demo for evaluating the positioning errors (see [evaluate.py](https://github.com/xiaohong-huang/RTK-Visual-Inertial-Navigation/blob/main/evaluate/evaluate.py)). 
## 4. Acknowledgements
The VIO framework is adapted from [VINS-Mono](https://github.com/HKUST-Aerial-Robotics/VINS-Mono). The Ceres-Solver-Modified is developed base on [Ceres-Solver](http://ceres-solver.org/)
## 5. Licence
The source code is released under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license.
