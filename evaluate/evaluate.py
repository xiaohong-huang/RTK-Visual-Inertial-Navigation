import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.spatial.transform import Rotation as s_R
#step:
#1. saving the results and ground files in the current folder.
#2. renaming the result file names according to the bag names. i.e. if the R1M1.bag is used, renaming the  
# result file name as R1M1_rtk.csv (for RTK positioning) or R1M1_spp.csv (for SPP positioning).
#3. run python evaluate.py



def plot_error(ground_filename, data_filename, align=False, label="",
               t_name="", is_rtk=False):
    try:
        data = pd.read_csv(data_filename)
    except Exception as ret:
        print(ret)
        return

    data[np.isnan(data)] = 1
    # calculating the prism position according to the antenna position and the body orientation.
    data[["px", "py", "pz"]] -= s_R.from_euler("zyx", data[["yaw", "pitch", "roll"]].values,
                                               degrees=True).as_matrix() @ ptg
    data_time = data["time"].values / 1e9
    data_values = data[["px", "py", "pz"]].values
    ground = pd.read_csv(ground_filename).iloc[10:]

    errors = []
    ground_values = ground[["px", "py", "pz"]].values
    diff_p = ground_values[1:] - ground_values[:-1]
    travel_distance = np.sqrt(diff_p[:, 0] ** 2 + diff_p[:, 1] ** 2 + diff_p[:, 2] ** 2)
    travel_distance = np.concatenate([np.array([0]), travel_distance])
    travel_distance = np.cumsum(travel_distance)
    travel_distance_error = []
    for i in range(len(ground)):
        g = ground[["px", "py", "pz"]].iloc[i].values
        ground_time = ground["time"].iloc[i]

        index = np.abs(ground_time - data_time).argmin()
        if np.abs(ground_time - data_time).min() > 1/400:
            print(ground_time, "not found")
            continue

        tmp = data_values[index] - g
        errors.append(tmp)
        travel_distance_error.append(travel_distance[i])
    errors = np.array(errors)
    if align:
        errors -= errors.mean(axis=0)

    plt.plot(np.array(travel_distance_error), errors[:, 0], label="x")
    plt.plot(np.array(travel_distance_error), errors[:, 1], label="y")
    plt.plot(np.array(travel_distance_error), errors[:, 2], label="z")
    plt.show()

    if not is_rtk:
        global RMSE_PLANE,RMSE_HEIGHT
        RMSE_HEIGHT[t_name][label] = np.sqrt((errors[:, 2] ** 2).mean())
        RMSE_PLANE[t_name][label] = np.sqrt((errors[:, 0] ** 2 + errors[:, 1] ** 2).mean())

    else:
        global MAE
        MAE[t_name][label] = np.sqrt((errors[:, 0] ** 2 + errors[:, 1] ** 2 + errors[:, 2] ** 2)).mean()



#imu-prism calibration
ptg = np.array([0.04128228786, -0.02040929358, -0.1396607903])

files = ["R1M1", "R1M2", "R2M1", "R2M2"]

RMSE_HEIGHT = {"R1M1":{},"R1M2":{},"R2M1":{},"R2M2":{}}
RMSE_PLANE = {"R1M1":{},"R1M2":{},"R2M1":{},"R2M2":{}}
for f in files:

    plot_error( f + "_ground.csv",
               f + "_spp.csv",
                align=True, label="Proposed",
               t_name=f)
print("                         plane RMSE:")
print(pd.DataFrame(RMSE_PLANE))
print("                         height RMSE:")
print(pd.DataFrame(RMSE_HEIGHT))

MAE={"R1M1":{},"R1M2":{},"R2M1":{},"R2M2":{}}
for f in files:
    plot_error( f + "_ground.csv",
               f + "_rtk.csv",
               is_rtk=True,
               t_name=f,
               label='Proposed')

print("                         MAE:")
print(pd.DataFrame(MAE))






