# C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\models\individual_pipes_run_2023_05_07_15_05_09

import json
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#BOUND TIGHTENING UNIFORM VS LINSPACE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Bound tightening"

start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\models\individual_pipes_run_2023_05_07_15_05_09"

data_dict = {}
MSE_list = []
# Iterate over subfolders
for file in os.listdir(start_directory):
    file_path = os.path.join(start_directory, file)

    f = open(file_path)
    data = json.load(f)
    f.close()

    data_dict[file] = data
    if "MSE" in data.keys():
        # print(file,data["MSE"])
        
        MSE_list.append(data["MSE"])
    # print(file)
    # print(f'max y_hat {np.max(data["y1_hat"])}, {np.max(data["y2_hat"])}')
    # print(f'min y_hat {np.min(data["y1_hat"])}, {np.min(data["y2_hat"])}\n')

    if "pipe7" in file:
        # data["y1_hat"] = np.array(data["y1_hat"])
        # data["y2_hat"] = np.array(data["y2_hat"])
        print(file)
        print(data["y1_hat"])
        print(data["y2_hat"])
    if "pipe20" in file:
        print(file)
        print(data["y1_hat"])
        print(data["y2_hat"])


    # folder_path = os.path.join(start_directory, foldername)
    
np.array(MSE_list).mean()