import pandas as pd
import matplotlib.pyplot as plt
from joypy import joyplot
import os
from matplotlib import cm
import numpy as np

# Assuming you have the absolute error values stored in a 2D NumPy array called 'error_grid'
# with shape (21, 24)

# Generate some random data for demonstration
start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\inputs\neural_network_par\data\ipopt\individual_pipes"

# all_data = pd.DataFrame()
pipe_knm = np.array([255.5169343, 255.5169343, 208.6287032,208.6287032,100.2219872,26.8636084,32.71143353,40.41306369,
                     68.90779278,114.2706469,13.94307458,102.2067737,12.47106503,78.85423786,80.8015493,
                    228.5412938,161.6030986,102.2067737,19.24327111,3.501408717,14.15077486])




# for file in os.listdir(start_directory):
#     data_file_path = os.path.join(start_directory, file)
#     print(data_file_path)
#     # Check if the item is a subfolder
#     # if os.path.isdir(folder_path):
#     # excel_file_path = os.path.join(folder_path, foldername +"_results"+ '.xlsx')
        
#     # Check if the corresponding Excel file exists
#         # if os.path.isfile(excel_file_path):
#     # Load the Excel file and read the 'feasibility_gap' sheet
#     data_df = pd.read_csv(data_file_path, names = ["q", "p1", "p2", "k"])


#     # all_data[file[30:-4]] = data_df["p1"]**2-data_df["p2"]**2
#     all_data[file[30:-4]] = data_df["q"]


# cs= ['pipe1', 'pipe2','pipe3', 'pipe4', 'pipe5', 'pipe6', 'pipe7', 'pipe8', 'pipe9',  'pipe10', 'pipe11', 'pipe12', 'pipe13', 'pipe14', 'pipe15',
#        'pipe16', 'pipe17', 'pipe18', 'pipe19', 'pipe20', 'pipe21']
# all_data[cs]

# all_data_s = all_data[cs]





save_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\feasible_operating_data"

# all_data.isna().sum()
# #Save the dataframe
# all_data_s.to_csv(os.path.join(save_path, "all_data_q.csv"), index=False)


all_data_k = pd.read_csv(os.path.join(save_path, "all_data_q.csv"))
all_data = all_data_k/(pipe_knm) #uuuups should have divided by knm not knm**2...

def ridge_plot(columns,   name, save_path= save_path,fs=(4,3),data=all_data, cmap=plt.cm.get_cmap('viridis_r'),alpha=0.5, color="green"):
    fig, axes = joyplot(
        data=data[columns],
        #  data=all_data["pipe6"],
        figsize=fs, 
        grid='y',
        # ylim='own',
        # overlap=1,
        linewidth=1,
        alpha=alpha,
        range_style='own',
        colormap=cmap
        # color= color,
    )
    plt.xlabel("mass flow [q/K]",fontsize=16)
    plt.ylabel("pipes",fontsize=16)
    for ax in axes:
        # label = ax.get_yticklabels()
        # ax.set_yticklabels(label, fontdict={'color': 'r'})
        ax.tick_params(axis='y', labelsize=14)
    
    # ax.tick_params(axis='y', labelsize=14) 
    # ax.tick_params(axis='x', labelsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(os.path.join(save_path, name),bbox_inches="tight",)

    return


cs= ['pipe1', 'pipe2','pipe3', 'pipe4', 'pipe5', 'pipe6', 'pipe7', 'pipe8', 'pipe9',  'pipe10', 'pipe11', 'pipe12', 'pipe13', 'pipe14', 'pipe15',
       'pipe16', 'pipe17', 'pipe18', 'pipe19', "pipe20", 'pipe21']
# cs_wide = ["pipe20","pipe6", "pipe7", "pipe8", "pipe11", "pipe13", "pipe19", "pipe21"]

# all_data[["pipe20"]]= all_data[["pipe20"]]/40
# all_data[["pipe10"]]= all_data[["pipe10"]]/3
# wide; 6,7,8,11,13,19,21

# all_data[["pipe8"]] <-0.4

# all_data["pipe8"] = all_data[all_data[["pipe8"]] <-0.4]["pipe8"]


# 1,2,3,4,16, 10, 11
cs_narrow = ["pipe1", "pipe2", "pipe3", "pipe4", "pipe10", "pipe11", "pipe16"]




ridge_plot(columns=cs, name="all_test_norm4.png",fs=(12,20))

# cs_wide = ["pipe6", "pipe7", "pipe8", "pipe11", "pipe13", "pipe19", "pipe21"]
cs_wide = ["pipe6", "pipe7", "pipe8", "pipe12", "pipe13", "pipe19", "pipe21"]
ridge_plot(columns=cs_wide, name="ridge_plot_wide_fix2.png")


# cs_mid = ["pipe1", "pipe2", "pipe3", "pipe4","pipe10","pipe16"]

# cs_mid = ["pipe1", "pipe2", "pipe3", "pipe4","pipe10","pipe16", "pipe20"]
cs_mid = ["pipe1", "pipe2", "pipe3", "pipe4","pipe10","pipe11","pipe16"]
ridge_plot(columns=cs_mid, name="ridge_plot_mid_fix2.png")


# cs_narrow = [ "pipe5", "pipe9", "pipe12", "pipe14", "pipe15", "pipe17", "pipe18"]
cs_narrow = [ "pipe5", "pipe9",  "pipe14", "pipe15", "pipe17", "pipe18", "pipe20"]
ridge_plot(columns=cs_narrow, name="ridge_plot_nar_fix2.png")

cs_20 = ["pipe20"]
ridge_plot(columns=cs_20, name="ridge_plot_20.png")



#plot 7 and 20
cs_7 = ["pipe7"]
cs_20 = ["pipe20"]
ridge_plot(columns=cs_7, name="ridge_plot_7.png", color="lightgreen", alpha  = 0.45)
ridge_plot(columns=cs_20, name="ridge_plot_20.png", color = "green",  alpha = 0.45)
# ridge_plot(columns=cs_narrow, name="ridge_plot_nar_fix2.png")

#yellow #9aadc5
#green #a4dfb5


#=========uniform data==================

save_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\Uniform_results"
uniform_data = pd.read_csv(os.path.join(save_path, "training_data_uniform_num_samples_100k_seed_42.csv"))

ridge_plot(["q"],   name = "uniform_den2", save_path= save_path,fs=(4,2),data=uniform_data ,alpha=0.5, color= "blue")



save_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\Linspace_results"
linspace_data = pd.read_csv(os.path.join(save_path, "training_data_linspace_num_samples_1000k_seed_42.csv"))

ridge_plot(["q"],   name = "uniform_den2", save_path= save_path,fs=(4,2),data=linspace_data ,alpha=0.5)

# ridge_plot(columns,   name, save_path= save_path,fs=(4,3),data=all_data, cmap=plt.cm.get_cmap('viridis_r'),alpha=0.5)



# uniform_data.hist(column=["q"], bins=100)
# plt.savefig(os.path.join(save_path, "uniform_hist_q.png"),bbox_inches="tight",)


# linspace_data.hist(column=["q"], bins=100)
# plt.savefig(os.path.join(save_path, "linspace_hist_q.png"),bbox_inches="tight",)


# uniform_data["p1"]**2-uniform_data["p2"]**2


# uniform_data["psq"] = uniform_data["p1"]**2-uniform_data["p2"]**2


# joyplot(
#     data=uniform_data[["q"]],
#     #  data=all_data["pipe6"],
#     figsize=fs, 
#     grid='y',
#     # ylim='own',
#     range_style='own',
#     overlap=0,
#     linewidth=1,
#     alpha=0.6,
#     hist=True,
#     bins=100,
#     # kind = 'lognorm'
#     # ymax = 10,
#     # colormap=cm.autumn_r,
# )
# plt.savefig(os.path.join(save_path, "uniform_hist_q.png"),bbox_inches="tight",)

# np.sqrt(2.228800)


uniform_data[uniform_data[["q"]]>0].mean()


