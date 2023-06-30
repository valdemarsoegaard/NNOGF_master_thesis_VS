import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
import copy
import seaborn as sns
# Specify the directory to start the iteration from
# start_directory = 'path/to/start/directory'

def weymouth_equation(ps):
    q = []
    for i,x in enumerate(ps):
            if x>=0:
                q.append(np.sqrt(ps[i]))
            else:     
                q.append(-np.sqrt(-ps[i]))  
    return  np.array(q)


#BENCHMARK
start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\weymout_opt_comp"

#FEASIBLE TRAINNG DATA/EACH PIPE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\each_pipe_nb"

#BOUND TIGHTENING UNIFORM VS LINSPACE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Bound tightening"

pipe_knm = np.array([255.5169343, 255.5169343, 208.6287032,208.6287032,100.2219872,26.8636084,32.71143353,40.41306369,
                     68.90779278,114.2706469,13.94307458,102.2067737,12.47106503,78.85423786,80.8015493,
                    228.5412938,161.6030986,102.2067737,19.24327111,3.501408717,14.15077486])

p_map= [[1,  1  ,2],
        [2,  1  ,2],
        [3,  2  ,3],
        [4,  2  ,3],
        [5,  3  ,4],
        [6,  5  ,6],
        [7,  6  ,7],
        [8,  7  ,4],
        [9,  4  ,14],
        [10, 9  ,10],
        [11, 9  ,10],
        [12, 10, 11],
        [13, 10, 11],
        [14, 11, 12],
        [15, 12, 13],
        [16, 13, 14],
        [17, 14, 15],
        [18, 15, 16],
        [19, 11, 17],
        [20, 18, 19],
        [21, 19, 20]]
# excel_theo_wey = pd.DataFrame()
# for p in p_map:
#     excel_theo_wey[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)*pipe_knm[p[0]-1]**2

q_p1_x = pd.DataFrame()
q_p1 = pd.DataFrame()


q_p7_x = pd.DataFrame()
q_p7 = pd.DataFrame()
# q_p7 = pd.DataFrame()
q_p20_x = pd.DataFrame()
q_p20 = pd.DataFrame()

# Iterate over subfolders
for foldername in os.listdir(start_directory):
    folder_path = os.path.join(start_directory, foldername)
    
    # Check if the item is a subfolder
    if os.path.isdir(folder_path):
        excel_file_path = os.path.join(folder_path, foldername +"_results"+ '.xlsx')
        
        # Check if the corresponding Excel file exists
        if os.path.isfile(excel_file_path):
            # Load the Excel file and read the 'feasibility_gap' sheet
            print(f"Reading {excel_file_path}")
            excel_qnm = pd.read_excel(excel_file_path, sheet_name='Q_nm')
            excel_parameters = pd.read_excel(excel_file_path, sheet_name='parameters')
            excel_pressure = pd.read_excel(excel_file_path, sheet_name='pr')

            excel_sqp = pd.DataFrame()
            for p in p_map:
                excel_sqp[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)

            p_theo = excel_sqp>=0
            n_theo = excel_sqp<0
            theo_qnm_values_pos = np.sqrt(excel_sqp[p_theo])
            theo_qnm_values_neg = -np.sqrt(-excel_sqp[n_theo])
            theo_qnm_values = np.transpose(theo_qnm_values_pos.fillna(theo_qnm_values_neg).values)

            sqp_values = np.transpose(excel_sqp.values)
            exp_qnm_values = np.transpose(excel_qnm.values/pipe_knm)


            q_p1_x[f"{foldername}"] = sqp_values[0]
            q_p1[f"{foldername}"] = exp_qnm_values[0]




            q_p7_x[f"{foldername}"] = sqp_values[6]
            # q_p7[f"{foldername}_x"] = sqp_values[6]
            q_p7[f"{foldername}"] = exp_qnm_values[6]
            # q_p7[f"{foldername}_q_theo"] = theo_qnm_values[6]
            



            q_p20_x[f"{foldername}"] = sqp_values[19]
            # q_p20[f"{foldername}_x"] = sqp_values[19]
            q_p20[f"{foldername}"] = exp_qnm_values[19]
            # q_p20[f"{foldername}_q_theo"] = theo_qnm_values[19]

q_p1_x.columns = ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
q_p1.columns = ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
columns = ["PLA (BP6)","PLA (BP10)","(U)NNOGF B6","(FOP)NNOGF"]
q_p1 = q_p1[columns]
q_p1_x = q_p1_x[columns]



q_p7_x.columns = ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
q_p7.columns = ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
columns = ["PLA (BP6)","PLA (BP10)","(U)NNOGF B6","(FOP)NNOGF"]
q_p7 = q_p7[columns]
q_p7_x = q_p7_x[columns]

q_p20_x.columns = ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
q_p20.columns = ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
columns = ["PLA (BP6)","PLA (BP10)","(U)NNOGF B6","(FOP)NNOGF"]
q_p20 = q_p20[columns]
q_p20_x = q_p20_x[columns]

size = 60
lw =2
def opt_theo_plot(df_x, df_y,pipe, fs = 18):
     
    # sorted_indices = np.argsort(sqp_values[i])
    # x_sorted = sqp_values[i][sorted_indices]
    # y_theo_sorted = theo_qnm_values[i][sorted_indices]
    # y_exp_sorted = exp_qnm_values[i][sorted_indices]
    if df_x.min().min() < 0:
        mn = df_x.min().min()*1.05
    else:
        mn = df_x.min().min()*0.95
    
    if df_x.max().max() < 0:
        mx = df_x.max().max()*0.95
    else:
        mx = df_x.max().max()*1.05

    x = np.linspace(mn, mx, 100)
    p = weymouth_equation(x)

    sns.set(style='ticks')
    # Create the plot
    plt.figure(figsize=(8, 6))  # Set the figure size

    # Plot the analytical method with sorted data
    sns.lineplot(x = x, y=  p, label='Weymouth equation', color='black', linewidth = lw)
    # plt.show()

    # markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
    markers = ['X','o', 'o', 'X' ]
    # ["(FOP)NNOGF","PLA (BP10)","PLA (BP6)","(U)NNOGF B6"]
    markers = dict(zip(df_y.columns,[ 'X', 'X','o', 'o' ]))

    for  i,c in enumerate(df_y.columns):
        ax = sns.scatterplot(x = df_x[c], y = df_y[c], label=c,  marker=markers[c], s = size)
    
    ax.tick_params(axis='y', labelsize=fs-2)
    ax.tick_params(axis='x', labelsize=fs-2)

    plt.ylabel('q', fontsize=fs)
    plt.xlabel('$\Delta$p$^2$', fontsize=fs)

    # plt.grid(True)  # Add gridlines
    plt.legend)  # Add legend
    sns.despine()


    legend = plt.legend(fontsize='large',loc='upper left')
    for handle in legend.legendHandles:
        handle._sizes = [100]
        handle.legend_loc = 'upper left'

    save_path = start_directory
    
    plt.savefig(save_path+fr"\opt_theo_pipe{pipe}.png",bbox_inches="tight",)            

    plt.close()

opt_theo_plot(q_p1_x, q_p1, "1")


opt_theo_plot(q_p7_x, q_p7, "7")

opt_theo_plot(q_p20_x, q_p20, "20")




# mn = q_p7_x.min().min()*1.1
# mx = q_p7_x.max().max()*0.9

# sorted_indices = np.argsort(sqp_values[i])
# x_sorted = sqp_values[i][sorted_indices]
# y_theo_sorted = theo_qnm_values[i][sorted_indices]
# y_exp_sorted = exp_qnm_values[i][sorted_indices]

# x = np.linspace(mn, mx, 100)
# p = weymouth_equation(x)

# sns.set(style='ticks')

# # Create the plot
# plt.figure(figsize=(8, 6))  # Set the figure size



# # Plot the analytical method with sorted data
# sns.lineplot(x = x, y=  p, label='Analytical Method', color='black')
# # plt.show()

# # markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
# markers = ['X','o', 'o', 'X' ]

# for  i,c in enumerate(q_p7.columns):
#      sns.scatterplot(x = q_p7_x[c], y = q_p7[c], label=c,  marker=markers[i])
     
# plt.ylabel('q')
# plt.xlabel('$\Delta$p$^2$')

# # Customize the plot appearance
# plt.grid(True)  # Add gridlines
# plt.legend()  # Add legend
# sns.despine()
# plt.show()

# save_path = os.path.join(start_directory, fr'weymouth_{foldername}')
# if not os.path.exists(save_path):
#     # Create the folder
#     os.makedirs(save_path)
#     print(f"Folder created at '{save_path}'")
# plt.savefig(save_path+f"\{foldername}_pipe{i+1}.png",bbox_inches="tight",)            
# # Close the plot to free up memory
# plt.close()




# p_theo = excel_sqp>=0
# n_theo = excel_sqp<0

# theo_qnm_values_pos = np.sqrt(excel_sqp[p_theo])
# theo_qnm_values_neg = -np.sqrt(-excel_sqp[n_theo])
# theo_qnm_values = np.transpose(theo_qnm_values_pos.fillna(theo_qnm_values_neg).values).flatten()

# feas_values = np.transpose(excel_feas.values).flatten()
# sqp_values = np.transpose(excel_sqp.values).flatten()
# exp_qnm_values = np.transpose(excel_qnm.values/pipe_knm).flatten()












