import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
# Specify the directory to start the iteration from
# start_directory = 'path/to/start/directory'

#BENCHMARK
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA"

#FEASIBLE TRAINNG DATA/EACH PIPE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\each_pipe_nb"


#BOUND TIGHTENING UNIFORM VS LINSPACE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Bound tightening"



#used pla vs nnogf:
start_directory =r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\in_report"

#used architecture test:
# start_directory =r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\arch_tests\used_models"


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


# p_map = [[1,1,2],
#         [2,	2,3],
#         [3,	2,4]]

# pipe_knm = np.array([2.591141949,2.591141949,6.346975626])








# excel_theo_wey = pd.DataFrame()
# for p in p_map:
#     excel_theo_wey[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)*pipe_knm[p[0]-1]**2

err_dict = {}

# Iterate over subfolders
for foldername in os.listdir(start_directory):
    folder_path = os.path.join(start_directory, foldername)
    
    # Check if the item is a subfolder
    if os.path.isdir(folder_path):
        excel_file_path = os.path.join(folder_path, foldername +"_results"+ '.xlsx')
        
        # Check if the corresponding Excel file exists
        if os.path.isfile(excel_file_path):
            # Load the Excel file and read the 'feasibility_gap' sheet

            excel_parameters = pd.read_excel(excel_file_path, sheet_name='parameters')
            m = excel_parameters["keys"] == "mae"
            print(foldername)
            print(excel_parameters[m])


            excel_feas = pd.read_excel(excel_file_path, sheet_name='feasibility_gap')


            excel_qnm = pd.read_excel(excel_file_path, sheet_name='Q_nm')
            excel_parameters = pd.read_excel(excel_file_path, sheet_name='parameters')
            excel_pressure = pd.read_excel(excel_file_path, sheet_name='pr')

            excel_theo_wey = pd.DataFrame()
            for p in p_map:
                excel_theo_wey[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)*pipe_knm[p[0]-1]**2




            mask = excel_parameters["keys"] =="obj_val"
            # Extract the data values from the DataFrame
            feas_values = np.transpose(excel_feas.values)
            qnm_sq_values = np.transpose(excel_qnm.values*abs(excel_qnm.values)) #should this be ^2 ???
            # theo_values = qnm_sq_values-feas_values
            theo_values = np.transpose(excel_theo_wey.values)
            rel_error = ((abs(qnm_sq_values-theo_values)+1e-9)/(theo_values+1e-9))*100


            # print(foldername)
            # print(f'max: {np.max(abs(feas_values)/1000)}, min: {np.min(abs(feas_values)/1000)}, mean: {np.mean(abs(feas_values)/1000)}') #maybe take square root of this...
            print(f'rel_error; max: {np.max(rel_error)}, min: {np.min(rel_error)}, mean: {np.mean(rel_error)}')
            # print(f'obj:{excel_parameters[mask]}\n')
            # print(f'rel_error; max: {np.max(rel_error)}, min')
            # scale = 10000


            nrmse_v2 = np.sqrt(np.mean(rel_error)**2) 
            nrmse = np.sqrt(np.mean(rel_error**2))
            print(f'nrmse: {nrmse}, nrmse_v2: {nrmse_v2}')
            err_dict[foldername] = nrmse

            h = 8; w = 6
            fs = 18
            fig, ax = plt.subplots(figsize=(h, w))
            im_ratio = h/w
            # heatmap = ax.imshow(abs(feas_values)/scale, cmap='gist_yarg', interpolation='nearest', vmin=0, vmax=50000//scale) #absolut error
            # heatmap = ax.imshow(rel_error, interpolation='nearest', vmin=-4,vmax = 4) #relativ error
            # heatmap = ax.imshow(rel_error, cmap='bwr', interpolation='nearest', vmin=-1,vmax = 1) #relativ error
            # heatmap = ax.imshow(rel_error, cmap='coolwarm', interpolation='nearest', vmin=-1,vmax = 1) #relativ error
            heatmap = ax.imshow(rel_error, cmap='coolwarm', interpolation='nearest', vmin=-100,vmax = 100) #relativ error
            
            # ax.set_title("Absolute Error Heat Map - {}".format(foldername))
            ax.set_xlabel(r"Time [hours]", fontsize=fs)
            ax.set_ylabel(r"Pipe $\mathcal{P}$", fontsize=fs)
            # cbar = plt.colorbar(heatmap, shrink = 0.89)
            cbar = plt.colorbar(heatmap, fraction=0.047*im_ratio)
            # cbar = plt.colorbar(heatmap, extend='both', shrink=0.8, ticks=range(1, feas_values.max()+1))


            # Format the y-axis ticks
            # ax.set_yticks(np.arange(feas_values.shape[0]), labels=np.arange(1,feas_values.shape[0]+1)[::-1])
            ax.set_yticks(np.arange(feas_values.shape[0]), labels=np.arange(1,feas_values.shape[0]+1))
            ax.set_xticks(np.arange(feas_values.shape[1]), labels=np.arange(1,feas_values.shape[1]+1))

            # Set the x-axis ticks to show every 3rd tick
            ax.xaxis.set_major_locator(MultipleLocator(3))

            # Set the y-axis ticks to show every 3rd tick
            ax.yaxis.set_major_locator(MultipleLocator(4))
            ax.tick_params(axis='y', labelsize=fs-2) 
            ax.tick_params(axis='x', labelsize=fs-2) 
            cbar.ax.tick_params(labelsize=fs-2)


            # legend_text = r'Absolute error [MNm$^3$/h]'#check the units...
            # legend_text = r'Absolute error [kg$/h]'
            legend_text = r'Relative error [%]' # [MNm$^3$/h]'
            # cbar.ax.text(0.5, -0.1, legend_text, va='center', ha='center')
            cbar.ax.set_ylabel(legend_text, rotation=90, fontsize=fs)
            # cbar.ax.yaxis.set_label_coords(-.1, .5)
            # Format colorbar values as integers

            # plt.show()
            # Save the plot in the associated subfolder
            save_path = os.path.join(start_directory, 'heatmap_'+foldername+'.png')
            plt.savefig(save_path,bbox_inches="tight",)            
            # Close the plot to free up memory
            plt.close()







# print(feas_values)
# excel_data
# 8.099855885404184




# {
# 'G_case_Gas_tree_compr_CB_f_mt_NN_constrained_ps1': 0.3574593125620862, 

# 'seq_1l_4n':  0.3859170302593904, 
# 'seq_1l_10n': 0.2792529415612169, 
# 'seq_2l_2n':  0.4695616425808823, 
# 'seq_2l_5n':  0.3574593125620862,
#  'seq_2l_8n': 0.3108519789541147, 
# 'seq_3l_3n':  0.32482734209948505, 
# 'seq_3l_5n':  0.2373738481148868
# }




# {
# 'uni_1m_b4_Av1':        0.5829626824080487,
# 'pla_Av1_V6':           0.5436736635261323, 
# 'pla_Av1_nt3_V10':      0.4028391245645183,
# 'uni_1m_b6_Av1':        0.45654220457504, 
# 'FOP_each_pipe_Av1_v2': 0.28239578548380734, 



# 'uni_1m_b3_Av1_nt3':    0.4512698212402897, 
# 'uni_1m_b4_Av1':        0.5829626824080487,
# 'uni_1m_b6_Av1':        0.45654220457504, 
# 'uni_1m_nb_Av1':        0.4593950775860708
# }