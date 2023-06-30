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

#BENCHMARK
start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA"

start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\test"

#FEASIBLE TRAINNG DATA/EACH PIPE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\each_pipe_nb"


#BOUND TIGHTENING UNIFORM VS LINSPACE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Bound tightening"

pipe_knm = np.array([255.5169343, 255.5169343, 208.6287032,208.6287032,100.2219872,26.8636084,32.71143353,40.41306369,
                     68.90779278,114.2706469,13.94307458,102.2067737,12.47106503,78.85423786,80.8015493,
                    228.5412938,161.6030986,102.2067737,19.24327111,3.501408717,14.15077486])

# top = cm.get_cmap('Oranges_r', 128)
# bottom = cm.get_cmap('Blues', 128)

# newcolors = np.vstack((top(np.linspace(0, 1, 128)),
#                        bottom(np.linspace(0, 1, 128))))
# newcmp = ListedColormap(newcolors, name='OrangeBlue')

# excel_pressure = pd.read_excel(excel_file_path, sheet_name='pr')

# excel_pressure

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


            excel_sqp = pd.DataFrame()
            for p in p_map:
                excel_sqp[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)

            


            p_theo = excel_sqp>=0
            n_theo = excel_sqp<0

            theo_qnm_values_pos = np.sqrt(excel_sqp[p_theo])
            theo_qnm_values_neg = -np.sqrt(-excel_sqp[n_theo])
            theo_qnm_values = np.transpose(theo_qnm_values_pos.fillna(theo_qnm_values_neg).values)


            mask = excel_parameters["keys"] =="obj_val"

            feas_values = np.transpose(excel_feas.values)
            sqp_values = np.transpose(excel_sqp.values)
            exp_qnm_values = np.transpose(excel_qnm.values/pipe_knm)

            for i,p in enumerate(sqp_values):
                sorted_indices = np.argsort(p)
                x_sorted = sqp_values[i][sorted_indices]
                y_theo_sorted = theo_qnm_values[i][sorted_indices]
                y_exp_sorted = exp_qnm_values[i][sorted_indices]

                sns.set(style='ticks')

                # Create the plot
                plt.figure(figsize=(8, 6))  # Set the figure size

                # Plot the analytical method with sorted data
                sns.lineplot(x = x_sorted, y=  y_theo_sorted, label='Analytical Method', color='blue')

                # Plot the approximation (using original x values)
                sns.lineplot(x = x_sorted, y =  y_exp_sorted , label='Approximation', color='red', linestyle='--', marker='o')

                # Set plot title and labels
                plt.title('Comparison of Analytical Method and Approximation')
                plt.xlabel('x')
                plt.ylabel('y')

                # Customize the plot appearance
                plt.grid(True)  # Add gridlines
                plt.legend()  # Add legend
                sns.despine()
                
                
                save_path = os.path.join(start_directory, fr'weymouth_{foldername}')
                if not os.path.exists(save_path):
                    # Create the folder
                    os.makedirs(save_path)
                    print(f"Folder created at '{save_path}'")
                plt.savefig(save_path+f"\{foldername}_pipe{i+1}.png",bbox_inches="tight",)            
                # Close the plot to free up memory
                plt.close()




            p_theo = excel_sqp>=0
            n_theo = excel_sqp<0

            theo_qnm_values_pos = np.sqrt(excel_sqp[p_theo])
            theo_qnm_values_neg = -np.sqrt(-excel_sqp[n_theo])
            theo_qnm_values = np.transpose(theo_qnm_values_pos.fillna(theo_qnm_values_neg).values).flatten()

            feas_values = np.transpose(excel_feas.values).flatten()
            sqp_values = np.transpose(excel_sqp.values).flatten()
            exp_qnm_values = np.transpose(excel_qnm.values/pipe_knm).flatten()


            
            sorted_indices = np.argsort(sqp_values)
            x_sorted = sqp_values[sorted_indices]
            y_theo_sorted = theo_qnm_values[sorted_indices]
            y_exp_sorted = exp_qnm_values[sorted_indices]



            sns.set(style='ticks')

            # Create the plot
            plt.figure(figsize=(8, 6))  # Set the figure size

            # Plot the analytical method with sorted data
            sns.lineplot(x = x_sorted, y=  y_theo_sorted, label='Analytical Method', color='blue')

            # Plot the approximation (using original x values)
            sns.lineplot(x = x_sorted, y =  y_exp_sorted , label='Approximation', color='red', linestyle='--', marker='o')

            # Set plot title and labels
            plt.title('Comparison of Analytical Method and Approximation')
            plt.xlabel('x')
            plt.ylabel('y')

            # Customize the plot appearance
            plt.grid(True)  # Add gridlines
            plt.legend()  # Add legend
            sns.despine()
            
            
            save_path = os.path.join(start_directory, 'weymouth'+foldername+'.png')
            plt.savefig(save_path,bbox_inches="tight",)            
            # Close the plot to free up memory
            plt.close()





# def weymouth_equation(p1, p2):
#     q = []
#     for i,x in enumerate(p1):
#             print(p1[i])
#             if (p1[i]-p2[i])>=0:        
#                 print(np.sqrt(p1[i]**2 - p2[i]**2))
#                 q.append(np.sqrt(p1[i]**2 - p2[i]**2))
#             else:
#                 q.append(-np.sqrt(p1[i]**2 - p2[i]**2))


#     return  np.array(q)

# p = np.linspace(0, 2 * np.pi, 100)  # Range for x values
# x = p**2-p[::-1]**2
# analytical_y = weymouth_equation(x, x[::-1])  # Analytical solution for y values




# Create the plot
# plt.figure(figsize=(8, 6))  # Set the figure size

# # Plot the analytical method
# plt.plot(x_sorted, y_theo_sorted, label='Analytical Method', color='blue')

# # Plot the approximation
# plt.plot(x_sorted, y_exp_sorted, label='Approximation', color='red', linestyle='--')

# # Set plot title and labels
# plt.title('Comparison of Analytical Method and Approximation')
# plt.xlabel('x')
# plt.ylabel('y')

# # Customize the plot appearance
# plt.grid(True)  # Add gridlines
# plt.legend()  # Add legend
# plt.tight_layout()  # Adjust spacing

# plt.show()  # Show the plot

