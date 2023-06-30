import os
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator

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



def get_excel_data(excel_file_path, pipe_knm = pipe_knm, p_map = p_map):
    excel_pressure = pd.read_excel(excel_file_path, sheet_name='pr')

    excel_theo_wey = pd.DataFrame()
    for p in p_map:
        excel_theo_wey[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)*pipe_knm[p[0]-1]**2

    excel_qnm = pd.read_excel(excel_file_path, sheet_name='Q_nm')
    excel_qnm_k = excel_qnm/pipe_knm

    p_theo = excel_theo_wey>=0
    n_theo = excel_theo_wey<0

    excel_theo_p_q = np.sqrt(excel_theo_wey[p_theo]/pipe_knm**2)
    excel_theo_n_q = -np.sqrt(-excel_theo_wey[n_theo]/pipe_knm**2)
    excel_theo_all_q = excel_theo_p_q.fillna(excel_theo_n_q)



    return excel_theo_all_q, excel_qnm_k


def get_training_data(data_path):
    training_data = pd.read_csv(data_path)
    return training_data

#==================================================================================================


unnogfb6_results_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\uni_1m_b6_Av1\uni_1m_b6_Av1_results.xlsx"
uniform_training_data_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\Uniform_results\training_data_uniform_num_samples_1000_p_range_3_7_k_seed_42.csv"
unnogfb6_theo_q, unnogfb6_opt_q = get_excel_data(unnogfb6_results_path)
u_train_data = get_training_data(uniform_training_data_path)



# ==================================================================================================

FOP_results_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\each_pipe_nb\each_pipe_nb_results.xlsx"

FOP_training_data_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\feasible_operating_data\all_data_q.csv"


FOP_theo_q, FOP_opt_q = get_excel_data(FOP_results_path)
FOP_train_data = get_training_data(FOP_training_data_path)/pipe_knm







FOP_p1_TD = pd.DataFrame()
FOP_p1_TD["q"] = FOP_train_data[["pipe1"]]
FOP_p1_TD["Model"] = "(FOP)NNOGF"
FOP_p1_TD["Data"] = "training data"  #maybe call this  "q" instead of "data"?


FOP_p1_opt = pd.DataFrame()
FOP_p1_opt["q"] = FOP_opt_q[["1"]]
FOP_p1_opt["Model"] = "(FOP)NNOGF"
FOP_p1_opt["Data"] = "optimal data"

uni_p1_TD = pd.DataFrame()
uni_p1_TD["q"] = u_train_data[["q"]]
uni_p1_TD["Model"] = "(U)NNOGF"
uni_p1_TD["Data"] = "training data"

uni_p1_opt = pd.DataFrame()
uni_p1_opt["q"] = unnogfb6_theo_q[["1"]]
uni_p1_opt["Model"] = "(U)NNOGF"
uni_p1_opt["Data"] = "optimal data"


comp_data_p1 = pd.concat([FOP_p1_TD, FOP_p1_opt, uni_p1_TD, uni_p1_opt], axis=0)


sb.set(style = 'whitegrid') 
sb.violinplot(data=comp_data_p1, x="Model", y="q",order=["(U)NNOGF", "FOP"],  hue="Data", split=True, inner  =None,scale = "width", width = 0.5) #myaybe even less width?
plt.show()



FOP_p20_TD = pd.DataFrame()
FOP_p20_TD["q"] = FOP_train_data[["pipe20"]]
FOP_p20_TD["Model"] = "(FOP)NNOGF p20"
FOP_p20_TD["Data"] = "training data"  #maybe call this  "q" instead of "data"?


FOP_p20_opt = pd.DataFrame()
FOP_p20_opt["q"] = FOP_opt_q[["20"]]
FOP_p20_opt["Model"] = "(FOP)NNOGF p20"
FOP_p20_opt["Data"] = "optimal data"

uni_p20_TD = pd.DataFrame()
uni_p20_TD["q"] = u_train_data[["q"]]
uni_p20_TD["Model"] = "(U)NNOGF p20"
uni_p20_TD["Data"] = "training data"

uni_p20_opt = pd.DataFrame()
uni_p20_opt["q"] = unnogfb6_theo_q[["20"]]
uni_p20_opt["Model"] = "(U)NNOGF p20"
uni_p20_opt["Data"] = "optimal data"


comp_data_p20 = pd.concat([FOP_p20_TD, FOP_p20_opt, uni_p20_TD, uni_p20_opt], axis=0)

sb.set(style = 'whitegrid') 
sb.violinplot(data=comp_data_p20, x="Model", y="q",order=["(U)NNOGF p20", "FOP p20"],  hue="Data", split=True, inner  =None,scale = "width", width = 0.5) #myaybe even less width?
plt.show()



comp_data_p1_p20 = pd.concat([comp_data_p1, comp_data_p20], axis=0)

sb.set(style = 'whitegrid') 
ax = sb.violinplot(data=comp_data_p1_p20, x="Model", y="q",order=["(U)NNOGF", "FOP","(U)NNOGF p20", "FOP p20"],
                hue="Data", split=True, inner  =None,scale = "width", width = 0.5,
                linewidth=1
                ) #myaybe even less width?
# plt.show()
# sb.move_legend(ax, "lower right")
# sb.move_legend(ax, "best")
sb.despine(left=True)
sb.move_legend(ax, "lower right")
# plt.legend(loc='lower right')
plt.savefig(r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\comp_data_p1_p20.png", bbox_inches="tight")
plt.close()

# ==================================================================================================


qs =[0, 0.005,.01,.05,.1, .5, .9,0.99, 0.995, 1]


qs= [0.67 ,0.75, 0.875]

u_train_data[["q"]].quantile(qs)






0.67-0.875








#==================================================================================================
# (U)NNOGF
#==================================================================================================

excel_file_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\uni_1m_b6_Av1\uni_1m_b6_Av1_results.xlsx"

excel_pressure = pd.read_excel(excel_file_path, sheet_name='pr')

excel_theo_wey = pd.DataFrame()
for p in p_map:
    excel_theo_wey[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)*pipe_knm[p[0]-1]**2

excel_qnm = pd.read_excel(excel_file_path, sheet_name='Q_nm')
# excel_theo_wey== 0

p_theo = excel_theo_wey>=0
n_theo = excel_theo_wey<0

excel_theo_p_q = np.sqrt(excel_theo_wey[p_theo]/pipe_knm**2)
excel_theo_n_q = -np.sqrt(-excel_theo_wey[n_theo]/pipe_knm**2)
excel_theo_all_q = excel_theo_p_q.fillna(excel_theo_n_q)

# excel_theo_p_q.quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])

abs(excel_theo_all_q).mean()
excel_theo_all_q.quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])

save_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\Uniform_results"
uniform_data = pd.read_csv(os.path.join(save_path, "training_data_uniform_num_samples_1000_p_range_3_7_k_seed_42.csv"))
p = uniform_data[["q"]]>=0
n = uniform_data[["q"]]<=0
qs =[0, 0.005,.01,.05,.1, .5, .9,0.99, 0.995, 1]
uniform_data[p]["q"].quantile(qs)
uniform_data[n]["q"].quantile(qs)
uniform_data[["q"]].quantile(qs)
q = [ 0.475,.4875, 0.499,0.4995,0.49975, 0.5,.5125,.515,.5155, 0.525] #, .5375,  0.55]
uniform_data[["q"]].quantile(q)
uniform_data[["q"]].quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])


#==================================================================================================
# FOP
#==================================================================================================

excel_file_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA\each_pipe_nb\each_pipe_nb_results.xlsx"


excel_pressure = pd.read_excel(excel_file_path, sheet_name='pr')

excel_theo_wey = pd.DataFrame()
for p in p_map:
    excel_theo_wey[str(p[0])] = (excel_pressure[str(p[1])]**2-excel_pressure[str(p[2])]**2)*pipe_knm[p[0]-1]**2

excel_qnm = pd.read_excel(excel_file_path, sheet_name='Q_nm')
excel_qnm_k = excel_qnm/pipe_knm
# excel_theo_wey== 0

p_theo = excel_theo_wey>=0
n_theo = excel_theo_wey<0

excel_theo_p_q = np.sqrt(excel_theo_wey[p_theo]/pipe_knm**2)
excel_theo_n_q = -np.sqrt(-excel_theo_wey[n_theo]/pipe_knm**2)
excel_theo_all_q = excel_theo_p_q.fillna(excel_theo_n_q)


rel_err = (excel_qnm_k-excel_theo_all_q)/excel_theo_all_q

# excel_theo_p_q.quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])

abs(excel_theo_all_q).mean()
excel_theo_all_q.quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])

save_path = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Neural network training\feasible_operating_data\all_data_q.csv"
FOPdata = pd.read_csv(save_path)/pipe_knm
p = FOPdata>=0
n = FOPdata<=0
qs =[0, 0.005,.01,.05,.1, .5, .9,0.99, 0.995, 1]
FOPdata[p].quantile(qs)
FOPdata[n].quantile(qs)
FOPdata[["q"]].quantile(qs)
q = [ 0.475,.4875, 0.499,0.4995,0.49975, 0.5,.5125,.515,.5155, 0.525] #, .5375,  0.55]
FOPdata[["q"]].quantile(q)

FOPdata.quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])


excel_theo_all_q.quantile([0, 0.005,.01,.05,.1, .5, .9,0.95, 0.99, 0.995, 1])

FOPdata.hist(bins=100)

FOPdata.plot(kind='box', showfliers=False)
plt.show()

excel_theo_all_q.hist(bins=24)
plt.show()



excel_theo_all_q.plot(kind='box', showfliers=False)
plt.show()


# excel_theo_all_q.boxplot()
# plt.show()


# # Create a DataFrame with your data
# data = {
#     'Group 1': [10, 12, 14, 16, 18],
#     'Group 2': [8, 9, 12, 14, 16]
# }
# df = pd.DataFrame(data)

# # Create a boxplot
# ax = df.boxplot()

# # Compute the percentile value
# percentile_value = df['Group 1'].quantile(0.75)  # Change the column and percentile as needed

# # Add a vertical line at the percentile value
# ax.axvline(x=percentile_value, color='red', linestyle='--')

# # Set labels and title
# ax.set_xlabel('Groups')
# ax.set_ylabel('Values')
# ax.set_title('Boxplot with Percentile Marker')

# # Show the plot
# plt.show()