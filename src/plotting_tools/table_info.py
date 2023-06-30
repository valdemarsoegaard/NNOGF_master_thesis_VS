import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#BOUND TIGHTENING UNIFORM VS LINSPACE
# start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Bound tightening"

start_directory = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\Belgian case study\Benchmark models\job_logs_nt6_NNOGF_vs_PLA"

data_dict = {}
# Iterate over subfolders
for foldername in os.listdir(start_directory):
    folder_path = os.path.join(start_directory, foldername)
    
    # Check if the item is a subfolder
    if os.path.isdir(folder_path):
        # excel_file_path = os.path.join(folder_path, foldername +"_results"+ '.xlsx')

        
        # Iterate over items in the subfolder
        items = os.listdir(folder_path)
        output_file_path_list = [x for x in items if '.out' in x]
        if len(output_file_path_list) > 0:
            # column_dict = {"column": output_file_path_list[0]}
            data = []
            print(output_file_path_list[0])
            output_file_path = os.path.join(folder_path, output_file_path_list[0])
            fp = open(output_file_path, "r")
            file_string = fp.read()
            if "Explored" in file_string:
                n_start = file_string.find("Explored")+9
                n_end = file_string.find("nodes")-1
                print(f"Nodes explored: {file_string[n_start:n_end]}")
                # column_dict["Nodes explored"] = file_string[n_start:n_end]
                data.append(file_string[n_start:n_end])
            if "Best objective" in file_string:
                n_start = file_string.find("Best objective")+15
                n_end = file_string.find("best bound")-2
                print(f"Best objective: {file_string[n_start:n_end]}")
                # column_dict["Best objective"] = file_string[n_start:n_end]
                data.append(file_string[n_start:n_end])
            if ") in" in file_string:
                n_start = file_string.find( ") in")+5
                inds = [i for i in range(len(file_string)) if file_string.startswith("seconds (", i)]
                n_end = inds[-2]
                print(f"Solution time: {file_string[n_start:n_end]}")
                # column_dict["Solution time"] = file_string[n_start:n_end]
                data.append(file_string[n_start:n_end])
            if "integer (" in file_string:
                n_start = file_string.find( "integer (")+9
                n_end = file_string.find("binary)")-1
                print(f"Number binaries: {file_string[n_start:n_end]}")
                # column_dict["Number binaries"] = file_string[n_start:n_end]
                data.append(file_string[n_start:n_end])
            print("\n")
            if len(data) == 4:
                data_dict[foldername] = data




# idx_dict = {"column": "column", "Nodes explored": "Nodes explored", "Best objective": "Best objective", "Solution time": "Solution time", "Number binaries": "Number binaries"}

df = pd.DataFrame.from_dict(data_dict, orient='index',
                       columns=['Nodes Explored', 'Expected Cost [€]', 'solution time', 'binaries'])
output_folder = r"C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\results\tables"
df.to_excel(os.path.join(output_folder, "table_info.xlsx"))
df.to_excel(os.path.join(start_directory, "table_info.xlsx"))