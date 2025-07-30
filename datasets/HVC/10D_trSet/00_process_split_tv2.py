import numpy as np
import math
import random
import sys
import subprocess
import scipy.stats as ss
import os
from os import path

import matplotlib.pyplot as plt

def dominates(a, b):
    return (np.asarray(a) <= b).all()

# Check number of parameters
n_params = len(sys.argv)
if n_params != 6:
    print("Use: "+sys.argv[0]+" DIM(2|3|4|...|10) N_SUBFOLDER OFFSET SUBFOLDER_SIZE SUBSUBFOLDERLETTER")
    exit(1)
str_dim = sys.argv[1]
subfolder = sys.argv[2]
n_subfolder = int(subfolder)
offset = int(sys.argv[3])
subfolder_size = int(sys.argv[4])
subsubfolder= sys.argv[5]

train_path = "ParetoFronts/SPLIT/" + str_dim + "D_T/"+ str_dim + "DT_" +subfolder+subsubfolder
output_path = train_path+"/FINAL"
file_prefix = "/0" + str_dim + "D"
aux_file_prefix = "/AUX/0" + str_dim +"D_R"

if not os.path.isdir(output_path):
    os.mkdir(output_path)
if not os.path.isdir(output_path+"/AUX"):
    os.mkdir(output_path+"/AUX")

# List all files within the directory
content = [ f for f in os.listdir(train_path) if os.path.isfile(os.path.join(train_path,f)) ]
content = sorted(content)
n_files = len(content)

# Process each file in the defined interval
for index in range(subfolder_size): # Recorre los archivos
    # Read input file
    file_index = offset + index
    input_filename = path.join(train_path,content[index])
    input_file = open(input_filename)
    lines = input_file.readlines()
    input_file.close()

    # Store data (assuming first line to be formatted as '# rows cols')
    n_lines = int(lines[0].split()[1])    # Retrieve the number of lines from the second token from the first line
    n_cols = int(lines[0].split()[2])    # Retrieve the number of columns from the third token from the first line
    data = np.zeros((n_lines, n_cols))
    max_vals = np.zeros(n_cols)
    min_vals = np.zeros(n_cols)
    for line_index in range(n_lines):
        for col_index in range(n_cols):
            line = lines[line_index+1].split()
            data[line_index, col_index] = float(line[col_index])
    #plt.scatter(data[:,0], data[:,1])
    #plt.show()

    # Remove duplicated rows
    data = np.unique(data, axis=0)
    n_lines = len(data)

    # Find non-dominated
    d_flag = np.zeros(n_lines)
    non_dominated = 0
    for i in range(n_lines):
        for j in range(n_lines):
            if( i!=j ):
                if dominates(data[j], data[i]):
                    d_flag[i] = 1
                    #print(str(j)+" dominates "+str(i))
        if d_flag[i] == 0:
            non_dominated += 1

    # Get rid of dominated individuals
    filtered_data = np.zeros((non_dominated, n_cols))
    aux_index = 0
    for i in range(n_lines):
        if d_flag[i] == 0:
            for j in range(n_cols):
                filtered_data[aux_index, j] = data[i,j]
            aux_index += 1
    data = filtered_data
    n_lines = non_dominated

    #plt.scatter(filtered_data[:,0], filtered_data[:,1])
    #plt.show()

    # Check if there are more than 2 non-dominated individuals
    if(n_lines > 2):
        # Normalize data
        max_vals = data.max(axis=0)
        min_vals = data.min(axis=0)
        for i in range(n_lines):
            for j in range(n_cols):
                data[i,j] = 0.01 + ( ((data[i,j]-min_vals[j]) * (0.99-0.01))/(max_vals[j]-min_vals[j]))
        #print(data)

        #plt.scatter(data[:,0], data[:,1])
        #plt.show()

        if file_index<9:
            out_filename = output_path+file_prefix + "_R0" + str(file_index+1) + ".hvc"
        else:
            out_filename = output_path+file_prefix + "_R" + str(file_index+1) + ".hvc"

        # Create auxiliary files without the i-th row line to obtain individual contributions
        for i in range(n_lines):
            if i<9:
                aux_filename = output_path+aux_file_prefix + "_R0" + str(i+1) + ".pof"
            else:
                aux_filename = output_path+aux_file_prefix + "_R" + str(i+1) + ".pof"

            aux_file = open(aux_filename, "w")
            aux_file.write("# "+str(n_lines-1)+" "+str(n_cols)+"\n")

            aux_lines = np.delete(data,i,0)

            for i in range(n_lines-1):
                for j in range(n_cols):
                    aux_file.write(str(aux_lines[i,j])+" ")
                aux_file.write("\n")
            aux_file.close()

        # Obtain the hypervolume of the auxiliary files
        try:
            # Obtain hypervolume using the cmd line:
            #./emo_indicator HV prefix n_lines ref_point
            # For instance:
            # ./emo_indicator HV ParetoFronts/M2/M2TRAINING_200_2102/AUX_FILES/AUX 7 1.01 1.01
            cmd = ["./emo_indicator", "HV", output_path+aux_file_prefix, str(n_lines)]
            for c in range(n_cols):
                cmd.append("1.01")
            #print(cmd)
            tmp=subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            print("signal error: "+str(e.output))

        # Retrieve hypervolume contributions
        hv_filename = output_path+aux_file_prefix+".hv"
        hv_file = open(hv_filename, "r")
        hv_lines = hv_file.readlines()
        hv_file.close()

        # Sort hypervolume contributions
        #print(hv_lines)
        hvc = np.zeros(n_lines)
        for i in range(n_lines):
            hvc[i] = float(hv_lines[i+1])
        hvc_ranked = ss.rankdata(hvc)

        # Create output file
        out_file = open(out_filename, "w")

        # Write first row
        out_file.write("# "+str(n_lines)+" "+str(n_cols)+" "+"\n")
        # Write remaining lines with their corresponding hv contribution
        for i in range(n_lines):
            line_str = ""
            for j in range(n_cols):
                line_str += str(data[i,j])+" "
            out_file.write(line_str+ hv_lines[i+1].replace("\n","")+" "+str(hvc_ranked[i])+"\n")
        out_file.close()
        print("output file created in "+out_filename)
    else:
        print("Not enough non dominated points in file: "+input_filename)
