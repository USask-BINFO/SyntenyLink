import numpy as np
import pandas as pd
import re
import warnings
import sys
import pickle
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import os
import wandb
import time
import SyntenyLink_mbp as mbp
import accuracy_heatmap_weight_comb as acc
import weight_comb_Calc as wcc

# os.environ["WANDB_NOTEBOOK_NAME"] = "wandb.ipynb"
os.environ["WANDB_SILENT"] = "true"
time_stamp = time.strftime("%m%d-%H%M")
# start a new wandb run to track this script
wandb.init(project="brassica_parameters_sinapis", name=str(time_stamp))

def main(gap_thresholds, min_block_lengths, y = 0):
    #get the input file as an argument (collinear file)

    input_file = sys.argv[sys.argv.index('-i') + 1]

    #take gap_threshold and min_block_length which are numbers as arguments
    # gap_threshold = int(sys.argv[sys.argv.index('-g') + 1])
    # min_block_length = int(sys.argv[sys.argv.index('-m') + 1])

    #Get the number of subgenomes as an argument
    n_subgenomes = int(sys.argv[sys.argv.index('-n') + 1])

    #convert the collinear file to a dataframe
    C_df_csv = pd.read_csv(input_file, sep = '\t', header=None)
    #make a copy of the dataframe
    C_df = C_df_csv.copy() 
    C_df_with_head = pd.read_csv(input_file, sep = '\t')
    C_df_head = C_df_with_head.iloc[0:, 1:-3]
    #convert all the 'x' to 0 and all the other entries which are not equal to 'x' to 1 omitting first column and last three columns and first row
    # C_df_head = C_df_head.replace(r'^x$', 0, regex=True)
    # #convert all the other entries starts with a letter to 1
    # C_df_head = C_df_head.replace(r'^[a-zA-Z]', 1, regex=True)

    C_df_updated = C_df.iloc[1:, 1:-3]
    #convert all the 'x' to 0 and all the other entries which are not equal to 'x' to 1 omitting first column and last three columns and first row
    C_df_updated = C_df_updated.replace(r'^x$', 0, regex=True)
    #convert all the other entries starts with a letter to 1
    C_df_updated = C_df_updated.replace(r'^[a-zA-Z]', 1, regex=True)
    #set first row index starts from 0
    C_df_updated.index = C_df_updated.index - 1
    C_df_updated_copy = C_df.iloc[1:, 1:-3]
    C_df_updated_copy.columns = C_df.iloc[0,1:-3]
    C_df_updated_copy.index = C_df_updated_copy.index - 1

    #If there exist a ground truth file, then compare the results with the ground truth
    GT = sys.argv[sys.argv.index('-gt') + 1]

    # Split the dataframe based on the first letter of column names
    dfs = {}
    for column in C_df_head.columns:
        print(column)
        first_letter = column[0]
        if first_letter not in dfs:
            dfs[first_letter] = pd.DataFrame()
        dfs[first_letter][column] = C_df_head[column]

    # Print the resulting dataframes
    for key, value in dfs.items():
        print(f"Dataframe with columns starting with '{key}':")
        first_letter_get = key
        print(first_letter_get)
        print(value)
        print()
        #make C_df_new empty first
        C_df_new = pd.DataFrame()
        C_df_new = C_df_with_head.iloc[:, [0]].join(value).join(C_df_with_head.iloc[:, -3:])
        value = value.replace(r'^x$', 0, regex=True)
        #convert all the other entries starts with a letter to 1
        value = value.replace(r'^[a-zA-Z]', 1, regex=True)
        # Convert the dataframe to a numpy array
        C = value.to_numpy()
        print(C)
        m, n = C.shape
        print(m, n)
        print(y)
        gap_threshold = gap_thresholds[y]
        print(f"Gap threshold for dataframe {y+1}: {gap_threshold}")
        min_block_length = min_block_lengths[y]
        print(f"Minimum block length for dataframe {y+1}: {min_block_length}")
        y += 1
        # print(m, n)

        break_point_indices = mbp.gap_calculation(C, gap_threshold, min_block_length, n, m)
        densities = mbp.get_densities(C, break_point_indices, n, m)
        mbp.create_excel_sheet(C_df_new, break_point_indices, densities, value)
        df_temp = mbp.get_subgenomes("Super_synteny_block_output.xlsx", n_subgenomes)
        df = mbp.assign_subgenomes(df_temp, f"Super_synteny_block_output.xlsx", n_subgenomes)
        df.to_excel(f"Super_synteny_bl_sub_placement_density.xlsx")
        # Get the column names from df_subgenome_density that start with N followed by a number
        column_names = df.filter(regex=r'^N\d+').columns.tolist()

        # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
        df_synteny = C_df_new.iloc[0:,1:-3].rename(columns={C_df_new.iloc[0:,1:-3].columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})
        #append the first column of C_df to the first column of df_synteny
        df_synteny.insert(0, "locus_id", C_df_new.iloc[0:,0])
        #update first row index starting from 0
        # df_synteny.index = df_synteny.index - 1
        # acc.subgenome_overlap(GT,"Super_synteny_bl_sub_placement_density.xlsx", df_synteny, 3, first_letter_get)
        #Input the file for chains generated in DagChainer as a command line argument
        chains_file = sys.argv[sys.argv.index('-c') + 1]
        #Input the file for blastn generated in DagChainer as a command line argument
        blastn_file = sys.argv[sys.argv.index('-bl') + 1] 

        num_blocks_main = len(pd.read_excel("Super_synteny_bl_sub_placement_density.xlsx"))      

        max_parameters = wcc.get_weight_accuracy(n_subgenomes,"Super_synteny_bl_sub_placement_density.xlsx", "Super_synteny_graph_nodes_sub.xlsx", GT, df_synteny, chains_file, blastn_file, C_df_new)
        print(max_parameters)

# Get the gap thresholds as command line arguments
gap_thresholds = []
if '-g' in sys.argv:
    index = sys.argv.index('-g')
    for i in range(index+1, len(sys.argv)):
        if sys.argv[i].startswith('-'):
            break
        gap_thresholds.append(float(sys.argv[i]))

# Get the block size thresholds as command line arguments
min_block_lengths = []
if '-m' in sys.argv:
    index = sys.argv.index('-m')
    for i in range(index+1, len(sys.argv)):
        if sys.argv[i].startswith('-'):
            break
        min_block_lengths.append(int(sys.argv[i]))

if __name__ == '__main__':
    main(gap_thresholds, min_block_lengths, y=0)