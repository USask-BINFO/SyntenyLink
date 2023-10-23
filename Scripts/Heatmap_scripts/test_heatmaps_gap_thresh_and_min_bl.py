import numpy as np
import pandas as pd
import re
import warnings
import sys
import pickle
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import time
import os
import wandb

# os.environ["WANDB_NOTEBOOK_NAME"] = "wandb.ipynb"
# os.environ["WANDB_SILENT"] = "true"
# time_stamp = time.strftime("%m%d-%H%M")
# # start a new wandb run to track this script
# wandb.init(project="brassica_parameters_gt_mbl", name=str(time_stamp))

# function for finding gaps in each columns
def find_gaps(col, gap_threshold, m):
    """
    Finds the gaps in the each column of the matrix
    col: column of the matrix
    gap_threshold: threshold for the gap
    m: number of rows in the matrix

    returns: indices of the beginning and end of the gaps
    """
    # indices for ones in given column
    indices_one = np.where(col == 1)[0]
    # print(indices_one)
    # calculate difference between consective indices
    d = np.diff(indices_one)
    # print(d)

    # if the difference is bigger than the gap threshold then store them
    initial_gaps = np.where(d > gap_threshold)[0]
    # print(initial_gaps)
    # number of gaps
    r = len(initial_gaps)
    # calculate gaps
    gaps = np.zeros((r, 2))
    # find the indices of beginning and end of gaps
    gaps[:, 0] = indices_one[initial_gaps] + 1
    gaps[:, 1] = indices_one[initial_gaps + 1] - 1
    # fix the end of gap index for the last gap
    if r > 0 and indices_one[-1] + 1 < m:
        gaps = np.append(gaps, np.array([[indices_one[-1] + 1, m]]), axis=0)
    return gaps

#get the input file as an argument (collinear file)

input_file = sys.argv[1]

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
GT = sys.argv[5]
#get the number of subgenomes as a comand line argument
n_subgenomes = int(sys.argv[4])

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
    # print(m, n)


    # gap_threshold is the threshold to determine if there is a gap or not
    # min_block_length is the minimum length of the signal blocks

    def gap_calculation(C, gap_threshold, min_block_length):
        # get all gaps in all columns
        gaps_cell = [find_gaps(C[:, i], gap_threshold, m) for i in range(n)]
        break_point_indices = np.zeros(10000)
        # for the first iteration, we only look at the best two columns
        col_ind = np.argpartition(np.sum(C[:min_block_length, :], axis=0), -2)[-2:]
        # print(col_ind)

        if gaps_cell[col_ind[0]][0][0] <= gaps_cell[col_ind[1]][0][0]:
            i = int(gaps_cell[col_ind[0]][0][0])
        else:
            i = int(gaps_cell[col_ind[1]][0][0])
        # the first gap is our first breakpoint
        break_point_indices[0] = i
        # we have one break point, so count=1
        count = 1
        while i <= m - 1:
            i = int(i)
            col_ind = np.argpartition(
                np.sum(C[i: min(m, i + min_block_length), :], axis=0), -3
            )[-3:]
            i_0 = gaps_cell[col_ind[0]][gaps_cell[col_ind[0]][:, 0] > i]
            i_1 = gaps_cell[col_ind[1]][gaps_cell[col_ind[1]][:, 0] > i]
            i_2 = gaps_cell[col_ind[2]][gaps_cell[col_ind[2]][:, 0] > i]
            if len(i_0) == 0:
                i_0 = [m]
            if len(i_1) == 0:
                i_1 = [m]
            if len(i_2) == 0:
                i_2 = [m]
            i = min(np.min(i_0), np.min(i_1), np.min(i_2))
            if abs(i - break_point_indices[count - 1]) > gap_threshold:
                break_point_indices[count] = i
                count += 1
        break_point_indices = np.resize(break_point_indices, count)

        if break_point_indices[-1]!= m-2 and break_point_indices[-1]!= m-1 and break_point_indices[-1]!= m:
            break_point_indices = np.append(break_point_indices, m-2)
        if break_point_indices[-1]== m-1 or break_point_indices[-1]== m:
            break_point_indices = break_point_indices[:-1] + 1
        if break_point_indices[-1]!= m-1 or break_point_indices[-1]!= m:
            break_point_indices = break_point_indices + 1
        return break_point_indices


    def get_densities(C, break_point_indices, n , m):
        # calculate densities of each block
        b = len(break_point_indices)  # number of breakpoints
        densities = np.zeros((b, n))  # densities for each blocks
        for i in range(b):
            if i == 0:  # for block one
                densities[i, :] = (
                    np.sum(C[: int(break_point_indices[i]), :], axis=0)
                    / int(break_point_indices[i])
                ).astype(float)
            elif i == b:  # for last block
                densities[i, :] = (
                    np.sum(C[int(break_point_indices[i - 1]+1): m, :], axis=0)
                    / int(m - break_point_indices[i - 1]+1)
                ).astype(float)
            else:  # for all else
                densities[i, :] = (
                    np.sum(
                        C[int(break_point_indices[i - 1]+1): int(break_point_indices[i]), :],
                        axis=0,
                    )
                    / int(break_point_indices[i] - break_point_indices[i - 1]+1)
                ).astype(float)
        return densities

    def create_excel_sheet(C_df, break_point_indices, densities):
        # calculate row_end indices as in excel sheet
        row_end = (break_point_indices[: len(break_point_indices)]
        ).astype(int)
        # row_end=np.append(row_end[:-1]+1)
        # calculate row_start indices as in excel sheet
        row_start = (
            np.append(0, break_point_indices[: len(break_point_indices)-1] + 1)
        ).astype(int)
        # calculate number of genes
        number_genes = (row_end - row_start+1).astype(int)

        # read gene IDs from the given file
        gene_id = C_df.iloc[:,0]
        gene_id = gene_id.to_numpy()
        # prepare data for gene_start and gene_end columns for excel sheet
        gene_start = []
        gene_end = []
        b = len(break_point_indices)
        for i in range(b):
            gene_start.append(gene_id[row_start[i]])
            gene_end.append(gene_id[row_end[i]])

        E = pd.DataFrame()
        # prepare data matrix for excel sheet
        E = np.column_stack(
            (np.arange(1, b + 1),
                gene_start,
                gene_end,
                row_start,
                row_end,
                number_genes,
                densities,
            )
        )

        # Initialize the header list with the fixed column names
        header = ["Block no.",
            "Block_start",
            "Block_end",
            "Row start #",
            "Row end #",
            "# genes in block",
        ]

        # Extract the prefix letter from the column names using regex
        prefixes = set((re.findall(r'^[A-Za-z]', col)[0], int(re.findall(r'\d+', col)[0])) for col in value.columns)
        prefixes = sorted(set([f'N{suffix}' if suffix % 2 == 1 else f'N{suffix}.r' for _, suffix in prefixes]))

        # Sort the columns by the letter prefix and the number following the letter
        cols = sorted(value.columns, key=lambda col: (re.findall(r'^[A-Za-z]', col)[0], int(re.findall(r'\d+', col)[0])))

        # Loop through the columns and update the header list
        for i in range(0, len(cols), 2):
            n = (i // 2) + 1
            header.extend([f"N{n}", f"N{n}.r"])

        # write matrix to excel sheet
        df = pd.DataFrame(E, columns=header)
        df.to_excel(f"{first_letter_get}_gene_matrix_output_Bra.xlsx",
                    sheet_name="Output", startrow = 0, startcol=0)



    def get_subgenomes(filepath, num_subgenomes):
        '''
        This function takes the filepath of the excel sheet generated by the previous function and the number of subgenomes as input and returns a list of subgenomes
        
        Parameters:
        filepath (str): The filepath of the excel sheet generated by the previous function
        num_subgenomes (int): The number of subgenomes
        
        Returns:
        subgenomes (list): A list of subgenomes
        '''
        # Read the generated break point output
        df_temp = pd.read_excel(filepath, usecols="H:AA")

        # Get the number of cells in each row where it doesn't contain a value equals to zero and store it as a value in the next column of df
        df_temp["Non_zero"] = df_temp.iloc[:, :].apply(lambda x: x[x != 0].count(), axis=1)

        # Iterate over each row and get the columnname with max and have it as a value in the next column of df disregarding the last column
        df_temp["subgenome1"] = df_temp.iloc[:, :-1].apply(lambda x: x.sort_values(ascending=False).index[0], axis=1)

        # Iterate over each row and get the second highest value's column name disregarding the last column of each row as the second highest value in the next column of df
        df_temp["subgenome2"] = df_temp.iloc[:, :-2].apply(lambda x: x.sort_values(ascending=False).index[1], axis=1)

        if num_subgenomes > 2:
            for i in range(3, num_subgenomes+1):
                # Iterate over each row and get the i-th highest value's column name disregarding the last (i-1) columns of each row as the i-th highest value in the next column of df
                df_temp[f"subgenome{i}"] = df_temp.iloc[:, :-i].apply(lambda x: x.sort_values(ascending=False).index[i-1], axis=1)

        return df_temp


    # df_temp = get_subgenomes("gene_matrix_output_Bra.xlsx", n_subgenomes)

    def assign_subgenomes(df_temp, file_path, num_subgenomes):
        '''
        Assigns the subgenomes to the blocks based on the gene density inside the blocks
        
        Parameters
        ----------
        df_temp: pd.DataFrame
            The dataframe with the gene density inside the blocks
        file_path: str
            The path of the file to be read
        num_subgenomes: int
            The number of subgenomes to be assigned to the blocks
        Returns
        -------
        df: pd.DataFrame
            The dataframe with the subgenomes assigned to the blocks
        '''
        # Columns with density lower than 0.1 inside a block has been removed
        for col in df_temp.columns[:-4]:
            for i in range(len(df_temp)):
                if df_temp[col][i] < 0.1:
                    df_temp.at[i, col] = 0
        warnings.filterwarnings("ignore")

        # Remove last four columns
        df = df_temp[df_temp.columns[:-4]]


        # Assignment of the columns into subgenomes depending on gene density inside blocks
        df.loc[:, "Non_zero"] = df.iloc[:, :].apply(
            lambda x: x[x != 0].count(), axis=1)
        for i in range(1, num_subgenomes+1):
            df.loc[:, f"subgenome{i}"] = df.iloc[:, :-i].apply(
                lambda x: x.sort_values(ascending=False).index[i-1], axis=1
            )

        # Read the generated break point output
        df_new = pd.read_excel(file_path, usecols="A:G")
        df = pd.concat([df_new.reset_index(drop=True), df], axis=1)

    # Assignment of the columns into subgenomes depending on gene density inside blocks
        for i in range(len(df)):
            if i == 0 and df["Non_zero"][i] == num_subgenomes - 1:
                if df["Non_zero"][i + 1] == num_subgenomes - 1:
                    df.at[i, f"subgenome{num_subgenomes}"] = df.at[i + 2, f"subgenome{num_subgenomes}"]
                if df["Non_zero"][i + 1] > num_subgenomes - 1:
                    df.at[i, f"subgenome{num_subgenomes}"] = df.at[i + 1, f"subgenome{num_subgenomes}"]
            elif df["Non_zero"][i] <= num_subgenomes - 1:
                if df["Non_zero"][i] == num_subgenomes - 1:
                    if df.at[i, f"subgenome{n_subgenomes-1}"] != df.at[i-1, f"subgenome{n_subgenomes}"] and df.at[i, f"subgenome{n_subgenomes-2}"] != df.at[i-1, f"subgenome{n_subgenomes}"]:
                        df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
                    elif df.at[i, f"subgenome{n_subgenomes-1}"] == df.at[i-1, f"subgenome{n_subgenomes}"] and df.at[i, f"subgenome{n_subgenomes-2}"] != df.at[i-1, f"subgenome{n_subgenomes}"]:
                        if df.at[i-1, f"subgenome{n_subgenomes-1}"] != df.at[i, f"subgenome{n_subgenomes-2}"]:
                            df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                            df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
                        elif df.at[i-1, f"subgenome{n_subgenomes-1}"] == df.at[i, f"subgenome{n_subgenomes-2}"]:
                            df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-2}"]
                            df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]

                    elif df.at[i, f"subgenome{n_subgenomes-2}"] == df.at[i-1, f"subgenome{n_subgenomes}"] and df.at[i, f"subgenome{n_subgenomes-1}"] != df.at[i-1, f"subgenome{n_subgenomes}"]:
                        if df.at[i-1, f"subgenome{n_subgenomes-2}"] != df.at[i, f"subgenome{n_subgenomes-1}"]:
                            df.at[i, f"subgenome{n_subgenomes-2}"] = df.at[i - 1, f"subgenome{n_subgenomes-2}"]
                            df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
                        elif df.at[i-1, f"subgenome{n_subgenomes-2}"] == df.at[i, f"subgenome{n_subgenomes-1}"]:
                            df.at[i, f"subgenome{n_subgenomes-2}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                            df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
            
                if df["Non_zero"][i] == num_subgenomes - 2:
                    if df.at[i, f"subgenome{n_subgenomes-2}"] != df.at[i-1, f"subgenome{n_subgenomes}"] and df.at[i, f"subgenome{n_subgenomes-2}"] != df.at[i-1, f"subgenome{n_subgenomes-1}"]:
                        df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                        df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
                    elif df.at[i, f"subgenome{n_subgenomes-2}"] == df.at[i-1, f"subgenome{n_subgenomes}"]:
                        df.at[i, f"subgenome{n_subgenomes-2}"] = df.at[i - 1, f"subgenome{n_subgenomes-2}"]
                        df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                        df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
                    elif df.at[i, f"subgenome{n_subgenomes-2}"] == df.at[i-1, f"subgenome{n_subgenomes-1}"]:
                        df.at[i, f"subgenome{n_subgenomes-2}"] = df.at[i - 1, f"subgenome{n_subgenomes-2}"]
                        df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                        df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
                    
                if df["Non_zero"][i] == num_subgenomes - 3:
                    df.at[i, f"subgenome{n_subgenomes-2}"] = df.at[i - 1, f"subgenome{n_subgenomes-2}"]
                    df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                    df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
        return df

        # df = assign_subgenomes(df_temp, "gene_matrix_output_Bra.xlsx", n_subgenomes)


        # # Saving the file
        # df.to_excel("subgenome_density_bra.xlsx")

        #get another input file for groundtruth as a command line argument

        # # Get the column names from df_subgenome_density that start with N followed by a number
        # column_names = df.filter(regex=r'^N\d+').columns.tolist()

        # # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
        # df_synteny = C_df_csv.iloc[1:,1:-3].rename(columns={C_df_csv.iloc[1:,1:-3].columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})
        # #append the first column of C_df to the first column of df_synteny
        # df_synteny.insert(0, "locus_id", C_df_csv.iloc[1:,0])
        # #update first row index starting from 0
        # df_synteny.index = df_synteny.index - 1

    def subgenome_overlap(GT_file_path, df_subgenome_density_file_path, df_synteny, num_subgenomes, sub_no):
        # Read the data
        df_subgenome_density = pd.read_excel(df_subgenome_density_file_path)
        # Create a list of subgenomes
        subgenomes = [[] for _ in range(num_subgenomes)]
        for i in range(len(df_subgenome_density["Row start #"])):
            for j in range(df_subgenome_density["Row start #"].values[i], df_subgenome_density["Row end #"].values[i] + 1):
                for k in range(num_subgenomes):
                    subgenomes[k].append(df_synteny[df_subgenome_density[f"subgenome{k+1}"].values[i]][j])

        # Create a dataframe of subgenomes
        df_subgenome = pd.DataFrame(subgenomes).transpose()
        # Change the column names
        df_subgenome.columns = [f"subgenome{i+1}" for i in range(num_subgenomes)]

        # Load the dataframes
        df_groundtruth = pd.read_excel(GT_file_path) # dataframe for groundtruth


        # Check for exact matches of gene ids and 'x' between subgenome1 of both dataframes
        exact_matches = 0
        total_genes = 0
        missing_genes = 0
        missing_genes_list = []

        for i, row_groundtruth in df_groundtruth.iterrows():
            if i < len(df_groundtruth) - 1:
                    gene_id_groundtruth = row_groundtruth[f"{first_letter_get}_subgenome{sub_no+1}"]
                    if gene_id_groundtruth != "x":
                        total_genes += 1
                        if i < len(df_subgenome):
                            row_subgenome = df_subgenome.iloc[i]
                            gene_id_subgenome = row_subgenome[f"subgenome{sub_no+1}"]
                            if gene_id_groundtruth == gene_id_subgenome:
                                exact_matches += 1
                            else:
                                missing_genes += 1
                                missing_genes_list.append(gene_id_groundtruth)

        # Calculate the overlap percentage for each subgenome
        overlaps = exact_matches / total_genes
        print(overlaps, sub_no+1)
        return overlaps

    #take gap_threshold and min_block_length which are numbers as arguments
    gap_threshold = int(sys.argv[2])
    min_block_length = int(sys.argv[3])

        # break_point_indices = gap_calculation(C, gap_threshold, min_block_length)
        # densities = get_densities(C, break_point_indices, n, m)
        # create_excel_sheet(C_df_new, break_point_indices, densities)
        # df_temp = get_subgenomes(f"{first_letter_get}_gene_matrix_output_Bra.xlsx", n_subgenomes)
        # df = assign_subgenomes(df_temp, f"{first_letter}_gene_matrix_output_Bra.xlsx", n_subgenomes)
        # df.to_excel(f"{first_letter_get}_subgenome_density_bra.xlsx")
        # # Get the column names from df_subgenome_density that start with N followed by a number
        # column_names = df.filter(regex=r'^N\d+').columns.tolist()

        # # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
        # df_synteny = C_df_new.iloc[0:,1:-3].rename(columns={C_df_new.iloc[0:,1:-3].columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})
        # #append the first column of C_df to the first column of df_synteny
        # df_synteny.insert(0, "locus_id", C_df_new.iloc[0:,0])
        # #update first row index starting from 0
        # df_synteny.index = df_synteny.index - 1

    #Input the file for chains generated in DagChainer as a command line argument
    chains_file = sys.argv[6]
    #Input the file for blastp generated in DagChainer as a command line argument
    blastp_file = sys.argv[7]       

    def generate_heatmap(acc_matrix, subgenome_index, gap_threshold, min_block_length):
        sub_acc_matrix = acc_matrix[subgenome_index]
        sns.set(font_scale=1.2)  # Increase the font scale to 1.2 or any desired value
        plt.figure(figsize=(8, 10))
        ax = sns.heatmap(sub_acc_matrix, cmap='YlGnBu', annot=True, fmt=".2f", cbar=False,
                        xticklabels=range(10, gap_threshold+1, 5), yticklabels=range(50, min_block_length+1, 10))
        plt.title("Subgenome " + str(subgenome_index + 1) + " Accuracy Matrix")
        plt.xlabel("Gap Threshold")
        plt.ylabel("Minimum Block Length")
        plt.savefig("subgenome_" + str(subgenome_index+1) + str(first_letter_get)+"_heatmap.png", dpi=300, bbox_inches='tight')



    def get_breakpoints(C, C_df_new, gap_threshold, min_block_length, n_subgenomes, GT, n, m):
        best_params = []
        accuracies = []
        for k in range(n_subgenomes):
            row_k = []
            for j in range(50, min_block_length+1, 10):
                row_j = []
                for i in range(10, gap_threshold+1, 5):
                    break_point_indices = gap_calculation(C, i, j)
                    densities = get_densities(C, break_point_indices, n, m)
                    create_excel_sheet(C_df_new, break_point_indices, densities)
                    df_temp = get_subgenomes(f"{first_letter_get}_gene_matrix_output_Bra.xlsx", n_subgenomes)
                    df = assign_subgenomes(df_temp, f"{first_letter_get}_gene_matrix_output_Bra.xlsx", n_subgenomes)
                    # Saving the file
                    df.to_excel(f"{first_letter_get}_subgenome_density_bra.xlsx")
                    
                    df = pd.read_excel(f"{first_letter_get}_subgenome_density_bra.xlsx")
                    # Get the column names from df_subgenome_density that start with N followed by a number
                    column_names = df.filter(regex=r'^N\d+').columns.tolist()

                    df_synteny = pd.DataFrame()
                    # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
                    df_synteny = C_df_new.iloc[0:,1:-3].rename(columns={C_df_new.iloc[0:,1:-3].columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})
                    #append the first column of C_df to the first column of df_synteny
                    df_synteny.insert(0, "locus_id", C_df_new.iloc[0:,0])
                    #update first row index starting from 0
                    # df_synteny.index = df_synteny.index - 1

                    acc = subgenome_overlap(GT, f"{first_letter_get}_subgenome_density_bra.xlsx", df_synteny, n_subgenomes, k)
                    print(acc)
                # if k == 0:
                #     wandb.log({"acc_sub1": acc, "bl_sub1": j, "gt_sub1": i, "subgenome": k+1})
                # elif k == 1:
                #     wandb.log({"acc_sub2": acc, "bl_sub2": j, "gt_sub2": i, "subgenome": k+1})
                # elif k == 2:
                #     wandb.log({"acc_sub3": acc, "bl_sub3": j, "gt_sub3": i, "subgenome": k+1})
                    row_j.append(acc)
                row_k.append(row_j)
            accuracies.append(row_k)
        print(accuracies)

        # Generate heatmap for subgenome k
        for l in range(n_subgenomes):
            generate_heatmap(accuracies, l, gap_threshold, min_block_length)
            print(accuracies)
            
            
        # Get the highest accuracy value and its corresponding gap threshold and minimum block length values for subgenome k
        max_acc = max([max(row) for row in accuracies[k]])
        max_indices = [(i, j) for i, row in enumerate(accuracies[k]) for j, val in enumerate(row) if val == max_acc]
        best_params.append((k, max_acc, max_indices))

        return best_params


    # df = pd.read_excel("subgenome_density_bra.xlsx")
    # # Get the column names from df_subgenome_density that start with N followed by a number
    # column_names = df.filter(regex=r'^N\d+').columns.tolist()

    # # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
    # df_synteny = C_df_csv.iloc[1:,1:-3].rename(columns={C_df_csv.iloc[1:,1:-3].columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})
    # #append the first column of C_df to the first column of df_synteny
    # df_synteny.insert(0, "locus_id", C_df_csv.iloc[1:,0])
    # #update first row index starting from 0
    # df_synteny.index = df_synteny.index - 1
    max_parameters = get_breakpoints(C, C_df_new, gap_threshold, min_block_length, n_subgenomes, GT, n , m)
    print(max_parameters)

#call the function to get breakpoints
# break_point_indices = get_breakpoints(C, gap_threshold, min_block_length)


# def get_blocks(df_subgenome_density, df_synteny, num_blocks, num_subgenomes):
#     '''
#     This function takes in the subgenome density dataframe, the synteny dataframe, the number of blocks and the number of subgenomes
#     and returns a dictionary with the blocks as the keys and the subgenomes as the values
    
#     Parameters
#     ----------
#     df_subgenome_density: pandas dataframe
#         The dataframe containing the subgenome density information
#         df_synteny: pandas dataframe
#         The dataframe containing the synteny information
#         num_blocks: int
#         The number of blocks
#         num_subgenomes: int
#         The number of subgenomes
        
#     Returns
#     -------
#     blocks: dictionary
#     The dictionary containing the blocks as the keys and the subgenomes as the values
#     '''
#     # store the blocks in a dictionary with the block number as the key
#     blocks = {}
#     for i in range(num_blocks):
#         blocks[i] = {}
#         for j in range(num_subgenomes):
#             subgenome_col = "subgenome" + str(j+1)
#             block_key = str(i) + "_" + df_subgenome_density[subgenome_col].values[i]
#             start_row = df_subgenome_density.loc[i, "Row start #"]
#             end_row = df_subgenome_density.loc[i, "Row end #"]
#             block_val = df_synteny.loc[start_row:end_row, df_subgenome_density.loc[i, subgenome_col]]
#             genes = df_synteny.loc[start_row:end_row, 'gene_id'].tolist()
#             block_val_dict = {}
#             for k in range(len(genes)):
#                 sub_key = genes[k]
#                 block_val_dict[sub_key] = block_val.iloc[k]
#             blocks[i][block_key] = block_val_dict
#     return blocks
    



# def create_chains_dict(df_chains):
#     """
#     Creates a dictionary with keys as the zeroth column in df_chains and values as the list of values belonging to the same
#     keys in the fourth column values in df_chains.
    
#     Args:
#     df_chains: pandas dataframe containing four columns
    
#     Returns:
#     chains: dictionary with keys as the zeroth column in df_chains and values as the list of values belonging to the same
#     keys in the fourth column values in df_chains.
#     """
#     chains = {}
#     for i in range(len(df_chains[0])):
#         chains[df_chains[0][i]] = []
#     for i in range(len(df_chains[0])):
#         chains[df_chains[0][i]].append(df_chains[2][i])
#     return chains

# def calculate_density(blocks_dict, num_blocks):
#     """
#     Calculates the density of each block in each subgenome.
    
#     Parameters:
#         blocks_dict (dict): A dictionary containing all the blocks and their genes for each subgenome.
#         num_blocks (int): The total number of blocks.
    
#     Returns:
#         dict: A dictionary containing the density of each block for each subgenome.
#     """
#     density_dict = {}
#     for i in range(num_blocks):
#         density_dict[i] = {}
#         for j in range(len(blocks_dict[i])):
#             subgenome_blocks = list(blocks_dict[i].values())[j]
#             num_genes = len([gene for gene in subgenome_blocks.values() if gene != "x"])
#             block_size = len(subgenome_blocks)
#             density_dict[i][list(blocks_dict[i].keys())[j]] = (
#                 num_genes / block_size if block_size > 0 else 0
#             )
#     return density_dict

# def calculate_pid(blocks, df, num_blocks, pickle_file):
#     '''
#     Calculates the PID of each block in each subgenome.
    
#     Parameters:
#         blocks (dict): A dictionary containing all the blocks and their genes for each subgenome.
#         df (pandas dataframe): A dataframe containing the PID values for each gene pair.
#         num_blocks (int): The total number of blocks.
#         pickle_file (str): The path to the pickle file containing the PID values for each gene pair.
        
#     Returns:
#     dict: A dictionary containing the PID of each block for each subgenome.
#     '''
#     # #Remove dot and number after dot in gene_id_AT
#     # df["gene_id_AT"] = [re.sub(r"\.\d+", "", gene_id) for gene_id in df["gene_id_AT"]]

#     # # Initialize the dictionary "pid"
#     # pid = {}
#     # for i in range(num_blocks):
#     #     pid[i] = {}
#     #     for j in range(len(blocks[i])):
#     #         key = list(blocks[i].keys())[j]
#     #         block = blocks[i][key]
#     #         pid_values = []
#     #         for sub_key in block.keys():
#     #             if sub_key in df["gene_id_AT"].tolist():
#     #                 val = block[sub_key]
#     #                 # Create a dictionary to store gene_id as keys and PID as values
#     #                 gene_id_to_pid = dict(zip(zip(df["gene_id_Brassica"], df["gene_id_AT"]), df["PID"]))
#     #                 if val!="x":
#     #                     brassica_id = df.loc[(df["gene_id_AT"] == sub_key) & (df["gene_id_Brassica"] == val), "gene_id_Brassica"].values
#     #                     pid_val = gene_id_to_pid.get((brassica_id[0], sub_key), "x") if len(brassica_id) > 0 and val != "x" else "x"
#     #                 else:
#     #                     pid_val = "x"
#     #                 pid_values.append(pid_val)
#     #         pid[i][key] = pid_values
#     # #save in pickle file
#     # with open(pickle_file, "wb") as f:
#     #     pickle.dump(pid, f)

#     #load from pickle file
#     with open(pickle_file, "rb") as f:
#         pid = pickle.load(f)
#     return pid

# def calculate_avg_pid(num_blocks, blocks, pid):
#     avg_pid = {}
#     for i in range(num_blocks):
#         avg_pid[i] = {}
#         for j in range(len(blocks[i])):
#             avg_pid[i][list(blocks[i].keys())[j]] = 0
#             if len([x for x in pid[i][list(blocks[i].keys())[j]] if x != "x"]) != 0:
#                 avg_pid[i][list(blocks[i].keys())[j]] = (
#                     sum([x for x in pid[i][list(blocks[i].keys())[j]] if x != "x"])
#                     / len([x for x in pid[i][list(blocks[i].keys())[j]]])
#                 ) / 100
#     return avg_pid


# def process_files_weight_cal(chains_file, blastp_file, df_subgenome_density_file_path, df_subgenomes_file_path, synteny_file, pickle_file):
#     '''
#     Processes the files and returns the required dataframes.
    
#     Parameters:
#         chains_file (str): Path to the chains file.
#         blastp_file (str): Path to the blastp file.
#         df_subgenome_density_file_path (str): Path to the subgenome density file.
#         synteny_file (str): Path to the synteny file.
#         pickle_file (str): Path to the pickle file.
        
#     Returns:
#     blocks_dict: A dictionary containing all the blocks and their genes for each subgenome.
#     density_dict: A dictionary containing the density of each block for each subgenome.
#     chains_dict: A dictionary containing the chains of each block for each subgenome.
#     pid_dict: A dictionary containing the PID of each block for each subgenome.
#     avg_pid_dict: A dictionary containing the average PID of each block for each subgenome.
#     '''

#     # read the blastp file
#     blastp = pd.read_csv(blastp_file, sep="\t", header=None)
#     df_blastp = blastp.iloc[:, 0:3]
#     #have column names for the blastp file as gene_id_Brassica, gene_id_AT, PID
#     df_blastp.columns = ["gene_id_Brassica", "gene_id_AT", "PID"]

#     # read the chains file
#     df_chains = pd.read_csv(chains_file, sep="\t", header=None)

#     # get the column names from df_subgenome_density that start with N followed by a number
#     df_subgenome = pd.read_excel(df_subgenome_density_file_path)
#     column_names = ['gene_id'] + df_subgenome.filter(regex=r'^N\d+').columns.tolist()
#     df_subgenome_density = pd.read_excel(df_subgenomes_file_path)

#     df_synteny = synteny_file.rename(columns={synteny_file.columns[i]: column_names[i] if i > 0 else "gene_id" for i in range(len(column_names))})

#     # Start the index from 0
#     df_synteny.index = df_synteny.index - 1

#     # number of blocks
#     num_blocks = len(df_subgenome_density)
#     print("===========================================")
#     print(f"Number of blocks: {num_blocks}")
#     print("===========================================")




#     blocks= get_blocks(df_subgenome_density, df_synteny, num_blocks, n_subgenomes)
#     density = calculate_density(blocks, num_blocks)
#     chains_dict = create_chains_dict(df_chains)
#     pid = calculate_pid(blocks, df_blastp, num_blocks, pickle_file)
#     avg_pid = calculate_avg_pid(num_blocks, blocks, pid)
    
#     return (blocks, density, chains_dict, pid, avg_pid)

# graph_input = process_files_weight_cal(chains_file, blastp_file, "subgenome_density_bra.xlsx", "subgenome_density_bra.xlsx", C_df_csv.iloc[1:,:-3], "pid_genes_removed_blastn.pickle")

# def create_block_graph(blocks, density, chains, avg_pid, num_subgenomes):
#     '''
#     Creates a graph of blocks.
    
#     Parameters:
#         blocks (dict): A dictionary containing all the blocks and their genes for each subgenome.
#         density (dict): A dictionary containing the density of each block for each subgenome.
#         chains (dict): A dictionary containing the chains of each block for each subgenome.
#         avg_pid (dict): A dictionary containing the average PID of each block for each subgenome.
#         num_subgenomes (int): Number of subgenomes.
    
#     Returns:
#     graph (dict): A dictionary containing the graph of blocks.
#     '''
    
#     graph = {}
#     num_blocks = len(blocks)
#     for i in range(num_blocks - 1):
#         for k in range(num_subgenomes):
#             graph[list(blocks[i].keys())[k]] = {}
#         for s in range(num_subgenomes):
#             graph[list(blocks[i + 1].keys())[s]] = {}
#         for j in range(len(chains)):
#             for l in range(num_subgenomes):
#                 for m in range(num_subgenomes):
#                     if (
#                         list(blocks[i].keys())[l].split("_")[1]
#                         == list(blocks[i + 1].keys())[m].split("_")[1]
#                         or list(blocks[i].keys())[l].split("_")[1].split(".")[0]
#                         == list(blocks[i + 1].keys())[l].split("_")[1].split(".")[0]
#                     ):
#                         graph[list(blocks[i].keys())[l]][
#                             list(blocks[i + 1].keys())[m]
#                         ] = (
#                             list(density[i].values())[l]
#                             - (
#                                 list(density[i].values())[l]
#                                 - list(density[i + 1].values())[m]
#                             )
#                         ) + (
#                             list(avg_pid[i].values())[l]
#                             - (
#                                 list(avg_pid[i].values())[l]
#                                 - list(avg_pid[i + 1].values())[m]
#                             )
#                         ) 
#                         new_list = (
#                             list(blocks[i][list(blocks[i].keys())[l]].values())
#                             + list(blocks[i + 1][list(blocks[i + 1].keys())[m]].values())
#                         )
#                         new_list = list(filter(lambda a: a != "x", new_list))
#                         # concatenate the two blocks

#                         if (
#                             list(blocks[i].keys())[l].split("_")[1].split(".")[0]
#                             == list(chains.keys())[j].split("_")[0]
#                             or list(blocks[i].keys())[l].split("_")[1]
#                             == list(chains.keys())[j].split("_")[0]
#                         ):
#                             if (
#                                 len(set(new_list) & set(
#                                     chains[list(chains.keys())[j]]))
#                                 / float(
#                                     len(set(new_list) | set(
#                                         chains[list(chains.keys())[j]]))
#                                 )
#                                 * 100
#                                 > 20
#                             ):
#                                 graph[list(blocks[i].keys())[l]][
#                                     list(blocks[i + 1].keys())[m]
#                                 ] = 2
#                             elif len(graph[list(blocks[i].keys())[l]]) != num_subgenomes:
#                                 graph[list(blocks[i].keys())[l]][
#                                     list(blocks[i + 1].keys())[m]
#                                 ] = (
#                                     list(density[i].values())[l]
#                                     - (
#                                         list(density[i].values())[l]
#                                         - list(density[i + 1].values())[m]
#                                     )
#                                 ) + (
#                                     list(avg_pid[i].values())[l]
#                                     - (
#                                         list(avg_pid[i].values())[l]
#                                         - list(avg_pid[i + 1].values())[m]
#                                     )
#                                 )

#                     elif (
#                         list(blocks[i].keys())[l].split("_")[1]
#                         == list(blocks[i + 1].keys())[m].split("_")[1].split(".")[0]
#                         or list(blocks[i + 1].keys())[l].split("_")[1]
#                         == list(blocks[i].keys())[m].split("_")[1].split(".")[0]
#                     ):
#                         graph[list(blocks[i].keys())[l]][
#                             list(blocks[i + 1].keys())[m]
#                         ] = (
#                             list(density[i].values())[l]
#                             - (
#                                 list(density[i].values())[l]
#                                 - list(density[i + 1].values())[m]
#                             )
#                         ) + (
#                             list(avg_pid[i].values())[l]
#                             - (
#                                 list(avg_pid[i].values())[l]
#                                 - list(avg_pid[i + 1].values())[m]
#                             )
#                         ) 
#                     else:
#                         graph[list(blocks[i].keys())[l]][list(blocks[i + 1].keys())[m]] = (
#                             list(density[i].values())[l]
#                             - (
#                                 list(density[i].values())[l]
#                                 - list(density[i + 1].values())[m]
#                             )
#                         ) + (
#                             list(avg_pid[i].values())[l]
#                             - (
#                                 list(avg_pid[i].values())[l]
#                                 - list(avg_pid[i + 1].values())[m]
#                             )
#                         )
                                
#     return graph


# # function for finding gaps in each columns

# def add_edge_weights(blocks, graph, weight, num_blocks):
#     '''
#     This function adds edge weights to the graph
    
#     Parameters
#     ----------
#     blocks : list
#         list of blocks
#     graph : dict
#         graph of blocks
#     weight : int
#         weight of the edge
    
#     Returns
#     -------
#     graph : dict
#         graph of blocks with edge weights
#     '''
#     for i in range(num_blocks - 1):
#         for node in blocks[i].keys():
#             max_value = -1
#             max_next_node = None
#             for next_node in blocks[i + 1].keys():
#                 if (
#                     node.split("_")[1] == next_node.split("_")[1]
#                     or node.split("_")[1] == next_node.split("_")[1].split(".")[0]
#                     or next_node.split("_")[1] == node.split("_")[1].split(".")[0]
#                 ):
#                     if graph[node][next_node] > max_value:
#                         max_value = graph[node][next_node]
#                         max_next_node = next_node
#             if max_next_node is not None:
#                 current_edge_value = max_value
#                 new_edge_value = current_edge_value + weight
#                 graph[node][max_next_node] = new_edge_value

#     return graph


# def get_nodes(start_node, end_node, num_blocks, graph):
#     '''
#     This function returns the nodes between the start and end node
    
#     Parameters
#     ----------
#     start_node : str
#         start node
#     end_node : str
#         end node
#     num_blocks : int
#         number of blocks
#     graph : dict
#         graph of blocks
    
#     Returns
#     -------
#     nodes_sub : list
#         list of nodes between the start and end node
#     '''
#     nodes_sub = []
#     for i in range(num_blocks - 1):
#         while start_node != end_node:
#             nodes_sub.append(start_node)
#             start_node = list(list(graph.values())[list(graph.keys()).index(start_node)])[
#                 list(
#                     list(graph.values())[
#                         list(graph.keys()).index(start_node)].values()
#                 ).index(
#                     max(
#                         list(
#                             list(graph.values())[
#                                 list(graph.keys()).index(start_node)
#                             ].values()
#                         )
#                     )
#                 )
#             ]
#             break
#     nodes_sub.append(end_node)
#     return nodes_sub


# def remove_nodes_from_graph(nodes_sub, graph, blocks):
#     '''
#     This function removes the nodes from the graph
    
#     Parameters
#     ----------
#     nodes_sub : list
#         list of nodes between the start and end node
#     graph : dict
#         graph of blocks
#     blocks : list
#         list of blocks
    
#     Returns
#     -------
#     graph : dict
#         graph of blocks with nodes removed
#     blocks : list
#         list of blocks with nodes removed
#     '''
#     for i in range(len(nodes_sub)):
#         graph.pop(nodes_sub[i], None)
#         for j in range(len(graph)):
#             graph[list(graph.keys())[j]].pop(nodes_sub[i], None)
#     # remove blocks in blocks[] matching the nodes in nodes[]
#     for i in range(len(nodes_sub)):
#         for j in range(len(blocks)):
#             blocks[j].pop(nodes_sub[i], None)
#     return graph, blocks


# def update_graph_edges(num_blocks, blocks, graph, weight):
#     '''
#     This function updates the graph edges
    
#     Parameters
#     ----------
#     num_blocks : int
#         number of blocks
#     blocks : list
#         list of blocks
#     graph : dict
#         graph of blocks
#     weight : int
#         weight of the edge
    
#     Returns
#     -------
#     graph : dict
#         graph of blocks with updated edges
#     '''
#     for i in range(num_blocks - 1):
#         for node in blocks[i].keys():
#             max_value = -1
#             max_next_node = None
#             for next_node in blocks[i + 1].keys():
#                 if (
#                     node.split("_")[1] == next_node.split("_")[1]
#                     or node.split("_")[1] == next_node.split("_")[1].split(".")[0]
#                     or next_node.split("_")[1] == node.split("_")[1].split(".")[0]
#                 ):
#                     if graph[node][next_node] > max_value:
#                         max_value = graph[node][next_node]
#                         max_next_node = next_node
#             if max_next_node is not None:
#                 current_edge_value = max_value
#                 new_edge_value = current_edge_value - weight
#                 graph[node][max_next_node] = new_edge_value
#     return graph

# def node_traverse_most_weighted_path(n_subgenomes, subgenome_density_file_path, nodes, W1, W2):
#     '''
#     This function traverses the graph to find the most weighted path
    
#     Parameters
#     ----------
#     n_subgenomes : int
#         number of subgenomes
        
#     Returns
#     -------
#     nodes_df : pandas dataframe
#         dataframe with the nodes_sub lists as columns
#     '''
#     # Get input from command line
#     weight_1 = float(W1)
#     weight_2 = float(W2)

#     # Get input from files
#     subgenome_density = pd.read_excel(subgenome_density_file_path)
#     # Create graph with edge weights and get nodes for first subgenome
#     graph = add_edge_weights(graph_input[0], create_block_graph(graph_input[0], graph_input[1], graph_input[2], graph_input[4], n_subgenomes), weight_1, len(graph_input[0]))
#     #get the first value in subgenome1 column (first row) in subgenome_density file as start node
#     start_node = f"0_{subgenome_density.iloc[0]['subgenome1']}"
#     #get the last value in subgenome1 column (last row) in subgenome_density file as end node
#     end_node = f"{len(subgenome_density)-1}_{subgenome_density.iloc[-1]['subgenome1']}"
#     nodes_sub1 = get_nodes(start_node, end_node, len(graph_input[0]), graph)

#     # Remove nodes from graph and get nodes for remaining subgenomes
#     blocks = graph_input[0]
#     nodes_df = pd.DataFrame({"nodes_sub1": nodes_sub1}) # Define the nodes_df variable here
#     for i in range(n_subgenomes - 1):
#         graph, blocks = remove_nodes_from_graph(nodes_sub1, graph, blocks)
#         graph = update_graph_edges(len(blocks), blocks, graph, weight_2)
#         start_node = f"0_{subgenome_density.iloc[0][f'subgenome{i+2}']}"
#         end_node = f"{len(subgenome_density)-1}_{subgenome_density.iloc[-1][f'subgenome{i+2}']}"
#         nodes_sub = get_nodes(start_node, end_node, len(blocks), graph)
#         nodes_sub_name = f"nodes_sub{i+2}"
#         nodes_df[nodes_sub_name] = nodes_sub
#         remove_nodes_from_graph(nodes_sub, graph, blocks)

#     # rename the columns as subgenomes
#     nodes_df.columns = [f"subgenome{i}" for i in range(1, n_subgenomes+1)]

#     #Append Row start # and Row end # columns from subgenome_density_bra as first two columns of the dataframe
#     nodes_df.insert(0, "Row start #", subgenome_density["Row start #"])
#     nodes_df.insert(1, "Row end #", subgenome_density["Row end #"])

#     #remove everything before _ in the nodes_df subgenome columns
#     for i in range(1, n_subgenomes+1):
#         column_name = f"subgenome{i}"
#         nodes_df[column_name] = nodes_df[column_name].str.split("_", n=1, expand=True)[1]

#     nodes_df.to_excel(nodes, index=False)

#     return nodes_df


# node_traverse_most_weighted_path(n_subgenomes, "subgenome_density_bra.xlsx", "nodes.xlsx", sys.argv[8], sys.argv[9])
# print("===========================================")
# print("Accuracy of subgenome assignment(Weighted_Graph_Approach(Main_BP)):")
# print("===========================================")
# subgenome_overlap(GT, "nodes.xlsx", df_synteny, 3)



# def find_small_gaps(col, gap_threshold, m, row_start):
#     """
#     Finds the gaps in the each column of the matrix
#     col: column of the matrix
#     gap_threshold: threshold for the gap
#     m: number of rows in the matrix

#     returns: indices of the beginning and end of the gaps
#     """
#     # indices for ones in given column
#     indices_one = np.where(col == 1)[0]

#     # calculate difference between consective indices
#     d = np.diff(indices_one)

#     # if the difference is bigger than the gap threshold then store them
#     initial_gaps = np.where(d > gap_threshold)[0]

#     # number of gaps
#     r = len(initial_gaps)

#     # calculate gaps
#     gaps = np.zeros((r, 2))

#     # find the indices of beginning and end of gaps
#     gaps[:, 0] = indices_one[initial_gaps] + 1 + row_start
#     gaps[:, 1] = indices_one[initial_gaps + 1] - 1 + row_start

#     # add row_start as the 1st column and indices_one[initial_gaps[0]] + 1 + row_start-1 as the second column for the first entry in gaps[]
#     if r == 1:
#         gaps[0, 0] = row_start
#         gaps[0, 1] = (indices_one[initial_gaps[0]] + 1 + row_start) - 2
#     if r > 1:
#         # add without replacing the first entry
#         gaps = np.insert(gaps, 0, np.array(
#             [[row_start, gaps[0, 0] - 2]]), axis=0)
#         gaps = np.insert(gaps, 1, np.array(
#             [[gaps[1, 0] - 1, gaps[1, 0] - 1]]), axis=0)

#     # fix the end of gap index for the last gap
#     if r > 0 and indices_one[-1] + 1 < m:
#         gaps = np.append(
#             gaps,
#             np.array([[indices_one[-1] + 1 + row_start, m + row_start - 1]]),
#             axis=0,
#         )
#     if r > 0 and indices_one[-1] + 1 == m:
#         gaps = np.append(
#             gaps, np.array([[m + row_start - 1, m + row_start - 1]]), axis=0
#         )

#     return gaps

# def get_gaps(breakpoints_file, sub_density_file, C, num_subgenomes, gap_threshold=1):
#     '''
#     Returns the gaps in the blocks
    
#     breakpoints_file: file containing the main breakpoints
#     sub_density_file: file containing the subgenome density and chromosome blocks
#     num_subgenomes: number of subgenomes
#     gap_threshold: threshold for the gap
    
#     returns: list of gaps
    
#     '''
#     breakpoints_main = pd.read_excel(breakpoints_file, usecols="H:AA")
#     sub_density = pd.read_excel(sub_density_file)
#     gaps = []
    
#     for i in range(len(breakpoints_main)):
#         for j in range(len(breakpoints_main.columns)):
#             if breakpoints_main.iloc[i, j] < 0.1 and breakpoints_main.iloc[i, j] != 0:
#                 # subgenome_names = [sub_density[f"subgenome{k+1}"][i].split(".")[0] for k in range(num_subgenomes)]
#                 # subgenome_combinations = [[subgenome_names[k], subgenome_names[l]] for k in range(num_subgenomes-1) for l in range(k+1, num_subgenomes)]
#                 if (
#                 (
#                     sub_density["subgenome1"][i].split(".")[0]
#                     != sub_density["subgenome2"][i]
#                 )
#                 and (
#                     sub_density["subgenome1"][i].split(".")[0]
#                     != sub_density["subgenome3"][i]
#                 )
#                 and (
#                     sub_density["subgenome3"][i].split(".")[0]
#                     != sub_density["subgenome2"][i]
#                 )
#                 and (
#                     sub_density["subgenome2"][i].split(".")[0]
#                     != sub_density["subgenome1"][i]
#                 )
#                 and (
#                     sub_density["subgenome3"][i].split(".")[0]
#                     != sub_density["subgenome1"][i]
#                 )
#                 and (
#                     sub_density["subgenome2"][i].split(".")[0]
#                     != sub_density["subgenome3"][i]
#                 )
#             ):

#                     # print("Block no.:", i+1, "Row no.:", sub_density['Row start #'][i], "Subgenome:", breakpoints_main.columns [j], "Density:", breakpoints_main.iloc[i,j])
#                     # gap_threshold is the threshold to determine if there is a gap or not
#                     # min_block_length is the minimum length of the signal blocks
#                     gap_threshold = 1
#                     # min_block_length = 10
#                     # get all gaps in all columns
#                     gaps_cell = [
#                         find_small_gaps(
#                             C[
#                                 sub_density["Row start #"][i]: sub_density["Row end #"][i]
#                                 + 1,
#                                 j,
#                             ],
#                             gap_threshold,
#                             sub_density["# genes in block"][i],
#                             sub_density["Row start #"][i],
#                         )
#                     ]
#                     if len(gaps_cell[0]) != 0:
#                         gaps.append((gaps_cell, i, j))
#     return gaps

# def merge_gaps(gaps):
#     '''
#     Merges gaps that are adjacent to each other.

#     Parameters
#     ----------
#     gaps : list
#         A list of tuples, where each tuple contains a list of arrays, the row number of the gap, and the column number of the gap.

#     Returns
#     -------
#     merged_gaps : list
#         A list of tuples, where each tuple contains a list of arrays, the row number of the gap, and the column number of the gap.
        
#     '''
#     merged_gaps = []
#     i = 0
#     while i < len(gaps) - 1:
#         if gaps[i][1] == gaps[i + 1][1]:
#             merged_row = []
#             j = 0
#             k = 0
#             while j < len(gaps[i][0][0]) and k < len(gaps[i + 1][0][0]):
#                 if j == 0 and k == 0:
#                     if gaps[i][0][0][0][1] > gaps[i + 1][0][0][0][1]:
#                         merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                         k += 1

#                     elif gaps[i][0][0][0][1] < gaps[i + 1][0][0][0][1]:
#                         merged_row.append(np.array(gaps[i][0][0][j]))
#                         j += 1

#                     elif gaps[i][0][0][0][1] == gaps[i + 1][0][0][0][1]:
#                         merged_row.append(np.array(gaps[i][0][0][j]))
#                         j += 1
#                         k += 1

#                 elif (
#                     j > 0
#                     and gaps[i + 1][0][0][k][0] == gaps[i][0][0][j][1] + 1
#                     and k > 0
#                 ):
#                     merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                     k += 1
#                 elif (
#                     j > 0
#                     and k > 0
#                     and gaps[i + 1][0][0][k][0] != gaps[i][0][0][j][1] + 1
#                     and gaps[i + 1][0][0][k][0] < gaps[i][0][0][j][0]
#                     and gaps[i + 1][0][0][k][1] < gaps[i][0][0][j][1]
#                 ):
#                     merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                     k += 1
#                 else:
#                     merged_row.append(np.array(gaps[i][0][0][j]))
#                     j += 1

#             while j < len(gaps[i][0][0]):
#                 merged_row.append(np.array(gaps[i][0][0][j]))
#                 j += 1
#             while k < len(gaps[i + 1][0][0]):
#                 merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                 k += 1
#             merged_gaps.append(([np.array(merged_row)], gaps[i][1]))
#             i += 3
#         elif (
#             gaps[i][1] == gaps[i + 1][1]
#             and gaps[i][1] != gaps[i + 2][1]
#             and gaps[i][1] != gaps[i - 1][1]
#         ):
#             merged_row = []
#             j = 0
#             k = 0
#             while j < len(gaps[i][0][0]) and k < len(gaps[i + 1][0][0]):
#                 if j == 0 and k == 0:
#                     if gaps[i][0][0][0][1] > gaps[i + 1][0][0][0][1]:
#                         merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                         k += 1

#                     elif gaps[i][0][0][0][1] < gaps[i + 1][0][0][0][1]:
#                         merged_row.append(np.array(gaps[i][0][0][j]))
#                         j += 1

#                     elif gaps[i][0][0][0][1] == gaps[i + 1][0][0][0][1]:
#                         merged_row.append(np.array(gaps[i][0][0][j]))
#                         j += 1
#                         k += 1

#                 elif (
#                     j > 0
#                     and gaps[i + 1][0][0][k][0] == gaps[i][0][0][j][1] + 1
#                     and k > 0
#                 ):
#                     merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                     k += 1
#                 elif (
#                     j > 0
#                     and k > 0
#                     and gaps[i + 1][0][0][k][0] != gaps[i][0][0][j][1] + 1
#                     and gaps[i + 1][0][0][k][0] < gaps[i][0][0][j][0]
#                     and gaps[i + 1][0][0][k][1] < gaps[i][0][0][j][1]
#                 ):
#                     merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                     k += 1
#                 else:
#                     merged_row.append(np.array(gaps[i][0][0][j]))
#                     j += 1

#             while j < len(gaps[i][0][0]):
#                 merged_row.append(np.array(gaps[i][0][0][j]))
#                 j += 1
#             while k < len(gaps[i + 1][0][0]):
#                 merged_row.append(np.array(gaps[i + 1][0][0][k]))
#                 k += 1
#             merged_gaps.append(
#                 ([np.array(merged_row)], gaps[i][1], gaps[i][2]))
#             i += 2
#         else:
#             merged_gaps.append(gaps[i])
#             i += 1
#     if i == len(gaps) - 1:
#         merged_gaps.append(gaps[i])
#     return merged_gaps



# def sort_gaps(gaps):
#     '''
#     Sort the gaps in each row of the gap list
    
#     Parameters
#     ----------
#     gaps : list
#         A list of tuples. Each tuple contains a list of arrays and a string.
#         The list of arrays contains the gaps in each row of the gene matrix.
#         The string is the name of the chromosome.
        
#     Returns
#     -------
#     sorted_gaps : list
#         A list of tuples. Each tuple contains a list of arrays and a string.
#         The list of arrays contains the gaps in each row of the gene matrix.
#         The string is the name of the chromosome.
#         The gaps in each row are sorted by the start position of the gap.
#     '''
#     sorted_gaps = []
#     for gap in gaps:
#         sorted_gap = []
#         for row in gap[0]:
#             sorted_gap.append(row[row[:, 0].argsort()])
#         sorted_gaps.append((sorted_gap, gap[1]))
#     return sorted_gaps

# def find_small_breakpoints(gaps, sub_density):
#     '''
#     Find the breakpoints that are smaller than the minimum gap size.
    
#     Parameters
#     ----------
#     gaps : list
#         A list of tuples. Each tuple contains a list of arrays and a string.
#         The list of arrays contains the gaps in each row of the gene matrix.
#         The string is the name of the chromosome.
        
#     Returns
#     -------
#     small_breakpoints : dict
#         A dictionary of dictionaries. The keys are the index of the gap
#         and the index of the row in the gap. The values are dictionaries
#         containing the start and end positions of the breakpoint and the
#         main breakpoint.
#     '''
#     sub_density = pd.read_excel(sub_density)
#     small_breakpoints = {}
#     for i in range(len(gaps)):
#         for j in range(len(gaps[i][0][0])):
#             if len(gaps[i][0][0]) == 2:
#                 rowstart = gaps[i][0][0][j][0]
#                 row_end = gaps[i][0][0][j][1]
#                 if j == 0:
#                     max_bp = sub_density["subgenome1"][gaps[i][1]]
#                     small_breakpoints[i, j] = {
#                         "Row start #": rowstart,
#                         "Row end #": row_end,
#                         "main_bp": gaps[i][1],
#                     }

#                 if j == 1 and rowstart > gaps[i][0][0][0][1] + 1:
#                     rowstart = gaps[i][0][0][j - 1][1] + 1
#                     row_end = gaps[i][0][0][j][0] - 1

#                     small_breakpoints[i, j] = {
#                         "Row start #": rowstart,
#                         "Row end #": row_end,
#                         "main_bp": gaps[i][1],
#                     }

#                 if j == len(gaps[i][0][0]) - 1:
#                     rowstart = gaps[i][0][0][j][0]
#                     row_end = gaps[i][0][0][j][1]
#                     small_breakpoints[i, j + 1] = {
#                         "Row start #": rowstart,
#                         "Row end #": row_end,
#                         "main_bp": gaps[i][1],
#                     }

#             if len(gaps[i][0][0]) > 2:
#                 rowstart = gaps[i][0][0][j][0]
#                 row_end = gaps[i][0][0][j][1]
#                 if j == 0:
#                     small_breakpoints[i, j] = {
#                         "Row start #": rowstart,
#                         "Row end #": row_end,
#                         "main_bp": gaps[i][1],
#                     }
#                 if j == 1 and rowstart == gaps[i][0][0][0][1] + 1:
#                     rowstart = gaps[i][0][0][j][1]
#                     row_end = gaps[i][0][0][j][0]
#                     small_breakpoints[i, j] = {
#                         "Row start #": rowstart,
#                         "Row end #": row_end,
#                         "main_bp": gaps[i][1],
#                     }

#                 if (j <= (len(gaps[i][0][0]) - 1)) and (j >= 1):
#                     if gaps[i][0][0][j][0] > gaps[i][0][0][j - 1][1] + 1:
#                         rowstart = gaps[i][0][0][j - 1][1] + 1
#                         row_end = gaps[i][0][0][j][0] - 1
#                         small_breakpoints[i, j] = {
#                             "Row start #": rowstart,
#                             "Row end #": row_end,
#                             "main_bp": gaps[i][1],
#                         }

#                     rowstart = gaps[i][0][0][j][0]
#                     row_end = gaps[i][0][0][j][1]
#                     small_breakpoints[i, j, 0] = {
#                         "Row start #": rowstart,
#                         "Row end #": row_end,
#                         "main_bp": gaps[i][1],
#                     }

#     return small_breakpoints


# def update_small_breakpoints(small_breakpoints):
#     '''
#     Update the small breakpoints.
    
#     Parameters
#     ----------
#     small_breakpoints : dict
#         A dictionary of dictionaries. The keys are the index of the gap
#         and the index of the row in the gap. The values are dictionaries
#         containing the start and end positions of the breakpoint and the
#         main breakpoint.
        
#     Returns
#     -------
#     small_breakpoints_new : dict
#         A dictionary of dictionaries. The keys are the index of the gap
#         and the index of the row in the gap. The values are dictionaries
#         containing the start and end positions of the breakpoint and the
#         main breakpoint.
#     '''
#     small_breakpoints = {
#         i: small_breakpoints[key] for i, key in enumerate(small_breakpoints)
#     }
#     small_breakpoints
#     # Check if the row start and row end numbers are like row start number = from the previous row end number + 1 and row end number = from the next row start number - 1, if not, update the row start and row end numbers.
#     small_breakpoints_new = small_breakpoints.copy()
#     for i in range(len(small_breakpoints) - 1):
#         if i == 0:
#             small_breakpoints_new[i] = small_breakpoints[i]
#         if i == len(small_breakpoints) - 1:
#             small_breakpoints_new[i] = small_breakpoints[i]
#         if i > 0 and i < len(small_breakpoints) - 1:
#             # if small_breakpoints[i]['Row end #'] == small_breakpoints[i+1]['Row end #'] and small_breakpoints[i]['Row start #'] == small_breakpoints[i+1]['Row start #']:
#             #     continue
#             if (
#                 small_breakpoints[i]["Row end #"]
#                 != small_breakpoints[i + 1]["Row start #"] - 1
#             ):
#                 # Check whether you can find this row end number in any of the row end numbers coming in upcoming rows.
#                 if (
#                     small_breakpoints[i]["Row start #"]
#                     != small_breakpoints[i + 1]["Row start #"]
#                     and small_breakpoints[i]["Row end #"]
#                     != small_breakpoints[i + 1]["Row end #"]
#                 ):
#                     row_end_list = []
#                     for s in range(i + 1, len(small_breakpoints)):
#                         row_end_list.append(small_breakpoints[s]["Row end #"])
#                     if small_breakpoints[i]["Row end #"] in row_end_list:
#                         # If you can find this row end number in any of the row end numbers coming in upcoming rows, then update the current row's end number.
#                         small_breakpoints_new[i] = {
#                             "Row start #": small_breakpoints[i]["Row start #"],
#                             "Row end #": small_breakpoints[i + 1]["Row start #"] - 1,
#                             "main_bp": small_breakpoints[i]["main_bp"],
#                         }

#                     elif (
#                         small_breakpoints[i]["Row end #"] not in row_end_list
#                         and small_breakpoints[i]["Row end #"]
#                         > small_breakpoints[i + 1]["Row start #"]
#                     ):
#                         small_breakpoints_new[i] = {
#                             "Row start #": small_breakpoints[i]["Row start #"],
#                             "Row end #": small_breakpoints[i + 1]["Row start #"] - 1,
#                             "main_bp": small_breakpoints[i]["main_bp"],
#                         }
#                         # If you can't find this row end number in any of the row end numbers coming in upcoming rows, then add it where row end number < row end number of next row but row end number > row end number of the row of focus and row start number of next row < row end number.
#                         for j in range(i + 1, len(small_breakpoints) - 1):
#                             if (
#                                 small_breakpoints[i]["Row end #"]
#                                 < small_breakpoints[j + 1]["Row start #"]
#                                 and small_breakpoints[i]["Row end #"]
#                                 > small_breakpoints[j]["Row start #"]
#                             ):
#                                 # Add this row end number and row start number of next row as a new row.
#                                 new_row = {
#                                     "Row start #": small_breakpoints[j]["Row start #"],
#                                     "Row end #": small_breakpoints[i]["Row end #"],
#                                     "main_bp": small_breakpoints[i]["main_bp"],
#                                 }
#                                 small_breakpoints_new[i + 1, j] = new_row

#                 if (
#                     small_breakpoints[i]["Row start #"]
#                     == small_breakpoints[i + 1]["Row start #"]
#                     and small_breakpoints[i]["Row end #"]
#                     != small_breakpoints[i + 1]["Row end #"]
#                     and small_breakpoints[i]["Row end #"]
#                     <= small_breakpoints[i + 1]["Row end #"]
#                 ):
#                     small_breakpoints_new[i] = small_breakpoints[i]
#                     small_breakpoints_new[i + 1] = {
#                         "Row start #": small_breakpoints[i]["Row end #"] + 1,
#                         "Row end #": small_breakpoints[i + 1]["Row end #"],
#                         "main_bp": small_breakpoints[i + 1]["main_bp"],
#                     }
#                 else:
#                     continue

#             if (
#                 small_breakpoints[i]["Row end #"]
#                 == small_breakpoints[i + 1]["Row start #"] - 1
#             ):
#                 if i < len(small_breakpoints) - 2:
#                     if (
#                         small_breakpoints[i + 1]["Row start #"]
#                         == small_breakpoints[i + 2]["Row start #"]
#                         and small_breakpoints[i + 1]["Row end #"]
#                         > small_breakpoints[i + 2]["Row start #"]
#                     ):
#                         small_breakpoints_new[i] = small_breakpoints[i]
#                         small_breakpoints_new[i + 1] = small_breakpoints[i + 2]

#     small_breakpoints = {
#         i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
#     }

#     return small_breakpoints

# def merge_main_small_breakpoints(small_breakpoints_new, sub_density_file):
#     '''
#     This function merges the main breakpoints and small breakpoints.
    
#     Parameters
#     ----------
#     small_breakpoints_new : dict
#         Dictionary of small breakpoints.
#         sub_density_file : str
#         Path to the sub density file.
        
#     Returns
#     -------
#     small_breakpoints_copy : dict
#     Dictionary of merged breakpoints.
#     '''
#     sub_density = pd.read_excel(sub_density_file)
#     small_breakpoints = {
#         i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
#     }

#     sub_density = (
#         sub_density.copy()
#     )  # make a copy to avoid modifying the original DataFrame
#     small_breakpoints_copy = {}  # make a copy to avoid modifying the original dictionary

#     for i in range(len(small_breakpoints) - 1):
#         if i == 0:
#             small_breakpoints_copy[i] = small_breakpoints[i]
#         if i == len(small_breakpoints) - 1:
#             small_breakpoints_copy[i] = small_breakpoints[i]
#         for j in range(len(sub_density)):
#             if i > 0 and i < len(small_breakpoints) - 1:
#                 if (
#                     small_breakpoints[i - 1]["Row end #"] + 1
#                     == small_breakpoints[i]["Row start #"]
#                 ):
#                     small_breakpoints_copy[i] = small_breakpoints[i]
#                 if (
#                     small_breakpoints[i - 1]["Row end #"] + 1
#                     != small_breakpoints[i]["Row start #"]
#                     and small_breakpoints[i - 1]["Row start #"]
#                     == sub_density["Row start #"][j]
#                 ):
#                     small_breakpoints_copy[i] = small_breakpoints[i]

#                 elif (
#                     small_breakpoints[i - 1]["Row end #"] + 1
#                     != small_breakpoints[i]["Row start #"]
#                     and small_breakpoints[i - 1]["Row end #"] + 1
#                     == sub_density["Row start #"][j]
#                 ):
#                     small_breakpoints_copy[i] = {
#                         "Row start #": sub_density["Row start #"][j],
#                         "Row end #": sub_density["Row end #"][j],
#                         "main_bp": j,
#                     }
#                     for k in range(j, small_breakpoints[i]["main_bp"]):
#                         if (
#                             sub_density["Row end #"][k] + 1
#                             == small_breakpoints[i]["Row start #"]
#                         ):
#                             small_breakpoints_copy[i, k] = small_breakpoints[i]
#                             break
#                         elif (
#                             sub_density["Row end #"][k] + 1
#                             != small_breakpoints[i]["Row start #"]
#                         ):
#                             if (
#                                 sub_density["Row end #"][k] + 1
#                                 == sub_density["Row start #"][k + 1]
#                             ):
#                                 small_breakpoints_copy[i, k] = {
#                                     "Row start #": sub_density["Row start #"][k + 1],
#                                     "Row end #": sub_density["Row end #"][k + 1],
#                                     "main_bp": k + 1,
#                                 }
#     small_breakpoints_copy = {
#         i: small_breakpoints_copy[key] for i, key in enumerate(small_breakpoints_copy)
#     }
#     small_breakpoints = small_breakpoints_copy
    
#     return small_breakpoints

# def modify_small_breakpoints(small_breakpoints):
#     '''
#     This function modifies the small breakpoints.
    
#     Parameters
#     ----------
#     small_breakpoints : dict
#         Dictionary of small breakpoints.
        
#     Returns
#     -------
#     small_breakpoints : dict
#     Dictionary of modified small breakpoints.
#     '''
#     # update row start of current row start # == row end of previous row end # +1 in small_breakpoints
#     for i in range(len(small_breakpoints) - 1):
#         if i == 0:
#             small_breakpoints[i]["Row start #"] = small_breakpoints[i]["Row start #"]
#         elif i > 0 and i < len(small_breakpoints) - 1:
#             small_breakpoints[i]["Row start #"] = small_breakpoints[i -
#                                                                     1]["Row end #"] + 1
#         else:
#             small_breakpoints[i]["Row start #"] = small_breakpoints[i]["Row start #"]

#     return small_breakpoints

# def remove_invalid_entries(small_breakpoints):
#     '''
#     This function removes invalid entries from the small breakpoints.

#     Parameters
#     ----------
#     small_breakpoints : dict
#         Dictionary of small breakpoints.

#     Returns
#     -------
#     small_breakpoints_new : dict
#     Dictionary of small breakpoints with invalid entries removed.
#     '''
#     # remove entries from small_breakpoints that have ['Row end #'] < ['Row start #']
#     small_breakpoints_new = {}
#     for i in range(len(small_breakpoints)):
#         if small_breakpoints[i]["Row end #"] >= small_breakpoints[i]["Row start #"]:
#             small_breakpoints_new[i] = small_breakpoints[i]
#             # rename keys in small_breakpoints_new
#             small_breakpoints_new = {
#                 i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
#             }
#     small_breakpoints = small_breakpoints_new
#     return small_breakpoints



# small_BP = remove_invalid_entries(modify_small_breakpoints(merge_main_small_breakpoints(update_small_breakpoints(find_small_breakpoints(sort_gaps(merge_gaps(merge_gaps(get_gaps( "gene_matrix_output_Bra.xlsx", "subgenome_density_bra.xlsx", C, n_subgenomes,1)))), "subgenome_density_bra.xlsx")), "subgenome_density_bra.xlsx")))
# # print(small_BP)

# def update_df_synteny(C_df,synteny_file, df_subgenome_density_file_path):
#     '''
#     This function updates the df_synteny dataframe.

#     Parameters
#     ----------
#     df_synteny : dataframe
#         Dataframe of synteny.

#     Returns
#     -------
#     df_synteny : dataframe
#         Dataframe of synteny with chromosome names where the gene belongs.
#     '''
#         # Read the data
#     df_subgenome_density = pd.read_excel(df_subgenome_density_file_path)

#     # Get the column names from df_subgenome_density that start with N followed by a number
#     column_names = df_subgenome_density.filter(regex=r'^N\d+').columns.tolist()

#     # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
#     df_synteny = synteny_file.rename(columns={synteny_file.columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})

#     # Start the index from 0
#     df_synteny.index = df_synteny.index - 1
#     for i in range(len(df_synteny)):
#         for j in range(len(df_synteny.columns)):
#             if df_synteny.iloc[i, j] == 1:
#                 df_synteny.iloc[i, j] = df_synteny.columns[j]
#     #append the first column of the C_df to the df_synteny as first column
#     df_synteny.insert(0, "Gene_id", C_df.iloc[1:, 0])
    
#     df_synteny.to_excel("modified_synteny_chromosome_names.xlsx")

#     return df_synteny

# update_df_synteny(C_df, C_df_updated, "subgenome_density_bra.xlsx")

# def check_gaps(final_gaps_df):
#     for i in range(len(final_gaps_df) - 2):
#         if (final_gaps_df.iloc[i]["Row end #"] + 1 != final_gaps_df.iloc[i + 1]["Row start #"]):
#             # print("error1", i)
#             if (final_gaps_df.iloc[i-1]["Row end #"] + 1 == final_gaps_df.iloc[i + 1]["Row start #"]):
#                 final_gaps_df.drop(i, inplace=True)
#                 final_gaps_df.reset_index(drop=True, inplace=True)
#             else:
#                 final_gaps_df.iloc[i, final_gaps_df.columns.get_loc("Row end #")] = final_gaps_df.iloc[i + 1]["Row start #"] - 1
#         if final_gaps_df.iloc[i]["Row end #"] < final_gaps_df.iloc[i]["Row start #"]:
#             print("error2", i)
#     return final_gaps_df


# def subgenome_assignment_all_BP(sub_density, main_breakpoints, small_breakpoints):
#         sub_density = pd.read_excel(sub_density)
#         synteny_df = pd.read_excel(main_breakpoints, usecols="C:V")
#         # print(synteny_df)
#         final_gaps_df = pd.DataFrame(columns=["Row start #", "Row end #"])
#         sub1 = sub2 = sub3 = 0
#         for i in range(len(small_breakpoints)):
#             sub1_found = False
#             sub2_found = False
#             sub3_found = False

#             # final_gaps_df = final_gaps_df.append({'Row start #': small_breakpoints [i]['Row start #'], 'Row end #': small_breakpoints [i]['Row end #'], 'Subgenome1': sub_density['subgenome1'][small_breakpoints[i]['main_bp']], 'Subgenome2': sub_density['subgenome2'][small_breakpoints[i]['main_bp']], 'Subgenome3': sub_density['subgenome3'][small_breakpoints[i]['main_bp']]}, ignore_index=True)
#             # break out of outer loop as well
#             # break
#             if small_breakpoints[i]["Row start #"] != small_breakpoints[i]["Row end #"]:
#                 gaps_list_sub1 = []
#                 gaps_list_sub2 = []
#                 gaps_list_sub3 = []

#                 for l in range(
#                     int(small_breakpoints[i]["Row start #"]),
#                     int(small_breakpoints[i]["Row end #"] + 1),
#                 ):
#                     if (
#                         synteny_df.iloc[l][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome1"]
#                         ]
#                         != 0
#                         and synteny_df.iloc[l][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome1"]
#                         ]
#                         == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome1"]
#                     ):
#                         gaps_list_sub1.append(
#                             synteny_df.iloc[l][
#                                 sub_density.iloc[small_breakpoints[i]
#                                                 ["main_bp"]]["subgenome1"]
#                             ]
#                         )
#                     if (
#                         synteny_df.iloc[l][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome2"]
#                         ]
#                         != 0
#                         and synteny_df.iloc[l][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome2"]
#                         ]
#                         == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome2"]
#                     ):
#                         gaps_list_sub2.append(
#                             synteny_df.iloc[l][
#                                 sub_density.iloc[small_breakpoints[i]
#                                                 ["main_bp"]]["subgenome2"]
#                             ]
#                         )
#                     if (
#                         synteny_df.iloc[l][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome3"]
#                         ]
#                         != 0
#                         and synteny_df.iloc[l][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome3"]
#                         ]
#                         == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome3"]
#                     ):
#                         gaps_list_sub3.append(
#                             synteny_df.iloc[l][
#                                 sub_density.iloc[small_breakpoints[i]
#                                                 ["main_bp"]]["subgenome3"]
#                             ]
#                         )

#                 if len(gaps_list_sub1) > 0:
#                     sub1_found = False
#                     if (
#                         sub1_found == False
#                         and gaps_list_sub1[0]
#                         == sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
#                     ):
#                         sub1 = sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
#                         # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                         # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                         if sub2_found == False and sub3_found == False:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == False and sub3_found == True:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == True and sub3_found == False:
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         sub1_found = True
#                         # if sub2_found and sub3_found:
#                         #     break
#                         # break

#                 if len(gaps_list_sub1) == 0:
#                     sub1_found = False
#                     for l in range(
#                     int(small_breakpoints[i]["Row start #"]),
#                     int(small_breakpoints[i]["Row end #"] + 1),
#                 ):
#                         for b in range(len(synteny_df.columns)):
#                             if (
#                                 sub1_found == False
#                                 and synteny_df.iloc[l][b] != 0
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 == sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                             ):
#                                 if sub2_found == False and sub3_found == False:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1 = synteny_df.iloc[l][b]
#                                     sub1_found = True

#                                 if sub2_found == False and sub3_found == True:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub3 != synteny_df.iloc[l][b]:
#                                         sub1 = synteny_df.iloc[l][b]
#                                         sub1_found = True
#                                     else:
#                                         sub1 = sub_density["subgenome1"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub1_found = True
#                                 if sub2_found == True and sub3_found == False:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub2 != synteny_df.iloc[l][b]:
#                                         sub1 = synteny_df.iloc[l][b]
#                                         sub1_found = True
#                                     else:
#                                         sub1 = sub_density["subgenome1"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub1_found = True
#                                 if sub2_found == True and sub3_found == True:
#                                     if (
#                                         sub2 != synteny_df.iloc[l][b]
#                                         and sub3 != synteny_df.iloc[l][b]
#                                     ):
#                                         sub1 = synteny_df.iloc[l][b]
#                                         sub1_found = True
#                                     else:
#                                         sub1 = sub_density["subgenome1"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub1_found = True
#                                 # if sub2_found and sub3_found:
#                                 #     break
#                                 # break
#                         for b in range(len(synteny_df.columns)):
#                             if (
#                                 sub1_found == False
#                                 and synteny_df.iloc[l][b] != 0
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                             ):

#                                 if sub2_found == False and sub3_found == False:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1 = synteny_df.iloc[l][b]
#                                     sub1_found = True
#                                 if sub2_found == False and sub3_found == True:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub3 != synteny_df.iloc[l][b]:
#                                         sub1 = synteny_df.iloc[l][b]
#                                         sub1_found = True
#                                     else:
#                                         sub1 = sub_density["subgenome1"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub1_found = True
#                                 if sub2_found == True and sub3_found == False:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub2 != synteny_df.iloc[l][b]:
#                                         sub1 = synteny_df.iloc[l][b]
#                                         sub1_found = True
#                                     else:
#                                         sub1 = sub_density["subgenome1"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub1_found = True
#                                 if sub2_found == True and sub3_found == True:
#                                     if (
#                                         sub2 != synteny_df.iloc[l][b]
#                                         and sub3 != synteny_df.iloc[l][b]
#                                     ):
#                                         sub1 = synteny_df.iloc[l][b]
#                                         sub1_found = True
#                                     else:
#                                         sub1 = sub_density["subgenome1"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub1_found = True
#                                 # if sub2_found and sub3_found:
#                                 #     break
#                                 # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                                 # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                                 # break

#                 if len(gaps_list_sub2) > 0:
#                     sub2_found = False
#                     if (
#                         sub2_found == False
#                         and gaps_list_sub2[0]
#                         == sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
#                     ):
#                         # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                         sub2 = sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
#                         # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                         if sub1_found == False and sub3_found == False:
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub1_found == False and sub3_found == True:
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub1_found == True and sub3_found == False:
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         sub2_found = True
#                         # if sub1_found and sub3_found:
#                         #     break
#                         # break
                
#                 elif len(gaps_list_sub2) == 0:
#                     sub2_found = False
#                     for l in range(
#                     int(small_breakpoints[i]["Row start #"]),
#                     int(small_breakpoints[i]["Row end #"] + 1),
#                 ):
#                         for b in range(len(synteny_df.columns)):
#                             if (
#                                 sub2_found == False
#                                 and synteny_df.iloc[l][b] != 0
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 == sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                             ):
#                                 # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]

#                                 # if sub1_found and sub3_found:
#                                 #     break
#                                 if sub1_found == False and sub3_found == False:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2 = synteny_df.iloc[l][b]
#                                     sub2_found = True

#                                 if sub1_found == False and sub3_found == True:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub3 != synteny_df.iloc[l][b]:
#                                         sub2 = synteny_df.iloc[l][b]
#                                         sub2_found = True
#                                     else:
#                                         sub2 = sub_density["subgenome2"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub2_found = True

#                                 if sub1_found == True and sub3_found == False:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub1 != synteny_df.iloc[l][b]:
#                                         sub2 = synteny_df.iloc[l][b]
#                                         sub2_found = True
#                                     else:
#                                         sub2 = sub_density["subgenome2"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub2_found = True
#                                 if sub1_found == True and sub3_found == True:
#                                     if (
#                                         sub1 != synteny_df.iloc[l][b]
#                                         and sub3 != synteny_df.iloc[l][b]
#                                     ):
#                                         sub2 = synteny_df.iloc[l][b]
#                                         sub2_found = True
#                                     else:
#                                         sub2 = sub_density["subgenome2"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub2_found = True

#                                 # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                                 # break
#                         for b in range(len(synteny_df.columns)):
#                             if (
#                                 sub2_found == False
#                                 and synteny_df.iloc[l][b] != 0
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                             ):
#                                 # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                                 # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                                 if sub1_found == False and sub3_found == False:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2 = synteny_df.iloc[l][b]
#                                     sub2_found = True
#                                 if sub1_found == False and sub3_found == True:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub3 != synteny_df.iloc[l][b]:
#                                         sub2 = synteny_df.iloc[l][b]
#                                         sub2_found = True
#                                     else:
#                                         sub2 = sub_density["subgenome2"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub2_found = True
#                                 if sub1_found == True and sub3_found == False:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub1 != synteny_df.iloc[l][b]:
#                                         sub2 = synteny_df.iloc[l][b]
#                                         sub2_found = True
#                                     else:
#                                         sub2 = sub_density["subgenome2"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub2_found = True
#                                 if sub1_found == True and sub3_found == True:
#                                     if (
#                                         sub1 != synteny_df.iloc[l][b]
#                                         and sub3 != synteny_df.iloc[l][b]
#                                     ):
#                                         sub2 = synteny_df.iloc[l][b]
#                                         sub2_found = True
#                                     else:
#                                         sub2 = sub_density["subgenome2"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub2_found = True
#                                 # if sub1_found and sub3_found:
#                                 #     break
#                                 # break

#                 if len(gaps_list_sub3) > 0:
#                     sub3_found = False
#                     if (
#                         sub3_found == False
#                         and gaps_list_sub3[0]
#                         == sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
#                     ):
#                         # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                         # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                         sub3 = sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
#                         if sub2_found == False and sub1_found == False:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == False and sub1_found == True:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == True and sub1_found == False:
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         sub3_found = True
#                         # if sub1_found and sub2_found:
#                         #     break
#                         # break

#                 if len(gaps_list_sub3) == 0:
#                     sub3_found = False
#                     for l in range(
#                     int(small_breakpoints[i]["Row start #"]),
#                     int(small_breakpoints[i]["Row end #"] + 1),
#                 ):
#                         for b in range(len(synteny_df.columns)):
#                             if (
#                                 sub3_found == False
#                                 and synteny_df.iloc[l][b] != 0
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 == sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                             ):
#                                 # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                                 # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                                 if sub2_found == False and sub1_found == False:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3 = synteny_df.iloc[l][b]
#                                     sub3_found = True
#                                 if sub2_found == False and sub1_found == True:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub1 != synteny_df.iloc[l][b]:
#                                         sub3 = synteny_df.iloc[l][b]
#                                         sub3_found = True
#                                     else:
#                                         sub3 = sub_density["subgenome3"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub3_found = True
#                                 if sub2_found == True and sub1_found == False:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub2 != synteny_df.iloc[l][b]:
#                                         sub3 = synteny_df.iloc[l][b]
#                                         sub3_found = True
#                                     else:
#                                         sub3 = sub_density["subgenome3"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub3_found = True
#                                 if sub2_found == True and sub1_found == True:
#                                     if (
#                                         sub1 != synteny_df.iloc[l][b]
#                                         and sub2 != synteny_df.iloc[l][b]
#                                     ):
#                                         sub3 = synteny_df.iloc[l][b]
#                                         sub3_found = True
#                                     else:
#                                         sub3 = sub_density["subgenome3"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub3_found = True
#                                 # if sub1_found and sub2_found:
#                                 #     break
#                                 # break
#                         for b in range(len(synteny_df.columns)):
#                             if (
#                                 sub3_found == False
#                                 and synteny_df.iloc[l][b] != 0
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                                 and synteny_df.iloc[l][b].split(".")[0]
#                                 != sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ].split(".")[0]
#                             ):
#                                 # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                                 # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]

#                                 if sub2_found == False and sub1_found == False:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3 = synteny_df.iloc[l][b]
#                                     sub3_found = True
#                                 if sub2_found == False and sub1_found == True:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub1 != synteny_df.iloc[l][b]:
#                                         sub3 = synteny_df.iloc[l][b]
#                                         sub3_found = True
#                                     else:
#                                         sub3 = sub_density["subgenome3"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub3_found = True
#                                 if sub2_found == True and sub1_found == False:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     if sub2 != synteny_df.iloc[l][b]:
#                                         sub3 = synteny_df.iloc[l][b]
#                                         sub3_found = True
#                                     else:
#                                         sub3 = sub_density["subgenome3"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub3_found = True
#                                 if sub2_found == True and sub1_found == True:
#                                     if (
#                                         sub1 != synteny_df.iloc[l][b]
#                                         and sub2 != synteny_df.iloc[l][b]
#                                     ):
#                                         sub3 = synteny_df.iloc[l][b]
#                                         sub3_found = True
#                                     else:
#                                         sub3 = sub_density["subgenome3"][
#                                             small_breakpoints[i]["main_bp"]
#                                         ]
#                                         sub3_found = True
#                                 # if sub1_found and sub2_found:
#                                 #     break
#                                 # break
#             if small_breakpoints[i]["Row start #"] == small_breakpoints[i]["Row end #"]:
#                 gaps_list_sub1 = []
#                 gaps_list_sub2 = []
#                 gaps_list_sub3 = []

#                 if (
#                     synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                         sub_density.iloc[small_breakpoints[i]
#                                         ["main_bp"]]["subgenome1"]
#                     ]
#                     != 0
#                     and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                         sub_density.iloc[small_breakpoints[i]
#                                         ["main_bp"]]["subgenome1"]
#                     ]
#                     == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome1"]
#                 ):
#                     gaps_list_sub1.append(
#                         synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome1"]
#                         ]
#                     )
#                 if (
#                     synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                         sub_density.iloc[small_breakpoints[i]
#                                         ["main_bp"]]["subgenome2"]
#                     ]
#                     != 0
#                     and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                         sub_density.iloc[small_breakpoints[i]
#                                         ["main_bp"]]["subgenome2"]
#                     ]
#                     == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome2"]
#                 ):
#                     gaps_list_sub2.append(
#                         synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome2"]
#                         ]
#                     )
#                 if (
#                     synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                         sub_density.iloc[small_breakpoints[i]
#                                         ["main_bp"]]["subgenome3"]
#                     ]
#                     != 0
#                     and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                         sub_density.iloc[small_breakpoints[i]
#                                         ["main_bp"]]["subgenome3"]
#                     ]
#                     == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome3"]
#                 ):
#                     gaps_list_sub3.append(
#                         synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
#                             sub_density.iloc[small_breakpoints[i]
#                                             ["main_bp"]]["subgenome3"]
#                         ]
#                     )

#                 if len(gaps_list_sub1) > 0:
#                     sub1_found = False
#                     if (
#                         sub1_found == False
#                         and gaps_list_sub1[0]
#                         == sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
#                     ):
#                         sub1 = sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
#                         # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                         # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                         if sub2_found == False and sub3_found == False:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == False and sub3_found == True:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == True and sub3_found == False:
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         sub1_found = True
#                         # if sub2_found and sub3_found:
#                         #     break
#                         # break

#                 if len(gaps_list_sub1) == 0:
#                     sub1_found = False

#                     for b in range(len(synteny_df.columns)):
#                         if (
#                             sub1_found == False
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             == sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                         ):
#                             if sub2_found == False and sub3_found == False:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 sub1_found = True

#                             if sub2_found == False and sub3_found == True:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub1_found = True
#                                 else:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1_found = True
#                             if sub2_found == True and sub3_found == False:
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub1_found = True
#                                 else:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1_found = True
#                             if sub2_found == True and sub3_found == True:
#                                 if (
#                                     sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 ):
#                                     sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub1_found = True
#                                 else:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1_found = True
#                             # if sub2_found and sub3_found:
#                             #     break
#                             # break
#                     for b in range(len(synteny_df.columns)):
#                         if (
#                             sub1_found == False
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                         ):

#                             if sub2_found == False and sub3_found == False:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 sub1_found = True
#                             if sub2_found == False and sub3_found == True:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub1_found = True
#                                 else:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1_found = True
#                             if sub2_found == True and sub3_found == False:
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub1_found = True
#                                 else:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1_found = True
#                             if sub2_found == True and sub3_found == True:
#                                 if (
#                                     sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 ):
#                                     sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub1_found = True
#                                 else:
#                                     sub1 = sub_density["subgenome1"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub1_found = True
#                             # if sub2_found and sub3_found:
#                             #     break
#                             # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                             # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                             # break

#                 if len(gaps_list_sub2) > 0:
#                     sub2_found = False
#                     if (
#                         sub2_found == False
#                         and gaps_list_sub2[0]
#                         == sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
#                     ):
#                         # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                         sub2 = sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
#                         # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                         if sub1_found == False and sub3_found == False:
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub1_found == False and sub3_found == True:
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub1_found == True and sub3_found == False:
#                             sub3 = sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         sub2_found = True
#                         # if sub1_found and sub3_found:
#                         #     break
#                         # break
                
#                 elif len(gaps_list_sub2) == 0:
#                     sub2_found = False

#                     for b in range(len(synteny_df.columns)):
#                         if (
#                             sub2_found == False
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             == sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                         ):
#                             # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]

#                             # if sub1_found and sub3_found:
#                             #     break
#                             if sub1_found == False and sub3_found == False:
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 sub2_found = True

#                             if sub1_found == False and sub3_found == True:
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub2_found = True
#                                 else:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2_found = True

#                             if sub1_found == True and sub3_found == False:
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub2_found = True
#                                 else:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2_found = True
#                             if sub1_found == True and sub3_found == True:
#                                 if (
#                                     sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 ):
#                                     sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub2_found = True
#                                 else:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2_found = True

#                             # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                             # break
#                     for b in range(len(synteny_df.columns)):
#                         if (
#                             sub2_found == False
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                         ):
#                             # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                             # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
#                             if sub1_found == False and sub3_found == False:
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub2 = synteny_df.iloc[l][b]
#                                 sub2_found = True
#                             if sub1_found == False and sub3_found == True:
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub2_found = True
#                                 else:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2_found = True
#                             if sub1_found == True and sub3_found == False:
#                                 sub3 = sub_density["subgenome3"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub2_found = True
#                                 else:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2_found = True
#                             if sub1_found == True and sub3_found == True:
#                                 if (
#                                     sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 ):
#                                     sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub2_found = True
#                                 else:
#                                     sub2 = sub_density["subgenome2"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub2_found = True
#                             # if sub1_found and sub3_found:
#                             #     break
#                             # break

#                 if len(gaps_list_sub3) > 0:
#                     sub3_found = False
#                     if (
#                         sub3_found == False
#                         and gaps_list_sub3[0]
#                         == sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
#                     ):
#                         # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                         # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                         sub3 = sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
#                         if sub2_found == False and sub1_found == False:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == False and sub1_found == True:
#                             sub2 = sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         if sub2_found == True and sub1_found == False:
#                             sub1 = sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ]
#                         sub3_found = True
#                         # if sub1_found and sub2_found:
#                         #     break
#                         # break

#                 if len(gaps_list_sub3) == 0:
#                     sub3_found = False
#                     for b in range(len(synteny_df.columns)):
#                         if (
#                             sub3_found == False
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             == sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                         ):
#                             # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                             # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
#                             if sub2_found == False and sub1_found == False:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 sub3_found = True
#                             if sub2_found == False and sub1_found == True:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub3_found = True
#                                 else:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3_found = True
#                             if sub2_found == True and sub1_found == False:
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub3_found = True
#                                 else:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3_found = True
#                             if sub2_found == True and sub1_found == True:
#                                 if (
#                                     sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     and sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 ):
#                                     sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub3_found = True
#                                 else:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3_found = True
#                             # if sub1_found and sub2_found:
#                             #     break
#                             # break
#                     for b in range(len(synteny_df.columns)):
#                         if (
#                             sub3_found == False
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome3"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome2"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                             and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
#                             != sub_density["subgenome1"][
#                                 small_breakpoints[i]["main_bp"]
#                             ].split(".")[0]
#                         ):
#                             # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
#                             # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]

#                             if sub2_found == False and sub1_found == False:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 sub3_found = True
#                             if sub2_found == False and sub1_found == True:
#                                 sub2 = sub_density["subgenome2"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub3_found = True
#                                 else:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3_found = True
#                             if sub2_found == True and sub1_found == False:
#                                 sub1 = sub_density["subgenome1"][
#                                     small_breakpoints[i]["main_bp"]
#                                 ]
#                                 if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
#                                     sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub3_found = True
#                                 else:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3_found = True
#                             if sub2_found == True and sub1_found == True:
#                                 if (
#                                     sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     and sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                 ):
#                                     sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
#                                     sub3_found = True
#                                 else:
#                                     sub3 = sub_density["subgenome3"][
#                                         small_breakpoints[i]["main_bp"]
#                                     ]
#                                     sub3_found = True
#                             # if sub1_found and sub2_found:
#                             #     break
#                             # break
#             final_gaps_df = final_gaps_df.append(
#                 {
#                     "Row start #": small_breakpoints[i]["Row start #"],
#                     "Row end #": small_breakpoints[i]["Row end #"],
#                     "subgenome1": sub1,
#                     "subgenome2": sub2,
#                     "subgenome3": sub3,
#                 },
#                 ignore_index=True,
#             )

#                 # break
#         #save the final dataframe to excel
#         final_gaps_df = check_gaps(final_gaps_df)
#         final_gaps_df.to_excel("final_gaps_df.xlsx")
#         print("=========================================")
#         print("Number of small + main blocks:" , len(final_gaps_df))
#         print("=========================================")
              
#         return final_gaps_df

# final_BP = subgenome_assignment_all_BP("nodes.xlsx", "modified_synteny_chromosome_names.xlsx", small_BP)
# # print(final_BP)

# print("===========================================")
# print("Accuracy of subgenome assignment(Main BP + Small BP):")
# print("===========================================")
# subgenome_overlap(GT, "final_gaps_df.xlsx", df_synteny, 3)

