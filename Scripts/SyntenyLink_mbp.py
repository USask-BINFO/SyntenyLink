import numpy as np
import pandas as pd
import re
import warnings
import sys
import pickle
import csv
import os

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

# gap_threshold is the threshold to determine if there is a gap or not
# min_block_length is the minimum length of the signal blocks

def gap_calculation(C, gap_threshold, min_block_length, n, m):
    # get all gaps in all columns
    gaps_cell = [find_gaps(C[:, i], gap_threshold, m) for i in range(n)]
    break_point_indices = np.zeros(1000)
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

def create_excel_sheet(C_df, break_point_indices, densities, C_df_updated_copy):
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
    gene_id = C_df.iloc[1:,0]
    gene_id = gene_id.to_numpy()
    # prepare data for gene_start and gene_end columns for excel sheet
    gene_start = []
    gene_end = []
    b = len(break_point_indices)
    for i in range(b):
        gene_start.append(gene_id[row_start[i]])
        gene_end.append(gene_id[row_end[i]])

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
    prefixes = set((re.findall(r'^[A-Za-z]', col)[0], int(re.findall(r'\d+', col)[0])) for col in C_df_updated_copy.columns)
    prefixes = sorted(set([f'N{suffix}' if suffix % 2 == 1 else f'N{suffix}.r' for _, suffix in prefixes]))

    # Sort the columns by the letter prefix and the number following the letter
    cols = sorted(C_df_updated_copy.columns, key=lambda col: (re.findall(r'^[A-Za-z]', col)[0], int(re.findall(r'\d+', col)[0])))

    # Loop through the columns and update the header list
    for i in range(0, len(cols), 2):
        n = (i // 2) + 1
        header.extend([f"N{n}", f"N{n}.r"])

    # write matrix to excel sheet
    df = pd.DataFrame(E, columns=header)
    df.to_excel(f"Super_synteny_block_output.xlsx",
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

def assign_subgenomes(df_temp, file_path, n_subgenomes):
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
    for i in range(1, n_subgenomes+1):
        df.loc[:, f"subgenome{i}"] = df.iloc[:, :-i].apply(
            lambda x: x.sort_values(ascending=False).index[i-1], axis=1
        )

    # Read the generated break point output
    df_new = pd.read_excel(file_path, usecols="A:G")
    df = pd.concat([df_new.reset_index(drop=True), df], axis=1)

    # Assignment of the columns into subgenomes depending on gene density inside blocks
    for i in range(len(df)):
        if i == 0 and df["Non_zero"][i] == n_subgenomes - 1:
            if df["Non_zero"][i + 1] == n_subgenomes - 1:
                df.at[i, f"subgenome{n_subgenomes}"] = df.at[i + 2, f"subgenome{n_subgenomes}"]
            if df["Non_zero"][i + 1] > n_subgenomes - 1:
                df.at[i, f"subgenome{n_subgenomes}"] = df.at[i + 1, f"subgenome{n_subgenomes}"]
        elif df["Non_zero"][i] <= n_subgenomes - 1:
            if df["Non_zero"][i] == n_subgenomes - 1:
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
        
            if df["Non_zero"][i] == n_subgenomes - 2:
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
                
            if df["Non_zero"][i] == n_subgenomes - 3:
                df.at[i, f"subgenome{n_subgenomes-2}"] = df.at[i - 1, f"subgenome{n_subgenomes-2}"]
                df.at[i, f"subgenome{n_subgenomes-1}"] = df.at[i - 1, f"subgenome{n_subgenomes-1}"]
                df.at[i, f"subgenome{n_subgenomes}"] = df.at[i - 1, f"subgenome{n_subgenomes}"]
    return df
