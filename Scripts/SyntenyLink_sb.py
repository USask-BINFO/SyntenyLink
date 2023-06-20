import numpy as np
import pandas as pd
import re
import warnings
import sys
import pickle
import csv

def find_small_gaps(col, gap_threshold, m, row_start):
    """
    Finds the gaps in the each column of the matrix
    col: column of the matrix
    gap_threshold: threshold for the gap
    m: number of rows in the matrix

    returns: indices of the beginning and end of the gaps
    """
    # indices for ones in given column
    indices_one = np.where(col == 1)[0]

    # calculate difference between consective indices
    d = np.diff(indices_one)

    # if the difference is bigger than the gap threshold then store them
    initial_gaps = np.where(d > gap_threshold)[0]

    # number of gaps
    r = len(initial_gaps)

    # calculate gaps
    gaps = np.zeros((r, 2))

    # find the indices of beginning and end of gaps
    gaps[:, 0] = indices_one[initial_gaps] + 1 + row_start
    gaps[:, 1] = indices_one[initial_gaps + 1] - 1 + row_start

    # add row_start as the 1st column and indices_one[initial_gaps[0]] + 1 + row_start-1 as the second column for the first entry in gaps[]
    if r == 1:
        gaps[0, 0] = row_start
        gaps[0, 1] = (indices_one[initial_gaps[0]] + 1 + row_start) - 2
    if r > 1:
        # add without replacing the first entry
        gaps = np.insert(gaps, 0, np.array(
            [[row_start, gaps[0, 0] - 2]]), axis=0)
        gaps = np.insert(gaps, 1, np.array(
            [[gaps[1, 0] - 1, gaps[1, 0] - 1]]), axis=0)

    # fix the end of gap index for the last gap
    if r > 0 and indices_one[-1] + 1 < m:
        gaps = np.append(
            gaps,
            np.array([[indices_one[-1] + 1 + row_start, m + row_start - 1]]),
            axis=0,
        )
    if r > 0 and indices_one[-1] + 1 == m:
        gaps = np.append(
            gaps, np.array([[m + row_start - 1, m + row_start - 1]]), axis=0
        )
    return gaps

def get_gaps(breakpoints_file, sub_density_file, C, num_subgenomes, gap_threshold=1):
    '''
    Returns the gaps in the blocks
    
    breakpoints_file: file containing the main breakpoints
    sub_density_file: file containing the subgenome density and chromosome blocks
    num_subgenomes: number of subgenomes
    gap_threshold: threshold for the gap
    
    returns: list of gaps
    
    '''
    breakpoints_main = pd.read_excel(breakpoints_file, usecols="H:AA")
    sub_density = pd.read_excel(sub_density_file)
    gaps = []
    
    for i in range(len(breakpoints_main)):
        for j in range(len(breakpoints_main.columns)):
            if breakpoints_main.iloc[i, j] < 0.1 and breakpoints_main.iloc[i, j] != 0:
                # subgenome_names = [sub_density[f"subgenome{k+1}"][i].split(".")[0] for k in range(num_subgenomes)]
                # subgenome_combinations = [[subgenome_names[k], subgenome_names[l]] for k in range(num_subgenomes-1) for l in range(k+1, num_subgenomes)]
                if (
                (
                    sub_density["subgenome1"][i].split(".")[0]
                    != sub_density["subgenome2"][i]
                )
                and (
                    sub_density["subgenome1"][i].split(".")[0]
                    != sub_density["subgenome3"][i]
                )
                and (
                    sub_density["subgenome3"][i].split(".")[0]
                    != sub_density["subgenome2"][i]
                )
                and (
                    sub_density["subgenome2"][i].split(".")[0]
                    != sub_density["subgenome1"][i]
                )
                and (
                    sub_density["subgenome3"][i].split(".")[0]
                    != sub_density["subgenome1"][i]
                )
                and (
                    sub_density["subgenome2"][i].split(".")[0]
                    != sub_density["subgenome3"][i]
                )
            ):

                    # print("Block no.:", i+1, "Row no.:", sub_density['Row start #'][i], "Subgenome:", breakpoints_main.columns [j], "Density:", breakpoints_main.iloc[i,j])
                    # gap_threshold is the threshold to determine if there is a gap or not
                    # min_block_length is the minimum length of the signal blocks
                    gap_threshold = 1
                    # min_block_length = 10
                    # get all gaps in all columns
                    gaps_cell = [
                        find_small_gaps(
                            C[
                                sub_density["Row start #"][i]: sub_density["Row end #"][i]
                                + 1,
                                j,
                            ],
                            gap_threshold,
                            sub_density["# genes in block"][i],
                            sub_density["Row start #"][i],
                        )
                    ]
                    if len(gaps_cell[0]) != 0:
                        gaps.append((gaps_cell, i, j))
    return gaps

def merge_gaps(gaps):
    '''
    Merges gaps that are adjacent to each other.

    Parameters
    ----------
    gaps : list
        A list of tuples, where each tuple contains a list of arrays, the row number of the gap, and the column number of the gap.

    Returns
    -------
    merged_gaps : list
        A list of tuples, where each tuple contains a list of arrays, the row number of the gap, and the column number of the gap.
        
    '''
    merged_gaps = []
    i = 0
    while i < len(gaps) - 1:
        if gaps[i][1] == gaps[i + 1][1]:
            merged_row = []
            j = 0
            k = 0
            while j < len(gaps[i][0][0]) and k < len(gaps[i + 1][0][0]):
                if j == 0 and k == 0:
                    if gaps[i][0][0][0][1] > gaps[i + 1][0][0][0][1]:
                        merged_row.append(np.array(gaps[i + 1][0][0][k]))
                        k += 1

                    elif gaps[i][0][0][0][1] < gaps[i + 1][0][0][0][1]:
                        merged_row.append(np.array(gaps[i][0][0][j]))
                        j += 1

                    elif gaps[i][0][0][0][1] == gaps[i + 1][0][0][0][1]:
                        merged_row.append(np.array(gaps[i][0][0][j]))
                        j += 1
                        k += 1

                elif (
                    j > 0
                    and gaps[i + 1][0][0][k][0] == gaps[i][0][0][j][1] + 1
                    and k > 0
                ):
                    merged_row.append(np.array(gaps[i + 1][0][0][k]))
                    k += 1
                elif (
                    j > 0
                    and k > 0
                    and gaps[i + 1][0][0][k][0] != gaps[i][0][0][j][1] + 1
                    and gaps[i + 1][0][0][k][0] < gaps[i][0][0][j][0]
                    and gaps[i + 1][0][0][k][1] < gaps[i][0][0][j][1]
                ):
                    merged_row.append(np.array(gaps[i + 1][0][0][k]))
                    k += 1
                else:
                    merged_row.append(np.array(gaps[i][0][0][j]))
                    j += 1

            while j < len(gaps[i][0][0]):
                merged_row.append(np.array(gaps[i][0][0][j]))
                j += 1
            while k < len(gaps[i + 1][0][0]):
                merged_row.append(np.array(gaps[i + 1][0][0][k]))
                k += 1
            merged_gaps.append(([np.array(merged_row)], gaps[i][1]))
            i += 3
        elif (
            gaps[i][1] == gaps[i + 1][1]
            and gaps[i][1] != gaps[i + 2][1]
            and gaps[i][1] != gaps[i - 1][1]
        ):
            merged_row = []
            j = 0
            k = 0
            while j < len(gaps[i][0][0]) and k < len(gaps[i + 1][0][0]):
                if j == 0 and k == 0:
                    if gaps[i][0][0][0][1] > gaps[i + 1][0][0][0][1]:
                        merged_row.append(np.array(gaps[i + 1][0][0][k]))
                        k += 1

                    elif gaps[i][0][0][0][1] < gaps[i + 1][0][0][0][1]:
                        merged_row.append(np.array(gaps[i][0][0][j]))
                        j += 1

                    elif gaps[i][0][0][0][1] == gaps[i + 1][0][0][0][1]:
                        merged_row.append(np.array(gaps[i][0][0][j]))
                        j += 1
                        k += 1

                elif (
                    j > 0
                    and gaps[i + 1][0][0][k][0] == gaps[i][0][0][j][1] + 1
                    and k > 0
                ):
                    merged_row.append(np.array(gaps[i + 1][0][0][k]))
                    k += 1
                elif (
                    j > 0
                    and k > 0
                    and gaps[i + 1][0][0][k][0] != gaps[i][0][0][j][1] + 1
                    and gaps[i + 1][0][0][k][0] < gaps[i][0][0][j][0]
                    and gaps[i + 1][0][0][k][1] < gaps[i][0][0][j][1]
                ):
                    merged_row.append(np.array(gaps[i + 1][0][0][k]))
                    k += 1
                else:
                    merged_row.append(np.array(gaps[i][0][0][j]))
                    j += 1

            while j < len(gaps[i][0][0]):
                merged_row.append(np.array(gaps[i][0][0][j]))
                j += 1
            while k < len(gaps[i + 1][0][0]):
                merged_row.append(np.array(gaps[i + 1][0][0][k]))
                k += 1
            merged_gaps.append(
                ([np.array(merged_row)], gaps[i][1], gaps[i][2]))
            i += 2
        else:
            merged_gaps.append(gaps[i])
            i += 1
    if i == len(gaps) - 1:
        merged_gaps.append(gaps[i])
    return merged_gaps



def sort_gaps(gaps):
    '''
    Sort the gaps in each row of the gap list
    
    Parameters
    ----------
    gaps : list
        A list of tuples. Each tuple contains a list of arrays and a string.
        The list of arrays contains the gaps in each row of the gene matrix.
        The string is the name of the chromosome.
        
    Returns
    -------
    sorted_gaps : list
        A list of tuples. Each tuple contains a list of arrays and a string.
        The list of arrays contains the gaps in each row of the gene matrix.
        The string is the name of the chromosome.
        The gaps in each row are sorted by the start position of the gap.
    '''
    sorted_gaps = []
    for gap in gaps:
        sorted_gap = []
        for row in gap[0]:
            sorted_gap.append(row[row[:, 0].argsort()])
        sorted_gaps.append((sorted_gap, gap[1]))
    return sorted_gaps

def find_small_breakpoints(gaps, sub_density):
    '''
    Find the breakpoints that are smaller than the minimum gap size.
    
    Parameters
    ----------
    gaps : list
        A list of tuples. Each tuple contains a list of arrays and a string.
        The list of arrays contains the gaps in each row of the gene matrix.
        The string is the name of the chromosome.
        
    Returns
    -------
    small_breakpoints : dict
        A dictionary of dictionaries. The keys are the index of the gap
        and the index of the row in the gap. The values are dictionaries
        containing the start and end positions of the breakpoint and the
        main breakpoint.
    '''
    sub_density = pd.read_excel(sub_density)
    small_breakpoints = {}
    for i in range(len(gaps)):
        for j in range(len(gaps[i][0][0])):
            if len(gaps[i][0][0]) == 2:
                rowstart = gaps[i][0][0][j][0]
                row_end = gaps[i][0][0][j][1]
                if j == 0:
                    max_bp = sub_density["subgenome1"][gaps[i][1]]
                    small_breakpoints[i, j] = {
                        "Row start #": rowstart,
                        "Row end #": row_end,
                        "main_bp": gaps[i][1],
                    }

                if j == 1 and rowstart > gaps[i][0][0][0][1] + 1:
                    rowstart = gaps[i][0][0][j - 1][1] + 1
                    row_end = gaps[i][0][0][j][0] - 1

                    small_breakpoints[i, j] = {
                        "Row start #": rowstart,
                        "Row end #": row_end,
                        "main_bp": gaps[i][1],
                    }

                if j == len(gaps[i][0][0]) - 1:
                    rowstart = gaps[i][0][0][j][0]
                    row_end = gaps[i][0][0][j][1]
                    small_breakpoints[i, j + 1] = {
                        "Row start #": rowstart,
                        "Row end #": row_end,
                        "main_bp": gaps[i][1],
                    }

            if len(gaps[i][0][0]) > 2:
                rowstart = gaps[i][0][0][j][0]
                row_end = gaps[i][0][0][j][1]
                if j == 0:
                    small_breakpoints[i, j] = {
                        "Row start #": rowstart,
                        "Row end #": row_end,
                        "main_bp": gaps[i][1],
                    }
                if j == 1 and rowstart == gaps[i][0][0][0][1] + 1:
                    rowstart = gaps[i][0][0][j][1]
                    row_end = gaps[i][0][0][j][0]
                    small_breakpoints[i, j] = {
                        "Row start #": rowstart,
                        "Row end #": row_end,
                        "main_bp": gaps[i][1],
                    }

                if (j <= (len(gaps[i][0][0]) - 1)) and (j >= 1):
                    if gaps[i][0][0][j][0] > gaps[i][0][0][j - 1][1] + 1:
                        rowstart = gaps[i][0][0][j - 1][1] + 1
                        row_end = gaps[i][0][0][j][0] - 1
                        small_breakpoints[i, j] = {
                            "Row start #": rowstart,
                            "Row end #": row_end,
                            "main_bp": gaps[i][1],
                        }

                    rowstart = gaps[i][0][0][j][0]
                    row_end = gaps[i][0][0][j][1]
                    small_breakpoints[i, j, 0] = {
                        "Row start #": rowstart,
                        "Row end #": row_end,
                        "main_bp": gaps[i][1],
                    }

    return small_breakpoints


def update_small_breakpoints(small_breakpoints):
    '''
    Update the small breakpoints.
    
    Parameters
    ----------
    small_breakpoints : dict
        A dictionary of dictionaries. The keys are the index of the gap
        and the index of the row in the gap. The values are dictionaries
        containing the start and end positions of the breakpoint and the
        main breakpoint.
        
    Returns
    -------
    small_breakpoints_new : dict
        A dictionary of dictionaries. The keys are the index of the gap
        and the index of the row in the gap. The values are dictionaries
        containing the start and end positions of the breakpoint and the
        main breakpoint.
    '''
    small_breakpoints = {
        i: small_breakpoints[key] for i, key in enumerate(small_breakpoints)
    }
    small_breakpoints
    # Check if the row start and row end numbers are like row start number = from the previous row end number + 1 and row end number = from the next row start number - 1, if not, update the row start and row end numbers.
    small_breakpoints_new = small_breakpoints.copy()
    for i in range(len(small_breakpoints) - 1):
        if i == 0:
            small_breakpoints_new[i] = small_breakpoints[i]
        if i == len(small_breakpoints) - 1:
            small_breakpoints_new[i] = small_breakpoints[i]
        if i > 0 and i < len(small_breakpoints) - 1:
            # if small_breakpoints[i]['Row end #'] == small_breakpoints[i+1]['Row end #'] and small_breakpoints[i]['Row start #'] == small_breakpoints[i+1]['Row start #']:
            #     continue
            if (
                small_breakpoints[i]["Row end #"]
                != small_breakpoints[i + 1]["Row start #"] - 1
            ):
                # Check whether you can find this row end number in any of the row end numbers coming in upcoming rows.
                if (
                    small_breakpoints[i]["Row start #"]
                    != small_breakpoints[i + 1]["Row start #"]
                    and small_breakpoints[i]["Row end #"]
                    != small_breakpoints[i + 1]["Row end #"]
                ):
                    row_end_list = []
                    for s in range(i + 1, len(small_breakpoints)):
                        row_end_list.append(small_breakpoints[s]["Row end #"])
                    if small_breakpoints[i]["Row end #"] in row_end_list:
                        # If you can find this row end number in any of the row end numbers coming in upcoming rows, then update the current row's end number.
                        small_breakpoints_new[i] = {
                            "Row start #": small_breakpoints[i]["Row start #"],
                            "Row end #": small_breakpoints[i + 1]["Row start #"] - 1,
                            "main_bp": small_breakpoints[i]["main_bp"],
                        }

                    elif (
                        small_breakpoints[i]["Row end #"] not in row_end_list
                        and small_breakpoints[i]["Row end #"]
                        > small_breakpoints[i + 1]["Row start #"]
                    ):
                        small_breakpoints_new[i] = {
                            "Row start #": small_breakpoints[i]["Row start #"],
                            "Row end #": small_breakpoints[i + 1]["Row start #"] - 1,
                            "main_bp": small_breakpoints[i]["main_bp"],
                        }
                        # If you can't find this row end number in any of the row end numbers coming in upcoming rows, then add it where row end number < row end number of next row but row end number > row end number of the row of focus and row start number of next row < row end number.
                        for j in range(i + 1, len(small_breakpoints) - 1):
                            if (
                                small_breakpoints[i]["Row end #"]
                                < small_breakpoints[j + 1]["Row start #"]
                                and small_breakpoints[i]["Row end #"]
                                > small_breakpoints[j]["Row start #"]
                            ):
                                # Add this row end number and row start number of next row as a new row.
                                new_row = {
                                    "Row start #": small_breakpoints[j]["Row start #"],
                                    "Row end #": small_breakpoints[i]["Row end #"],
                                    "main_bp": small_breakpoints[i]["main_bp"],
                                }
                                small_breakpoints_new[i + 1, j] = new_row

                if (
                    small_breakpoints[i]["Row start #"]
                    == small_breakpoints[i + 1]["Row start #"]
                    and small_breakpoints[i]["Row end #"]
                    != small_breakpoints[i + 1]["Row end #"]
                    and small_breakpoints[i]["Row end #"]
                    <= small_breakpoints[i + 1]["Row end #"]
                ):
                    small_breakpoints_new[i] = small_breakpoints[i]
                    small_breakpoints_new[i + 1] = {
                        "Row start #": small_breakpoints[i]["Row end #"] + 1,
                        "Row end #": small_breakpoints[i + 1]["Row end #"],
                        "main_bp": small_breakpoints[i + 1]["main_bp"],
                    }
                else:
                    continue

            if (
                small_breakpoints[i]["Row end #"]
                == small_breakpoints[i + 1]["Row start #"] - 1
            ):
                if i < len(small_breakpoints) - 2:
                    if (
                        small_breakpoints[i + 1]["Row start #"]
                        == small_breakpoints[i + 2]["Row start #"]
                        and small_breakpoints[i + 1]["Row end #"]
                        > small_breakpoints[i + 2]["Row start #"]
                    ):
                        small_breakpoints_new[i] = small_breakpoints[i]
                        small_breakpoints_new[i + 1] = small_breakpoints[i + 2]

    small_breakpoints = {
        i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
    }

    return small_breakpoints

def merge_main_small_breakpoints(small_breakpoints_new, sub_density_file):
    '''
    This function merges the main breakpoints and small breakpoints.
    
    Parameters
    ----------
    small_breakpoints_new : dict
        Dictionary of small breakpoints.
        sub_density_file : str
        Path to the sub density file.
        
    Returns
    -------
    small_breakpoints_copy : dict
    Dictionary of merged breakpoints.
    '''
    sub_density = pd.read_excel(sub_density_file)
    small_breakpoints = {
        i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
    }

    sub_density = (
        sub_density.copy()
    )  # make a copy to avoid modifying the original DataFrame
    small_breakpoints_copy = {}  # make a copy to avoid modifying the original dictionary

    for i in range(len(small_breakpoints) - 1):
        if i == 0:
            small_breakpoints_copy[i] = small_breakpoints[i]
        if i == len(small_breakpoints) - 1:
            small_breakpoints_copy[i] = small_breakpoints[i]
        for j in range(len(sub_density)):
            if sub_density["Row start #"][j] < small_breakpoints[0]["Row start #"]:
                small_breakpoints_copy[j] = {
                    "Row start #": sub_density["Row start #"][j],
                    "Row end #": sub_density["Row end #"][j],
                    "main_bp": j,
                }
            if sub_density["Row start #"][j] > small_breakpoints[0]["Row start #"]:
                if i > 0 and i < len(small_breakpoints) - 1:
                    if (
                        small_breakpoints[i - 1]["Row end #"] + 1
                        == small_breakpoints[i]["Row start #"]
                    ):
                        small_breakpoints_copy[i] = small_breakpoints[i]
                    if (
                        small_breakpoints[i - 1]["Row end #"] + 1
                        != small_breakpoints[i]["Row start #"]
                        and small_breakpoints[i - 1]["Row start #"]
                        == sub_density["Row start #"][j]
                    ):
                        small_breakpoints_copy[i] = small_breakpoints[i]

                    elif (
                        small_breakpoints[i - 1]["Row end #"] + 1
                        != small_breakpoints[i]["Row start #"]
                        and small_breakpoints[i - 1]["Row end #"] + 1
                        == sub_density["Row start #"][j]
                    ):
                        small_breakpoints_copy[i] = {
                            "Row start #": sub_density["Row start #"][j],
                            "Row end #": sub_density["Row end #"][j],
                            "main_bp": j,
                        }
                        for k in range(j, small_breakpoints[i]["main_bp"]):
                            if (
                                sub_density["Row end #"][k] + 1
                                == small_breakpoints[i]["Row start #"]
                            ):
                                small_breakpoints_copy[i, k] = small_breakpoints[i]
                                break
                            elif (
                                sub_density["Row end #"][k] + 1
                                != small_breakpoints[i]["Row start #"]
                            ):
                                if (
                                    sub_density["Row end #"][k] + 1
                                    == sub_density["Row start #"][k + 1]
                                ):
                                    small_breakpoints_copy[i, k] = {
                                        "Row start #": sub_density["Row start #"][k + 1],
                                        "Row end #": sub_density["Row end #"][k + 1],
                                        "main_bp": k + 1,
                                    }
    last_small_bp = list(small_breakpoints_copy.keys())[-1]
    last_sub_density_row_end = sub_density.iloc[-1]['Row end #']
    if small_breakpoints_copy[last_small_bp]['Row end #'] < last_sub_density_row_end:
        for j in range(len(sub_density)):
            if sub_density['Row start #'][j] > small_breakpoints_copy[last_small_bp]['Row end #']:
                small_breakpoints_copy[last_small_bp+1+j] = {
                    'Row start #': sub_density['Row start #'][j],
                    'Row end #': sub_density['Row end #'][j],
                    'main_bp': j
                }

    small_breakpoints_copy = {
        i: small_breakpoints_copy[key] for i, key in enumerate(small_breakpoints_copy)
    }
    small_breakpoints = small_breakpoints_copy
    return small_breakpoints


def modify_small_breakpoints(small_breakpoints):
    '''
    This function modifies the small breakpoints.
    
    Parameters
    ----------
    small_breakpoints : dict
        Dictionary of small breakpoints.
        
    Returns
    -------
    small_breakpoints : dict
    Dictionary of modified small breakpoints.
    '''
    # update row start of current row start # == row end of previous row end # +1 in small_breakpoints
    for i in range(len(small_breakpoints) - 1):
        if i == 0:
            small_breakpoints[i]["Row start #"] = small_breakpoints[i]["Row start #"]
        elif i > 0 and i < len(small_breakpoints) - 1:
            small_breakpoints[i]["Row start #"] = small_breakpoints[i -
                                                                    1]["Row end #"] + 1
            small_breakpoints[i]["Row end #"] = small_breakpoints[i +
                                                                    1]["Row start #"] -1
        else:
            small_breakpoints[i]["Row start #"] = small_breakpoints[i]["Row start #"]

    return small_breakpoints

def update_df_synteny(C_df,synteny_file, df_subgenome_density_file_path, modified_synteny):
    '''
    This function updates the df_synteny dataframe.

    Parameters
    ----------
    df_synteny : dataframe
        Dataframe of synteny.

    Returns
    -------
    df_synteny : dataframe
        Dataframe of synteny with chromosome names where the gene belongs.
    '''
        # Read the data
    df_subgenome_density = pd.read_excel(df_subgenome_density_file_path)

    # Get the column names from df_subgenome_density that start with N followed by a number
    column_names = df_subgenome_density.filter(regex=r'^N\d+').columns.tolist()

    # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
    df_synteny = synteny_file.rename(columns={synteny_file.columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})

    # Start the index from 0
    df_synteny.index = df_synteny.index - 1
    for i in range(len(df_synteny)):
        for j in range(len(df_synteny.columns)):
            if df_synteny.iloc[i, j] == 1:
                df_synteny.iloc[i, j] = df_synteny.columns[j]
    #append the first column of the C_df to the df_synteny as first column
    df_synteny.insert(0, "Gene_id", C_df.iloc[1:, 0])
    
    df_synteny.to_excel(modified_synteny)

    return df_synteny


def remove_invalid_entries(small_breakpoints):
    '''
    This function removes invalid entries from the small breakpoints.

    Parameters
    ----------
    small_breakpoints : dict
        Dictionary of small breakpoints.

    Returns
    -------
    small_breakpoints_new : dict
        Dictionary of small breakpoints with invalid entries removed.
    '''
    small_breakpoints_new = {}
    for i in range(len(small_breakpoints)):
        if small_breakpoints[i]["Row end #"] >= small_breakpoints[i]["Row start #"]:
            small_breakpoints_new[i] = small_breakpoints[i]
    
    i = 1
    while i < len(small_breakpoints_new):
        small_breakpoints_new = {
        i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
    }
        if small_breakpoints_new[i-1]["Row end #"] == small_breakpoints_new[i]["Row start #"]:
            del small_breakpoints_new[i-1]
        else:
            i += 1

    small_breakpoints_new = {
        i: small_breakpoints_new[key] for i, key in enumerate(small_breakpoints_new)
    }
    
    return small_breakpoints_new

# small_BP = remove_invalid_entries(modify_small_breakpoints(merge_main_small_breakpoints(update_small_breakpoints(find_small_breakpoints(sort_gaps(merge_gaps(merge_gaps(get_gaps( "gene_matrix_output_Bra.xlsx", "subgenome_density_bra.xlsx", C, n_subgenomes,1)))), "subgenome_density_bra.xlsx")), "subgenome_density_bra.xlsx")))
# print(small_BP)

def check_gaps(final_gaps_df):
    for i in range(len(final_gaps_df)- 1):
        if i < len(final_gaps_df)- 1:
            if (final_gaps_df.iloc[i]["Row end #"] + 1 != final_gaps_df.iloc[i + 1]["Row start #"]):
                # print("error1", i)
                if (final_gaps_df.iloc[i-1]["Row end #"] + 1 == final_gaps_df.iloc[i + 1]["Row start #"]):
                    final_gaps_df.drop(i, inplace=True)
                    final_gaps_df.reset_index(drop=True, inplace=True)
                else:
                    final_gaps_df.at[i, "Row end #"] = final_gaps_df.iloc[i + 1]["Row start #"] - 1
            if final_gaps_df.iloc[i]["Row end #"] < final_gaps_df.iloc[i]["Row start #"]:
                final_gaps_df.at[i-1, "Row end #"] = final_gaps_df.iloc[i-1]["Row start #"] + 1
                final_gaps_df.at[i, "Row start #"] = final_gaps_df.iloc[i-1]["Row end #"] + 1
            if final_gaps_df.iloc[i]["Row end #"] < final_gaps_df.iloc[i]["Row start #"]:
                print("error2", i)
    return final_gaps_df


def subgenome_assignment_all_BP(sub_density, main_breakpoints, small_breakpoints, final_gaps_file):
        sub_density = pd.read_excel(sub_density)
        synteny_df = pd.read_excel(main_breakpoints, usecols="C:V")
        # print(synteny_df)
        final_gaps_df = pd.DataFrame(columns=["Row start #", "Row end #"])
        sub1 = sub2 = sub3 = 0
        for i in range(len(small_breakpoints)):
            sub1_found = False
            sub2_found = False
            sub3_found = False
            if i == 0:
                sub1 = sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome1"]
                sub2 = sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome2"]
                sub3 = sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome3"]
            # final_gaps_df = final_gaps_df.append({'Row start #': small_breakpoints [i]['Row start #'], 'Row end #': small_breakpoints [i]['Row end #'], 'Subgenome1': sub_density['subgenome1'][small_breakpoints[i]['main_bp']], 'Subgenome2': sub_density['subgenome2'][small_breakpoints[i]['main_bp']], 'Subgenome3': sub_density['subgenome3'][small_breakpoints[i]['main_bp']]}, ignore_index=True)
            # break out of outer loop as well
            # break
            if small_breakpoints[i]["Row start #"] != small_breakpoints[i]["Row end #"]:
                gaps_list_sub1 = []
                gaps_list_sub2 = []
                gaps_list_sub3 = []

                for l in range(
                    int(small_breakpoints[i]["Row start #"]),
                    int(small_breakpoints[i]["Row end #"] + 1),
                ):
                    if (
                        synteny_df.iloc[l][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome1"]
                        ]
                        != 0
                        and synteny_df.iloc[l][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome1"]
                        ]
                        == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome1"]
                    ):
                        gaps_list_sub1.append(
                            synteny_df.iloc[l][
                                sub_density.iloc[small_breakpoints[i]
                                                ["main_bp"]]["subgenome1"]
                            ]
                        )
                    if (
                        synteny_df.iloc[l][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome2"]
                        ]
                        != 0
                        and synteny_df.iloc[l][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome2"]
                        ]
                        == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome2"]
                    ):
                        gaps_list_sub2.append(
                            synteny_df.iloc[l][
                                sub_density.iloc[small_breakpoints[i]
                                                ["main_bp"]]["subgenome2"]
                            ]
                        )
                    if (
                        synteny_df.iloc[l][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome3"]
                        ]
                        != 0
                        and synteny_df.iloc[l][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome3"]
                        ]
                        == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome3"]
                    ):
                        gaps_list_sub3.append(
                            synteny_df.iloc[l][
                                sub_density.iloc[small_breakpoints[i]
                                                ["main_bp"]]["subgenome3"]
                            ]
                        )

                if len(gaps_list_sub1) > 0:
                    sub1_found = False
                    if (
                        sub1_found == False
                        and gaps_list_sub1[0]
                        == sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
                    ):
                        sub1 = sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
                        # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                        # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                        if sub2_found == False and sub3_found == False:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == False and sub3_found == True:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == True and sub3_found == False:
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        sub1_found = True
                        # if sub2_found and sub3_found:
                        #     break
                        # break

                if len(gaps_list_sub1) == 0:
                    sub1_found = False
                    for l in range(
                    int(small_breakpoints[i]["Row start #"]),
                    int(small_breakpoints[i]["Row end #"] + 1),
                ):
                        for b in range(len(synteny_df.columns)):
                            if (
                                sub1_found == False
                                and synteny_df.iloc[l][b] != 0
                                and synteny_df.iloc[l][b].split(".")[0]
                                == sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                            ):
                                if sub2_found == False and sub3_found == False:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1 = synteny_df.iloc[l][b]
                                    sub1_found = True

                                if sub2_found == False and sub3_found == True:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub3 != synteny_df.iloc[l][b]:
                                        sub1 = synteny_df.iloc[l][b]
                                        sub1_found = True
                                    else:
                                        sub1 = sub_density["subgenome1"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub1_found = True
                                if sub2_found == True and sub3_found == False:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub2 != synteny_df.iloc[l][b]:
                                        sub1 = synteny_df.iloc[l][b]
                                        sub1_found = True
                                    else:
                                        sub1 = sub_density["subgenome1"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub1_found = True
                                if sub2_found == True and sub3_found == True:
                                    if (
                                        sub2 != synteny_df.iloc[l][b]
                                        and sub3 != synteny_df.iloc[l][b]
                                    ):
                                        sub1 = synteny_df.iloc[l][b]
                                        sub1_found = True
                                    else:
                                        sub1 = sub_density["subgenome1"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub1_found = True
                                # if sub2_found and sub3_found:
                                #     break
                                # break
                        for b in range(len(synteny_df.columns)):
                            if (
                                sub1_found == False
                                and synteny_df.iloc[l][b] != 0
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                            ):

                                if sub2_found == False and sub3_found == False:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1 = synteny_df.iloc[l][b]
                                    sub1_found = True
                                if sub2_found == False and sub3_found == True:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub3 != synteny_df.iloc[l][b]:
                                        sub1 = synteny_df.iloc[l][b]
                                        sub1_found = True
                                    else:
                                        sub1 = sub_density["subgenome1"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub1_found = True
                                if sub2_found == True and sub3_found == False:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub2 != synteny_df.iloc[l][b]:
                                        sub1 = synteny_df.iloc[l][b]
                                        sub1_found = True
                                    else:
                                        sub1 = sub_density["subgenome1"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub1_found = True
                                if sub2_found == True and sub3_found == True:
                                    if (
                                        sub2 != synteny_df.iloc[l][b]
                                        and sub3 != synteny_df.iloc[l][b]
                                    ):
                                        sub1 = synteny_df.iloc[l][b]
                                        sub1_found = True
                                    else:
                                        sub1 = sub_density["subgenome1"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub1_found = True
                                # if sub2_found and sub3_found:
                                #     break
                                # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                                # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                                # break

                if len(gaps_list_sub2) > 0:
                    sub2_found = False
                    if (
                        sub2_found == False
                        and gaps_list_sub2[0]
                        == sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
                    ):
                        # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                        sub2 = sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
                        # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                        if sub1_found == False and sub3_found == False:
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub1_found == False and sub3_found == True:
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub1_found == True and sub3_found == False:
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        sub2_found = True
                        # if sub1_found and sub3_found:
                        #     break
                        # break
                
                elif len(gaps_list_sub2) == 0:
                    sub2_found = False
                    for l in range(
                    int(small_breakpoints[i]["Row start #"]),
                    int(small_breakpoints[i]["Row end #"] + 1),
                ):
                        for b in range(len(synteny_df.columns)):
                            if (
                                sub2_found == False
                                and synteny_df.iloc[l][b] != 0
                                and synteny_df.iloc[l][b].split(".")[0]
                                == sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                            ):
                                # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]

                                # if sub1_found and sub3_found:
                                #     break
                                if sub1_found == False and sub3_found == False:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2 = synteny_df.iloc[l][b]
                                    sub2_found = True

                                if sub1_found == False and sub3_found == True:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub3 != synteny_df.iloc[l][b]:
                                        sub2 = synteny_df.iloc[l][b]
                                        sub2_found = True
                                    else:
                                        sub2 = sub_density["subgenome2"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub2_found = True

                                if sub1_found == True and sub3_found == False:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub1 != synteny_df.iloc[l][b]:
                                        sub2 = synteny_df.iloc[l][b]
                                        sub2_found = True
                                    else:
                                        sub2 = sub_density["subgenome2"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub2_found = True
                                if sub1_found == True and sub3_found == True:
                                    if (
                                        sub1 != synteny_df.iloc[l][b]
                                        and sub3 != synteny_df.iloc[l][b]
                                    ):
                                        sub2 = synteny_df.iloc[l][b]
                                        sub2_found = True
                                    else:
                                        sub2 = sub_density["subgenome2"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub2_found = True

                                # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                                # break
                        for b in range(len(synteny_df.columns)):
                            if (
                                sub2_found == False
                                and synteny_df.iloc[l][b] != 0
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                            ):
                                # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                                # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                                if sub1_found == False and sub3_found == False:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2 = synteny_df.iloc[l][b]
                                    sub2_found = True
                                if sub1_found == False and sub3_found == True:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub3 != synteny_df.iloc[l][b]:
                                        sub2 = synteny_df.iloc[l][b]
                                        sub2_found = True
                                    else:
                                        sub2 = sub_density["subgenome2"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub2_found = True
                                if sub1_found == True and sub3_found == False:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub1 != synteny_df.iloc[l][b]:
                                        sub2 = synteny_df.iloc[l][b]
                                        sub2_found = True
                                    else:
                                        sub2 = sub_density["subgenome2"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub2_found = True
                                if sub1_found == True and sub3_found == True:
                                    if (
                                        sub1 != synteny_df.iloc[l][b]
                                        and sub3 != synteny_df.iloc[l][b]
                                    ):
                                        sub2 = synteny_df.iloc[l][b]
                                        sub2_found = True
                                    else:
                                        sub2 = sub_density["subgenome2"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub2_found = True
                                # if sub1_found and sub3_found:
                                #     break
                                # break

                if len(gaps_list_sub3) > 0:
                    sub3_found = False
                    if (
                        sub3_found == False
                        and gaps_list_sub3[0]
                        == sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
                    ):
                        # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                        # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                        sub3 = sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
                        if sub2_found == False and sub1_found == False:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == False and sub1_found == True:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == True and sub1_found == False:
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        sub3_found = True
                        # if sub1_found and sub2_found:
                        #     break
                        # break

                if len(gaps_list_sub3) == 0:
                    sub3_found = False
                    for l in range(
                    int(small_breakpoints[i]["Row start #"]),
                    int(small_breakpoints[i]["Row end #"] + 1),
                ):
                        for b in range(len(synteny_df.columns)):
                            if (
                                sub3_found == False
                                and synteny_df.iloc[l][b] != 0
                                and synteny_df.iloc[l][b].split(".")[0]
                                == sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                            ):
                                # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                                # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                                if sub2_found == False and sub1_found == False:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3 = synteny_df.iloc[l][b]
                                    sub3_found = True
                                if sub2_found == False and sub1_found == True:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub1 != synteny_df.iloc[l][b]:
                                        sub3 = synteny_df.iloc[l][b]
                                        sub3_found = True
                                    else:
                                        sub3 = sub_density["subgenome3"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub3_found = True
                                if sub2_found == True and sub1_found == False:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub2 != synteny_df.iloc[l][b]:
                                        sub3 = synteny_df.iloc[l][b]
                                        sub3_found = True
                                    else:
                                        sub3 = sub_density["subgenome3"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub3_found = True
                                if sub2_found == True and sub1_found == True:
                                    if (
                                        sub1 != synteny_df.iloc[l][b]
                                        and sub2 != synteny_df.iloc[l][b]
                                    ):
                                        sub3 = synteny_df.iloc[l][b]
                                        sub3_found = True
                                    else:
                                        sub3 = sub_density["subgenome3"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub3_found = True
                                # if sub1_found and sub2_found:
                                #     break
                                # break
                        for b in range(len(synteny_df.columns)):
                            if (
                                sub3_found == False
                                and synteny_df.iloc[l][b] != 0
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                                and synteny_df.iloc[l][b].split(".")[0]
                                != sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ].split(".")[0]
                            ):
                                # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                                # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]

                                if sub2_found == False and sub1_found == False:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3 = synteny_df.iloc[l][b]
                                    sub3_found = True
                                if sub2_found == False and sub1_found == True:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub1 != synteny_df.iloc[l][b]:
                                        sub3 = synteny_df.iloc[l][b]
                                        sub3_found = True
                                    else:
                                        sub3 = sub_density["subgenome3"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub3_found = True
                                if sub2_found == True and sub1_found == False:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    if sub2 != synteny_df.iloc[l][b]:
                                        sub3 = synteny_df.iloc[l][b]
                                        sub3_found = True
                                    else:
                                        sub3 = sub_density["subgenome3"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub3_found = True
                                if sub2_found == True and sub1_found == True:
                                    if (
                                        sub1 != synteny_df.iloc[l][b]
                                        and sub2 != synteny_df.iloc[l][b]
                                    ):
                                        sub3 = synteny_df.iloc[l][b]
                                        sub3_found = True
                                    else:
                                        sub3 = sub_density["subgenome3"][
                                            small_breakpoints[i]["main_bp"]
                                        ]
                                        sub3_found = True
                                # if sub1_found and sub2_found:
                                #     break
                                # break
            if small_breakpoints[i]["Row start #"] == small_breakpoints[i]["Row end #"]:
                gaps_list_sub1 = []
                gaps_list_sub2 = []
                gaps_list_sub3 = []

                if (
                    synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                        sub_density.iloc[small_breakpoints[i]
                                        ["main_bp"]]["subgenome1"]
                    ]
                    != 0
                    and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                        sub_density.iloc[small_breakpoints[i]
                                        ["main_bp"]]["subgenome1"]
                    ]
                    == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome1"]
                ):
                    gaps_list_sub1.append(
                        synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome1"]
                        ]
                    )
                if (
                    synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                        sub_density.iloc[small_breakpoints[i]
                                        ["main_bp"]]["subgenome2"]
                    ]
                    != 0
                    and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                        sub_density.iloc[small_breakpoints[i]
                                        ["main_bp"]]["subgenome2"]
                    ]
                    == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome2"]
                ):
                    gaps_list_sub2.append(
                        synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome2"]
                        ]
                    )
                if (
                    synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                        sub_density.iloc[small_breakpoints[i]
                                        ["main_bp"]]["subgenome3"]
                    ]
                    != 0
                    and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                        sub_density.iloc[small_breakpoints[i]
                                        ["main_bp"]]["subgenome3"]
                    ]
                    == sub_density.iloc[small_breakpoints[i]["main_bp"]]["subgenome3"]
                ):
                    gaps_list_sub3.append(
                        synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][
                            sub_density.iloc[small_breakpoints[i]
                                            ["main_bp"]]["subgenome3"]
                        ]
                    )

                if len(gaps_list_sub1) > 0:
                    sub1_found = False
                    if (
                        sub1_found == False
                        and gaps_list_sub1[0]
                        == sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
                    ):
                        sub1 = sub_density["subgenome1"][small_breakpoints[i]["main_bp"]]
                        # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                        # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                        if sub2_found == False and sub3_found == False:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == False and sub3_found == True:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == True and sub3_found == False:
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        sub1_found = True
                        # if sub2_found and sub3_found:
                        #     break
                        # break

                if len(gaps_list_sub1) == 0:
                    sub1_found = False

                    for b in range(len(synteny_df.columns)):
                        if (
                            sub1_found == False
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            == sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                        ):
                            if sub2_found == False and sub3_found == False:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                sub1_found = True

                            if sub2_found == False and sub3_found == True:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub1_found = True
                                else:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1_found = True
                            if sub2_found == True and sub3_found == False:
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub1_found = True
                                else:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1_found = True
                            if sub2_found == True and sub3_found == True:
                                if (
                                    sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                ):
                                    sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub1_found = True
                                else:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1_found = True
                            # if sub2_found and sub3_found:
                            #     break
                            # break
                    for b in range(len(synteny_df.columns)):
                        if (
                            sub1_found == False
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                        ):

                            if sub2_found == False and sub3_found == False:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                sub1_found = True
                            if sub2_found == False and sub3_found == True:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub1_found = True
                                else:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1_found = True
                            if sub2_found == True and sub3_found == False:
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub1_found = True
                                else:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1_found = True
                            if sub2_found == True and sub3_found == True:
                                if (
                                    sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                ):
                                    sub1 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub1_found = True
                                else:
                                    sub1 = sub_density["subgenome1"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub1_found = True
                            # if sub2_found and sub3_found:
                            #     break
                            # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                            # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                            # break

                if len(gaps_list_sub2) > 0:
                    sub2_found = False
                    if (
                        sub2_found == False
                        and gaps_list_sub2[0]
                        == sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
                    ):
                        # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                        sub2 = sub_density["subgenome2"][small_breakpoints[i]["main_bp"]]
                        # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                        if sub1_found == False and sub3_found == False:
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub1_found == False and sub3_found == True:
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub1_found == True and sub3_found == False:
                            sub3 = sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        sub2_found = True
                        # if sub1_found and sub3_found:
                        #     break
                        # break
                
                elif len(gaps_list_sub2) == 0:
                    sub2_found = False

                    for b in range(len(synteny_df.columns)):
                        if (
                            sub2_found == False
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            == sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                        ):
                            # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]

                            # if sub1_found and sub3_found:
                            #     break
                            if sub1_found == False and sub3_found == False:
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                sub2_found = True

                            if sub1_found == False and sub3_found == True:
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub2_found = True
                                else:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2_found = True

                            if sub1_found == True and sub3_found == False:
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub2_found = True
                                else:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2_found = True
                            if sub1_found == True and sub3_found == True:
                                if (
                                    sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                ):
                                    sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub2_found = True
                                else:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2_found = True

                            # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                            # break
                    for b in range(len(synteny_df.columns)):
                        if (
                            sub2_found == False
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                        ):
                            # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                            # sub3=sub_density['subgenome3'][small_breakpoints[i]['main_bp']]
                            if sub1_found == False and sub3_found == False:
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub2 = synteny_df.iloc[l][b]
                                sub2_found = True
                            if sub1_found == False and sub3_found == True:
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub2_found = True
                                else:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2_found = True
                            if sub1_found == True and sub3_found == False:
                                sub3 = sub_density["subgenome3"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub2_found = True
                                else:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2_found = True
                            if sub1_found == True and sub3_found == True:
                                if (
                                    sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    and sub3 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                ):
                                    sub2 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub2_found = True
                                else:
                                    sub2 = sub_density["subgenome2"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub2_found = True
                            # if sub1_found and sub3_found:
                            #     break
                            # break

                if len(gaps_list_sub3) > 0:
                    sub3_found = False
                    if (
                        sub3_found == False
                        and gaps_list_sub3[0]
                        == sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
                    ):
                        # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                        # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                        sub3 = sub_density["subgenome3"][small_breakpoints[i]["main_bp"]]
                        if sub2_found == False and sub1_found == False:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == False and sub1_found == True:
                            sub2 = sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        if sub2_found == True and sub1_found == False:
                            sub1 = sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ]
                        sub3_found = True
                        # if sub1_found and sub2_found:
                        #     break
                        # break

                if len(gaps_list_sub3) == 0:
                    sub3_found = False
                    for b in range(len(synteny_df.columns)):
                        if (
                            sub3_found == False
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            == sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                        ):
                            # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                            # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]
                            if sub2_found == False and sub1_found == False:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                sub3_found = True
                            if sub2_found == False and sub1_found == True:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub3_found = True
                                else:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3_found = True
                            if sub2_found == True and sub1_found == False:
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub3_found = True
                                else:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3_found = True
                            if sub2_found == True and sub1_found == True:
                                if (
                                    sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    and sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                ):
                                    sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub3_found = True
                                else:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3_found = True
                            # if sub1_found and sub2_found:
                            #     break
                            # break
                    for b in range(len(synteny_df.columns)):
                        if (
                            sub3_found == False
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b] != 0
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome3"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome2"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                            and synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b].split(".")[0]
                            != sub_density["subgenome1"][
                                small_breakpoints[i]["main_bp"]
                            ].split(".")[0]
                        ):
                            # sub1=sub_density['subgenome1'][small_breakpoints[i]['main_bp']]
                            # sub2=sub_density['subgenome2'][small_breakpoints[i]['main_bp']]

                            if sub2_found == False and sub1_found == False:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                sub3_found = True
                            if sub2_found == False and sub1_found == True:
                                sub2 = sub_density["subgenome2"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub3_found = True
                                else:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3_found = True
                            if sub2_found == True and sub1_found == False:
                                sub1 = sub_density["subgenome1"][
                                    small_breakpoints[i]["main_bp"]
                                ]
                                if sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]:
                                    sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub3_found = True
                                else:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3_found = True
                            if sub2_found == True and sub1_found == True:
                                if (
                                    sub1 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    and sub2 != synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                ):
                                    sub3 = synteny_df.iloc[int(small_breakpoints[i]["Row start #"])][b]
                                    sub3_found = True
                                else:
                                    sub3 = sub_density["subgenome3"][
                                        small_breakpoints[i]["main_bp"]
                                    ]
                                    sub3_found = True
                            # if sub1_found and sub2_found:
                            #     break
                            # break
            final_gaps_df = final_gaps_df.append(
                {
                    "Row start #": small_breakpoints[i]["Row start #"],
                    "Row end #": small_breakpoints[i]["Row end #"],
                    "subgenome1": sub1,
                    "subgenome2": sub2,
                    "subgenome3": sub3,
                },
                ignore_index=True,
            )

                # break
        #save the final dataframe to excel
        final_gaps_df = check_gaps(final_gaps_df)
        final_gaps_df.to_excel(final_gaps_file)
              
        return final_gaps_df

def extract_subgenome_columns(final_gaps_files, n_subgenomes):
    sub_columns = []
    for i in range(1, n_subgenomes+1):
        final_gaps_df_sub = pd.read_excel(final_gaps_files[i-1])
        sub_col = final_gaps_df_sub[f"subgenome{i}"]
        sub_columns.append(sub_col)
        if i == 1:
            row_start = final_gaps_df_sub["Row start #"]
            row_end = final_gaps_df_sub["Row end #"]
    sub_columns.insert(0, row_start)
    sub_columns.insert(1, row_end)
    final_gaps_df = pd.concat(sub_columns, axis=1)
    final_gaps_df.to_excel("subgenome_placement_blocks.all.xlsx")

def check_subgenome(df_subgenome, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3):
    for i in range(1, len(df_subgenome)-1):
        if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i, 'subgenome2']:
            if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i-1, 'subgenome1'] or df_subgenome.loc[i, 'subgenome1'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome1'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub2.loc[i, 'subgenome1']
            elif df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i-1, 'subgenome2'] or df_subgenome.loc[i, 'subgenome2'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome2'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub1.loc[i, 'subgenome2']
            elif df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i+1, 'subgenome1'] or df_subgenome.loc[i, 'subgenome1'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome1'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub2.loc[i, 'subgenome1']
            elif df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i+1, 'subgenome2'] or df_subgenome.loc[i, 'subgenome2'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome2'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub1.loc[i, 'subgenome2']
            elif df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i-2, 'subgenome1'] or df_subgenome.loc[i, 'subgenome1'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome1'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub2.loc[i, 'subgenome1']
            elif df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i-2, 'subgenome2'] or df_subgenome.loc[i, 'subgenome2'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome2'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub1.loc[i, 'subgenome2']
            elif df_subgenome.loc[i-1, 'subgenome1'] == df_subgenome_sub2.loc[i, 'subgenome1']:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub2.loc[i, 'subgenome1']
            elif df_subgenome.loc[i-1, 'subgenome2'] == df_subgenome_sub1.loc[i, 'subgenome2']:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub1.loc[i, 'subgenome2']
        if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i, 'subgenome3']:
            if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i-1, 'subgenome1'] or df_subgenome.loc[i, 'subgenome1'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome1'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub3.loc[i, 'subgenome1']
            elif df_subgenome.loc[i, 'subgenome3'] == df_subgenome.loc[i-1, 'subgenome3'] or df_subgenome.loc[i, 'subgenome3'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome3'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub1.loc[i, 'subgenome3']
            elif df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i+1, 'subgenome1'] or df_subgenome.loc[i, 'subgenome1'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome1'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub3.loc[i, 'subgenome1']
            elif df_subgenome.loc[i, 'subgenome3'] == df_subgenome.loc[i+1, 'subgenome3'] or df_subgenome.loc[i, 'subgenome3'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome3'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub1.loc[i, 'subgenome3']
            elif df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i-2, 'subgenome1'] or df_subgenome.loc[i, 'subgenome1'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome1'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub3.loc[i, 'subgenome1']
            elif df_subgenome.loc[i, 'subgenome3'] == df_subgenome.loc[i-2, 'subgenome3'] or df_subgenome.loc[i, 'subgenome3'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome3'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub1.loc[i, 'subgenome3']
            elif df_subgenome.loc[i-1, 'subgenome1'] == df_subgenome_sub3.loc[i, 'subgenome1']:
                df_subgenome.loc[i, 'subgenome1'] = df_subgenome_sub3.loc[i, 'subgenome1']
            elif df_subgenome.loc[i-1, 'subgenome3'] == df_subgenome_sub1.loc[i, 'subgenome3']:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub1.loc[i, 'subgenome3']
        if df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i, 'subgenome3']:
            if df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i-1, 'subgenome2'] or df_subgenome.loc[i, 'subgenome2'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome2'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub3.loc[i, 'subgenome2']
            elif df_subgenome.loc[i, 'subgenome3'] == df_subgenome.loc[i-1, 'subgenome3'] or df_subgenome.loc[i, 'subgenome3'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome3'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub2.loc[i, 'subgenome3']
            elif df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i+1, 'subgenome2'] or df_subgenome.loc[i, 'subgenome2'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome2'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub3.loc[i, 'subgenome2']
            elif df_subgenome.loc[i, 'subgenome3'] == df_subgenome.loc[i+1, 'subgenome3'] or df_subgenome.loc[i, 'subgenome3'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome3'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub2.loc[i, 'subgenome3']
            elif df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i-2, 'subgenome2'] or df_subgenome.loc[i, 'subgenome2'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome2'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub3.loc[i, 'subgenome2']
            elif df_subgenome.loc[i, 'subgenome3'] == df_subgenome.loc[i-2, 'subgenome3'] or df_subgenome.loc[i, 'subgenome3'].split(".")[0] == df_subgenome.loc[i-1, 'subgenome3'].split(".")[0]:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub2.loc[i, 'subgenome3']
            elif df_subgenome.loc[i-1, 'subgenome2'] == df_subgenome_sub3.loc[i, 'subgenome2']:
                df_subgenome.loc[i, 'subgenome2'] = df_subgenome_sub3.loc[i, 'subgenome2']
            elif df_subgenome.loc[i-1, 'subgenome3'] == df_subgenome_sub2.loc[i, 'subgenome3']:
                df_subgenome.loc[i, 'subgenome3'] = df_subgenome_sub2.loc[i, 'subgenome3']

    df_subgenome.to_excel("subgenome_placement_blocks.all.xlsx")
    return df_subgenome

#check whether the subgenome columns in finalgaps_df have duplicates in between them in same row
def error_check(output_final_df, df_subgenome, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3):
    for i in range(1, len(output_final_df)):
            if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i, 'subgenome2']:
                check_subgenome(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)
            if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i, 'subgenome3']:
                check_subgenome(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)
            if df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i, 'subgenome3']:
                check_subgenome(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)


def check_neighbourhoods(df_subgenome, window_size_sub1, window_size_sub2, window_size_sub3):
    subgenome1 = df_subgenome["subgenome1"]
    subgenome2 = df_subgenome["subgenome2"]
    subgenome3 = df_subgenome["subgenome3"]
    # window_size_sub1 = 4
    # window_size_sub2 = 15
    # window_size_sub3 = 15

    for i in range(len(df_subgenome)):
        if i > window_size_sub1:
            substring = subgenome1[i - window_size_sub1 : i]
            match_percentage_sub1 = sum(substring == subgenome1[i]) / len(substring) * 100
            match_percentage_sub2 = sum(substring == subgenome2[i]) / len(substring) * 100
            match_percentage_sub3 = sum(substring == subgenome3[i]) / len(substring) * 100
            # Calculate the percentage of matching between substring and subgenome1[i]

            if match_percentage_sub2 > 90 and match_percentage_sub3 < 70 and match_percentage_sub1 < 90:
                # Perform your desired operations based on the match_percentage value
                # For example, you can swap subgenome1 and subgenome2
                df_subgenome.at[i, "subgenome1"], df_subgenome.at[i, "subgenome2"] = df_subgenome.at[i, "subgenome2"], df_subgenome.at[i, "subgenome1"]
            if match_percentage_sub3 > 90 and match_percentage_sub2 < 70 and match_percentage_sub1 < 90:
                # Perform your desired operations based on the match_percentage value
                # For example, you can swap subgenome1 and subgenome2
                df_subgenome.at[i, "subgenome1"], df_subgenome.at[i, "subgenome3"] = df_subgenome.at[i, "subgenome3"], df_subgenome.at[i, "subgenome1"]
        if i > window_size_sub2:
            substring = subgenome2[i - window_size_sub2 : i]
            match_percentage_sub1 = sum(substring == subgenome1[i]) / len(substring) * 100
            match_percentage_sub2 = sum(substring == subgenome2[i]) / len(substring) * 100
            match_percentage_sub3 = sum(substring == subgenome3[i]) / len(substring) * 100
            # Calculate the percentage of matching between substring and subgenome1[i]

            if match_percentage_sub1 > 90 and match_percentage_sub3 < 70 and match_percentage_sub2 < 70:
                # Perform your desired operations based on the match_percentage value
                # For example, you can swap subgenome1 and subgenome2
                df_subgenome.at[i, "subgenome1"], df_subgenome.at[i, "subgenome2"] = df_subgenome.at[i, "subgenome2"], df_subgenome.at[i, "subgenome1"]
            if match_percentage_sub3 > 90 and match_percentage_sub2 < 70 and match_percentage_sub1 < 70:
                # Perform your desired operations based on the match_percentage value
                # For example, you can swap subgenome1 and subgenome2
                df_subgenome.at[i, "subgenome2"], df_subgenome.at[i, "subgenome3"] = df_subgenome.at[i, "subgenome3"], df_subgenome.at[i, "subgenome2"]
        if i > window_size_sub3:
            substring = subgenome3[i - window_size_sub3 : i]
            match_percentage_sub1 = sum(substring == subgenome1[i]) / len(substring) * 100
            match_percentage_sub2 = sum(substring == subgenome2[i]) / len(substring) * 100
            match_percentage_sub3 = sum(substring == subgenome3[i]) / len(substring) * 100
            # Calculate the percentage of matching between substring and subgenome1[i]

            if match_percentage_sub1 > 90 and match_percentage_sub3 < 90 and match_percentage_sub2 < 70:
                # Perform your desired operations based on the match_percentage value
                # For example, you can swap subgenome1 and subgenome2
                df_subgenome.at[i, "subgenome1"], df_subgenome.at[i, "subgenome3"] = df_subgenome.at[i, "subgenome3"], df_subgenome.at[i, "subgenome1"]
            if match_percentage_sub2 > 90 and match_percentage_sub2 < 70 and match_percentage_sub3 < 70:
                # Perform your desired operations based on the match_percentage value
                # For example, you can swap subgenome1 and subgenome2
                df_subgenome.at[i, "subgenome2"], df_subgenome.at[i, "subgenome3"] = df_subgenome.at[i, "subgenome3"], df_subgenome.at[i, "subgenome2"]
            
    return df_subgenome

