
import numpy as np
import pandas as pd
import re
import warnings
import sys
import pickle
import csv
import os
import SyntenyLink_mbp as mbp
import SyntenyLink_wg as wg
import SyntenyLink_sb as sb
import SyntenyLink_mn as mn
import SyntenyLink_acc as acc
import math

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
        acc.subgenome_overlap(GT,"Super_synteny_bl_sub_placement_density.xlsx", df_synteny, 3, first_letter_get)
        #Input the file for chains generated in DagChainer as a command line argument
        chains_file = sys.argv[sys.argv.index('-c') + 1]
        #Input the file for blastn generated in DagChainer as a command line argument
        blastn_file = sys.argv[sys.argv.index('-bl') + 1] 

        num_blocks_main = len(pd.read_excel("Super_synteny_bl_sub_placement_density.xlsx"))
        wg.node_traverse_most_weighted_path(GT, n_subgenomes,df_synteny, chains_file, blastn_file, C_df_new, num_blocks_main, first_letter_get)
        print("===========================================")
        print(f"Number of blocks: {num_blocks_main}")
        print("===========================================")
        num_blocks = len(pd.read_excel("Super_synteny_bl_sub_placement_density.xlsx"))

        for k in range(n_subgenomes):
            #Focus on rows having more than 3 non_zero density columns
            main_bp = pd.read_excel("Super_synteny_block_output_non_zero.xlsx")
            temp_file = pd.read_excel("Super_synteny_bl_sub_placement_density.xlsx")
            main_bp_density_rmv = pd.read_excel(f"Super_synteny_graph_nodes_sub{k+1}.xlsx")
            #Append non_zero column from temp file to main_bp_density_rmv
            main_bp_density_rmv.insert(2, "Non_zero", temp_file["Non_zero"])
            main_bp_density_rmv = main_bp_density_rmv.loc[:, ~main_bp_density_rmv.columns.str.contains('^Unnamed')]
            #remove any column starts as Unnamed
            main_bp = main_bp.loc[:, ~main_bp.columns.str.contains('^Unnamed')]
            #append the first 6 columns of main_bp_density_rmv to main_bp
            main_bp.insert(0, "Row start #", main_bp_density_rmv["Row start #"])
            main_bp.insert(1, "Row end #", main_bp_density_rmv["Row end #"])

            # Create a new dataframe
            df = pd.DataFrame(columns=main_bp.columns)

            for i in range(len(main_bp)):
                for j in range(len(main_bp.columns)):
                    if (main_bp.iloc[i]['Non_zero'] > 3) or (main_bp.iloc[i][j] != 0 and main_bp.iloc[i][j] < 0.1):
                            df = df.append(main_bp.iloc[i])

            # Drop duplicate rows
            df = df.drop_duplicates()

            # Reset the index of the dataframe
            df.reset_index(drop=True, inplace=True)
            sb.update_df_synteny(C_df_new, value, "Super_synteny_bl_sub_placement_density.xlsx", f"modified_chr_names{k+1}.xlsx", main_bp_density_rmv)

            #extract only the columns starting with N and keep Row start # and Row end # two columns as well
            df_focus = df.filter(regex=r'^N\d+|Row start #|Row end #|Non_zero')
            chr_names_C = pd.read_excel(f"modified_chr_names{k+1}.xlsx")
            # Extract corresponding rows from C numpy array geting Row start # and Row end # values and combine them to a new numpy array
            C_new = pd.DataFrame()
            for i in range(len(df_focus)):
                row_start = int(df_focus.iloc[i]['Row start #'])
                row_end = int(df_focus.iloc[i]['Row end #'])
                #extract the rows from C
                C_new = C_new.append(chr_names_C[row_start:row_end+1])
            
            # C_new.to_excel(f"low_density_synteny{k+1}.xlsx")

            new_df_focus = sb.small_blk(df_focus, main_bp_density_rmv, chr_names_C)
            #if row start # is greater than row end #, remove the row
            new_df_focus = new_df_focus[new_df_focus['Row start #'] <= new_df_focus['Row end #']]
            new_df_focus = new_df_focus.reset_index(drop=True)

            # replace the nan values in new_df_focus with the values in previous cell
            for i in range(len(new_df_focus)):
                if pd.isnull(new_df_focus.iloc[i]['subgenome1']):
                    new_df_focus.at[i, 'subgenome1'] = new_df_focus.iloc[i - 1]['subgenome1']
                if pd.isnull(new_df_focus.iloc[i]['subgenome2']):
                    new_df_focus.at[i, 'subgenome2'] = new_df_focus.iloc[i - 1]['subgenome2']
                if pd.isnull(new_df_focus.iloc[i]['subgenome3']):
                    new_df_focus.at[i, 'subgenome3'] = new_df_focus.iloc[i - 1]['subgenome3']

            #sort new_df_focus by row start #
            new_df_focus = new_df_focus.sort_values(by=['Row end #'])

            new_df_focus = sb.multi_col(new_df_focus)

            #remove duplicated row end # values rows keeping first occurence
            new_df_focus = new_df_focus.drop_duplicates(subset=['Row end #'], keep='first')

            new_sbp = sb.all_blk(new_df_focus, main_bp_density_rmv)
            new_sbp.to_excel(f'subgenome_placement_blocks.all{k+1}.xlsx', index=False)
        
        print("=========================================")
        print("Number of small + main blocks:" , len(pd.read_excel("subgenome_placement_blocks.all1.xlsx")))
        

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

        final_gaps_files = ["subgenome_placement_blocks.all1.xlsx", "subgenome_placement_blocks.all2.xlsx", "subgenome_placement_blocks.all3.xlsx"]
        extract_subgenome_columns(final_gaps_files, n_subgenomes)
        # Read the data
        df_subgenome_density = pd.read_excel("subgenome_placement_blocks.all.xlsx")

        # Create a list of subgenomes
        subgenomes = [[] for _ in range(n_subgenomes)]

        for i in range(len(df_subgenome_density["Row start #"])):
            for j in range(df_subgenome_density["Row start #"].values[i], df_subgenome_density["Row end #"].values[i] + 1):
                for k in range(n_subgenomes):
                    subgenomes[k].append(df_synteny[df_subgenome_density[f"subgenome{k+1}"].values[i]][j])

        # Create a dataframe of subgenomes
        df_subgenome = pd.DataFrame(subgenomes).transpose()

        # Change the column names
        df_subgenome.columns = [f"subgenome{i+1}" for i in range(n_subgenomes)]
        #saving the subgenome assignment to excel
        # df_subgenome.to_excel("subgenome_assignment.xlsx")
        df_subgenome_sub1 = pd.read_excel("subgenome_placement_blocks.all1.xlsx")
        df_subgenome_sub2 = pd.read_excel("subgenome_placement_blocks.all2.xlsx")
        df_subgenome_sub3 = pd.read_excel("subgenome_placement_blocks.all3.xlsx")

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
                


        output_final_df = check_subgenome(df_subgenome_density, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)

        #check whether the subgenome columns in finalgaps_df have duplicates in between them in same row
        def error_check(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3):
            for i in range(1, len(output_final_df)):
                    if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i, 'subgenome2']:
                        check_subgenome(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)
                    if df_subgenome.loc[i, 'subgenome1'] == df_subgenome.loc[i, 'subgenome3']:
                        check_subgenome(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)
                    if df_subgenome.loc[i, 'subgenome2'] == df_subgenome.loc[i, 'subgenome3']:
                        check_subgenome(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)

        error_check(output_final_df, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)

        print("===========================================")
        print(f"Accuracy of subgenome assignment(Main BP + Small BP):")
        print("===========================================")
        acc.subgenome_overlap(GT,"subgenome_placement_blocks.all.xlsx", df_synteny, 3, first_letter_get)
        # Convert the DataFrame to a NumPy array
        # Read the specified range of cells into a pandas DataFrame
        df = pd.read_excel("subgenome_placement_blocks.all.xlsx", header= None, usecols="E:G", skiprows=1)
        subgenomes = np.array(df)
        len_subgenomes = len(subgenomes[:,0])
        print(len_subgenomes)

        window_up= 3
        window_down= 19
        arg_max_count_2=np.array([-1,-1,-1])
        for k in range(len_subgenomes):
            if k<window_up:
                start_ind_2=0
                end_ind_2=k+window_down
            elif k< len_subgenomes-window_down:
                start_ind_2=k-window_up
                end_ind_2=k+window_down
            else:
                start_ind_2=k-window_up
                end_ind=len_subgenomes-1
            for i in range(3):
                #find the column that has the maximum count for all_head[i]
                count1_2=mn.count_occurrences(subgenomes[start_ind_2:end_ind_2,0],subgenomes[k,i])
                count2_2=mn.count_occurrences(subgenomes[start_ind_2:end_ind_2,1],subgenomes[k,i])
                count3_2=mn.count_occurrences(subgenomes[start_ind_2:end_ind_2,2],subgenomes[k,i])
                counts_2=np.array([count1_2,count2_2,count3_2])
                arg_max_count_2[i]=np.argmax(counts_2)
            
            arg_max_count_2=np.array(arg_max_count_2) 
        
            if not(np.array_equal(arg_max_count_2,np.array([0,1,2]))):
            #if the entries do not belong to their higher dense columns
                for i in range(3):
                    elements_2=np.where(arg_max_count_2==i)
                    if len(elements_2[0])==1:
                            temp=subgenomes[k,elements_2[0][0]]
                            subgenomes[k,i]=subgenomes[k,elements_2[0][0]]
                            subgenomes[k,elements_2[0][0]]=temp
                    elif len(elements_2[0])>1:
                        #elements that are competing for column i, that is, there is a tie
                        for j in range(len(elements_2[0])):
                            #loop over the elements that are competing for column i
                            if k == 0:
                                previous_2=1
                            else: 
                                previous_2=k-1
                            if mn.compare_subgenomes(subgenomes[k-1,i],subgenomes[k,elements_2[0][j]]) and subgenomes[k,i]!=subgenomes[k,elements_2[0][j]]:
                                #if the next entry in the column i is equal to entry j 
                                # and that entry is already not in column i 
                                # then swap them
                                temp=subgenomes[k,i]
                                subgenomes[k,i]=subgenomes[k,elements_2[0][j]]
                                subgenomes[k,elements_2[0][j]]=temp

        cols = pd.read_excel("subgenome_placement_blocks.all.xlsx", header= None, usecols="C:D", skiprows=1)
        df_new_2 = pd.DataFrame()
        df_new_2['Row start #'] = cols[2]
        df_new_2['Row end #'] = cols[3]
        df_new_2['subgenome1']=subgenomes[:,0]
        df_new_2['subgenome2']=subgenomes[:,1]
        df_new_2['subgenome3']=subgenomes[:,2]
        df_new_2.to_excel("Final_subgenome_placement_result.xlsx")

        # Read the data
        df_subgenome_density = pd.read_excel("Final_subgenome_placement_result.xlsx")

        # Create a list of subgenomes
        subgenomes = [[] for _ in range(n_subgenomes)]

        for i in range(len(df_subgenome_density["Row start #"])):
            for j in range(df_subgenome_density["Row start #"].values[i], df_subgenome_density["Row end #"].values[i] + 1):
                for k in range(n_subgenomes):
                    subgenomes[k].append(df_synteny[df_subgenome_density[f"subgenome{k+1}"].values[i]][j])

        # Create a dataframe of subgenomes
        df_subgenome = pd.DataFrame(subgenomes).transpose()

        # Change the column names
        df_subgenome.columns = [f"subgenome{i+1}" for i in range(n_subgenomes)]
        df_subgenome.to_excel("Final_result.xlsx")


        acc.subgenome_overlap(GT,"Final_subgenome_placement_result.xlsx", df_synteny, 3, first_letter_get)

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