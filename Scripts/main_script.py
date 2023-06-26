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

def main():
    #get the input file as an argument (collinear file)

    input_file = sys.argv[sys.argv.index('-i') + 1]

    #take gap_threshold and min_block_length which are numbers as arguments
    gap_threshold = int(sys.argv[sys.argv.index('-g') + 1])
    min_block_length = int(sys.argv[sys.argv.index('-m') + 1])

    #Get the number of subgenomes as an argument
    n_subgenomes = int(sys.argv[sys.argv.index('-n') + 1])

    #convert the collinear file to a dataframe
    C_df_csv = pd.read_csv(input_file, sep = '\t', header=None)
    #make a copy of the dataframe
    C_df = C_df_csv.copy()
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


    #convert the dataframe to a numpy array
    C = C_df_updated.to_numpy()
    m, n = C.shape
    print(m, n)


    break_point_indices = mbp.gap_calculation(C, gap_threshold, min_block_length, n, m)
    densities = mbp.get_densities(C, break_point_indices, n, m)
    mbp.create_excel_sheet(C_df, break_point_indices, densities, C_df_updated_copy)
    df_temp = mbp.get_subgenomes("Super_synteny_block_output.xlsx", n_subgenomes)
    df = mbp.assign_subgenomes(df_temp, "Super_synteny_block_output.xlsx", n_subgenomes)
    #Saving the file
    df.to_excel("Super_synteny_bl_sub_placement_density.xlsx")
    # Get the column names from df_subgenome_density that start with N followed by a number
    column_names = df.filter(regex=r'^N\d+').columns.tolist()

    # Replace the columns starting from the second column until the end of C_df_updated with the selected column names
    df_synteny = C_df_csv.iloc[1:,1:-3].rename(columns={C_df_csv.iloc[1:,1:-3].columns[i-1]: column_names[i-1] for i in range(1, len(column_names)+1)})
    #append the first column of C_df to the first column of df_synteny
    df_synteny.insert(0, "locus_id", C_df_csv.iloc[1:,0])
    #update first row index starting from 0
    df_synteny.index = df_synteny.index - 1

    #Input the file for chains generated in DagChainer as a command line argument
    chains_file = sys.argv[sys.argv.index('-c') + 1]
    #Input the file for blastn generated in DagChainer as a command line argument
    blastn_file = sys.argv[sys.argv.index('-bl') + 1] 

    num_blocks_main = len(pd.read_excel("Super_synteny_bl_sub_placement_density.xlsx"))
    wg.node_traverse_most_weighted_path(n_subgenomes,df_synteny, chains_file, blastn_file, C_df_csv, num_blocks_main)
    print("===========================================")
    print(f"Number of blocks: {num_blocks_main}")
    print("===========================================")

    num_blocks = len(pd.read_excel("Super_synteny_bl_sub_placement_density.xlsx"))

    # Determine file names based on number of subgenomes
    subgenome_density_files = "Super_synteny_bl_sub_placement_density.xlsx" 
    gene_matrix_output_files = "Super_synteny_block_output.xlsx" 

    # Define the small_BP_subi list to hold small breakpoints for each subgenome
    small_BP_sub = []
    for i in range(n_subgenomes):
        # Generate small breakpoints
        small_BP_subi = sb.remove_invalid_entries(
            sb.modify_small_breakpoints(
                sb.merge_main_small_breakpoints(
                    sb.update_small_breakpoints(
                        sb.find_small_breakpoints(
                            sb.sort_gaps(
                                sb.merge_gaps(
                                    sb.merge_gaps(
                                        sb.get_gaps(
                                            gene_matrix_output_files, 
                                            subgenome_density_files, 
                                            C, 
                                            n_subgenomes, 
                                            1
                                        )
                                    )
                                )
                            ), 
                            subgenome_density_files
                        )
                    ), 
            subgenome_density_files
                )
            )
        )
        small_BP_sub.append(small_BP_subi)
        sb.update_df_synteny(C_df, C_df_updated, subgenome_density_files, f"abc_synteny_chromosome_names.success.colinear{i+1}.xlsx")

    for k in range(n_subgenomes):
        sb.subgenome_assignment_all_BP(f"Super_synteny_graph_nodes_sub{k+1}.xlsx", f"abc_synteny_chromosome_names.success.colinear{k+1}.xlsx", small_BP_sub[k], f"subgenome_placement_blocks.all.sub{k+1}.xlsx")
    print("=========================================")
    print("Number of small + main blocks:" , len(pd.read_excel("subgenome_placement_blocks.all.sub1.xlsx")))
    print("=========================================")

    final_gaps_files = ["subgenome_placement_blocks.all.sub1.xlsx", "subgenome_placement_blocks.all.sub2.xlsx", "subgenome_placement_blocks.all.sub3.xlsx"]
    sb.extract_subgenome_columns(final_gaps_files, n_subgenomes)
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
    df_subgenome_sub1 = pd.read_excel("subgenome_placement_blocks.all.sub1.xlsx")
    df_subgenome_sub2 = pd.read_excel("subgenome_placement_blocks.all.sub2.xlsx")
    df_subgenome_sub3 = pd.read_excel("subgenome_placement_blocks.all.sub3.xlsx")

    output_final_df = sb.check_subgenome(df_subgenome_density, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)
    sb.error_check(output_final_df, df_subgenome, df_subgenome_sub1, df_subgenome_sub2, df_subgenome_sub3)

    window_size_sub1 = int(sys.argv[sys.argv.index('-ws1') + 1])
    window_size_sub2 = int(sys.argv[sys.argv.index('-ws2') + 1])
    window_size_sub3 = int(sys.argv[sys.argv.index('-ws3') + 1])

    df_subgenome = sb.check_neighbourhoods(output_final_df, window_size_sub1, window_size_sub2, window_size_sub3)
    df_subgenome.to_excel("subgenome_placement_blocks.all.xlsx")

    # # Read the specified range of cells into a pandas DataFrame
    df = pd.read_excel("subgenome_placement_blocks.all.xlsx", header= None, usecols="E:G", skiprows=1)

    # Convert the DataFrame to a NumPy array
    subgenomes = np.array(df)
    len_subgenomes = len(subgenomes[:,0])
    print(len_subgenomes)

    window_up= int(sys.argv[sys.argv.index('-wup') + 1])
    window_down= int(sys.argv[sys.argv.index('-wdwn') + 1])
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
    df_subgenome.to_excel("Final_result.xlsx")

    # #If there exist a ground truth file, then compare the results with the ground truth
    # GT = sys.argv[sys.argv.index('-gt') + 1]
    # acc.subgenome_overlap(GT,"subgenome_placement_blocks.all.xlsx", df_synteny, 3)

if __name__ == '__main__':
    main()