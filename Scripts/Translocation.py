import graph_ADT as G
import pandas as pd
import numpy as np

class Translocation:
    def __init__(self):
        # Create a breakpoint class object
        self.synteny = G.Graph()

    def main(self):
        # # Run the main function
        adjusted_subgenomes_stranded_segmental_original, dag, GT, n_sub, header, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns, original_ploidy_numpy, last_column, ploidy_numpy_overlap, subgenome_letters = self.synteny.main()
        # print("Adjusted subgenomes stranded segmental original:", adjusted_subgenomes_stranded_segmental_original)

        # Parse synteny blocks
        block_dict = self.synteny.parse_synteny_blocks(dag)
        # print("Block dict:", block_dict)

        synteny_block = self.synteny.process_synteny_blocks(block_dict)

        gene_combinations = self.synteny.weight_graph(adjusted_subgenomes_stranded_segmental_original)
        # print("Gene combinations:", gene_combinations)
        #save gene_combinations as an excel file
        # gene_combinations_df = pd.DataFrame(gene_combinations)
        # gene_combinations_df.to_excel('gene_combinations.xlsx')

        
        graphs = self.synteny.create_graphs(gene_combinations, synteny_block)
        formatted_dfs = {}
        for key, G in graphs.items():
            # print(f"Graph for {key}:")
            # print(G.nodes(data=True))
            # print(G.edges(data=True))

            traversed = self.synteny.traverse_highest_weighted_path(G, key, n_sub)
            #convert the traversed paths to a dataframe which is a list
            traversed_dict = {}
            for k, v in traversed.items():
                traversed_dict[k] = pd.DataFrame(v)

            
            formatted_dfs[key] = traversed_dict
            
        return formatted_dfs, header, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns, GT, n_sub, original_ploidy_numpy, last_column, ploidy_numpy_overlap, subgenome_letters
    
    def get_subgenome_chr(self, adjusted_subgenome, header, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns):
        """
        Purpose: Assigns each gene to a chromosome based on the subgenome assignments.
        Args:
            adjusted_subgenome: Subgenome assignments, header: Headers from the collinear file., ploidy_status: The ploidy status of the organism., chrsets: The chromosome sets., file: The collinear file.
        Returns:
            A dictionary with subgenome-chromosome assignments.
        """
        # Calculate the total number of columns needed for each chrset based on its value.
        total_columns_needed = {key: int(value) * 2 for key, value in chrsets.items()}
        chr_names = {}

        # Create a cumulative total to keep track of the column ranges for each chrset.
        cumulative_total_columns = {}
        cumulative_sum = 0
        for key, value in total_columns_needed.items():
            cumulative_total_columns[key] = (cumulative_sum, cumulative_sum + value)
            cumulative_sum += value

        # Iterate over each key in adjusted_subgenome to map each subgenome to its chromosomes.
        for key, subgenome_dict in adjusted_subgenome.items():
            gene_id_matrix_df = gene_id_matrix_df_[key].reset_index(drop=True)
            # Concatenate the last two columns to the gene_id_matrix_df
            gene_id_matrix_df_extended = pd.concat([gene_id_matrix_df, last_two_columns], axis=1)

            # Determine the chrset key corresponding to the current adjusted_subgenome key
            chrset_key = None
            for k in chrsets:
                if k in key:  # or any other logic that matches the keys
                    chrset_key = k
                    break

            if not chrset_key:
                continue  # Skip if no corresponding chrset is found

            # Apply the relevant column names based on chrset.
            start_col, end_col = cumulative_total_columns[chrset_key]
            relevant_columns = header[start_col:end_col]
            # Add last two columns to the relevant_columns
            relevant_columns = relevant_columns + header[-3:-1]
            gene_id_matrix_df_extended.columns = relevant_columns

            # Create a DataFrame to hold the subgenome-chromosome assignments.
            subgenome_dict_new = {subgenome_num: [] for subgenome_num in subgenome_dict.keys()}
            # Replace gene IDs in df with chromosome names.
            for subgenome_num, df in subgenome_dict.items():
                for i, row in df.iterrows():
                    for j in range(1, len(row) + 1):
                        gene_id = row[j-1]  # Correctly get the gene_id from the row
                        if gene_id != 'x':
                            found = False
                            for col in gene_id_matrix_df_extended.columns:
                                # Find the column where the gene ID is located
                                if gene_id in gene_id_matrix_df_extended[col].values:
                                    # Replace the gene ID in the adjusted_subgenome with the column name (chromosome)
                                    subgenome_dict_new[subgenome_num].append(col)
                                    found = True
                                    break  # Move to the next gene ID once replaced

                            if not found:
                                subgenome_dict_new[subgenome_num].append('overlap')
                            
                        else:
                            subgenome_dict_new[subgenome_num].append('x')

            # Convert the subgenome-chromosome dictionary to a DataFrame
            subgenome_dict_df = pd.DataFrame.from_dict(subgenome_dict_new, orient='index').transpose()

            # Save to Excel and add to the dictionary
            excel_filename = f'{key}_all_subgenome_Chr.xlsx'
            subgenome_dict_df.to_excel(excel_filename, index=False)
        
            chr_names[key] = subgenome_dict_new

        # After replacing, 'adjusted_subgenome' will have chromosome names instead of gene IDs.
        return chr_names
    
    def update_chr_x_with_previous_or_next_chr(self, df):
        """
        Purpose: This function updates 'x' values in subgenome columns with the previous or next chromosome name.
        Args: df
        Returns: df
        """
        for key, df_sub in df.items():
            # Convert the subgenome_df to a DataFrame
            df_sub = pd.DataFrame(df_sub)
            for i, row in df_sub.iterrows():
                for j in range(len(row)):
                    gene_id = row[j]
                    if gene_id == 'x':
                        # Find the previous non-'x' value in the same column
                        previous_gene_id = None
                        for k in range(i - 1, -1, -1):
                            previous_gene_id = df_sub.iloc[k, j]
                            if previous_gene_id != 'x':
                                df_sub.iloc[i, j] = previous_gene_id
                                break
                        # If no previous non-'x' value is found, find the next non-'x' value in the same column
                        if previous_gene_id == 'x' or previous_gene_id is None:
                            for k in range(i + 1, len(df_sub)):
                                next_gene_id = df_sub.iloc[k, j]
                                if next_gene_id != 'x':
                                    df_sub.iloc[i, j] = next_gene_id
                                    break
            df[key] = df_sub
            # Save the updated DataFrame to an Excel file
            excel_filename = f'{key}_all_subgenome_Chr_x_updated.xlsx'
            df_sub.to_excel(excel_filename, index=False)

        return df
    
    def remove_inconsistencies_in_subgenome_assignments(self, df):
        """
        Purpose: This function removes inconsistencies in subgenome assignments.
        Args: df
        Returns: df
        """
        for key, df_sub in df.items():
            # Convert the subgenome_df to a DataFrame
            df_sub = pd.DataFrame(df_sub)
            num_rows, num_cols = df_sub.shape

            for col in range(num_cols):
                chromosome_count = 1
                last_chromosome = df_sub.iloc[0, col]
                for row in range(1, num_rows):
                    current_chromosome = df_sub.iloc[row, col]
                    if current_chromosome == last_chromosome:
                        chromosome_count += 1
                    else:
                        if chromosome_count < 3:
                            consistent_chromosome = None
                            for other_col in range(num_cols):
                                if other_col == col:
                                    continue
                                other_chromosome = df_sub.iloc[row, other_col]
                                if other_chromosome == last_chromosome:
                                    consistent_chromosome = other_chromosome
                                    break
                            if consistent_chromosome:
                                df_sub.iloc[row-chromosome_count:row, col] = consistent_chromosome
                            else:
                                df_sub.iloc[row-chromosome_count:row, col] = 'x'
                        chromosome_count = 1
                    last_chromosome = current_chromosome

                if chromosome_count < 3:
                    df_sub.iloc[num_rows-chromosome_count:num_rows, col] = 'x'

            df[key] = df_sub
            # Save the updated DataFrame to an Excel file
            excel_filename = f'{key}_all_subgenome_Chr_inconsistencies_removed.xlsx'
            df_sub.to_excel(excel_filename, index=False)

        return df




    # def update_chr_names_overlap(self, df):
    #     num_key = 0    
    #     condition_met = False
    #     while not condition_met:
    #         condition_met = True
    #         for key, df_sub in df.items():
    #             for i, row in df_sub.iterrows():
    #                 for j in range(1, len(row) + 1):
    #                     subgenome = f'subgenome_{j}'
    #                     gene_id = row[subgenome]
    #                     if gene_id != 'x':
    #                         # If length of gene id is less than 10
    #                         if len(gene_id) >= 10:
    #                             condition_met = False
    #                             # Get the three characters before character 'g' in the gene id.
    #                             chr_id = gene_id[gene_id.find('g')-3:gene_id.find('g')]
    #                             # If there is a 0 between letter and number, remove it
    #                             chr_id = chr_id[0] + str(int(chr_id[2])+(num_key*10))
    #                             chr_id_reverse = chr_id + '.r'
                                
    #                             # Check if chr_id or chr_id_reverse is already in another subgenome in the same row
    #                             if chr_id in row.values or chr_id_reverse in row.values:
    #                                 df_sub.at[i, subgenome] = 'x'
    #                             else:
    #                                 # Check if chr_id is in any column of the previous or next 10 rows
    #                                 if (i - 10 >= 0 and chr_id in df_sub.iloc[i-10:i].values) or (i + 10 < len(df_sub) and chr_id in df_sub.iloc[i+1:i+11].values):
    #                                     df_sub.at[i, subgenome] = chr_id
    #                                 # Check if chr_id_reverse is in any column of the previous or next 10 rows
    #                                 elif (i - 10 >= 0 and chr_id_reverse in df_sub.iloc[i-10:i].values) or (i + 10 < len(df_sub) and chr_id_reverse in df_sub.iloc[i+1:i+11].values):
    #                                     df_sub.at[i, subgenome] = chr_id_reverse
    #                                 else:
    #                                     df_sub.at[i, subgenome] = chr_id
                

    #             num_key += 1
    #             # Save the updated dataframe to an Excel file
    #             excel_filename = f'{key}_all_subgenome_Chr_overlap.xlsx'
    #             df_sub.to_excel(excel_filename, index=False)
        
    #             df[key] = df_sub
            
    #     return df


    def chr_block_subgenome_assignment(self, adjusted_subgenome_dict):
        """
        Purpose: This function assigns subgenomes to each block based on the adjusted subgenome assignments.
        Args: adjusted_subgenome_dict
        Returns: final_dict
        """
        final_dict = {}

        # Iterate through each key's subgenome data
        for key, subgenome_df in adjusted_subgenome_dict.items():
            all_blocks_df = pd.DataFrame()  # Initialize here to reset for each key
            
            #convert the subgenome_df to a dataframe
            subgenome_df = pd.DataFrame(subgenome_df)

            # Iterate through each column (subgenome)
            for subgenome in subgenome_df.columns:
                blocks = []
                current_chr = None
                block_start = None
                
                # Iterate through rows to create blocks
                for index, row in subgenome_df.iterrows():
                    chr_name = row[subgenome]
                    if chr_name != 'x':  # Ignore 'x', focus on chromosome names
                        if chr_name != current_chr:
                            # End the previous block if there was one
                            if current_chr is not None:
                                blocks.append({'start': block_start, 'end': index - 1, 'chr': current_chr})
                            # Start a new block
                            current_chr = chr_name
                            block_start = index
                    # Capture the last block if we're at the end and there is an ongoing block
                    if index == len(subgenome_df) - 1 and current_chr is not None:
                        blocks.append({'start': block_start, 'end': index, 'chr': current_chr})

                # Convert blocks to DataFrame for this subgenome
                if blocks:
                    blocks_df = pd.DataFrame(blocks)
                    blocks_df['subgenome'] = subgenome
                    # Append to the all_blocks_df for this key
                    all_blocks_df = pd.concat([all_blocks_df, blocks_df], ignore_index=True)

            # Now we have all blocks for this key, pivot this DataFrame to have subgenomes as columns
            if not all_blocks_df.empty:
                final_df = all_blocks_df.pivot_table(index=['start', 'end'], columns='subgenome', values='chr', aggfunc='first').reset_index()
                final_df.columns.name = None  # Remove the pivot table multi-index naming
                # final_df.fillna('x', inplace=True)  # Replace NaN with 'x' for uniformity
                #replacing NaN with previous value
                final_df.fillna(method='ffill', inplace=True)
                final_df.sort_values(by=['start', 'end'], inplace=True)
                final_df.reset_index(drop=True, inplace=True)

                # Process to check continuity between ends and starts
                end_adjustments = []
                for i in range(1, len(final_df)):
                    if final_df.at[i, 'start'] == final_df.at[i - 1, 'end']:
                        for col1 in [col for col in final_df.columns if 'subgenome' in col]:
                            for col2 in [col for col in final_df.columns if 'subgenome' in col]:
                                if final_df.at[i, col1] == final_df.at[i - 1, col2] and final_df.at[i, col1] != 'x' and subgenome_df.at[final_df.at[i - 1, 'end'], col2] == 'x':
                                    end_adjustments.append((i - 1, final_df.at[i - 1, 'end'] - 1))

                # Apply end adjustments after all checks are complete to avoid conflicts
                for i, new_end in end_adjustments:
                    final_df.at[i, 'end'] = new_end
                
            else:
                # If no blocks, create an empty DataFrame with the correct columns
                final_df = pd.DataFrame(columns=['start', 'end'] + [f'subgenome_{i}' for i in range(1, len(subgenome_df.columns) + 5)])

            # Save to Excel and add to the dictionary
            excel_filename = f'{key}_chr_subgenomes.xlsx'
            final_df.to_excel(excel_filename, index=False)
            final_dict[key] = final_df

        return final_dict
    
    def remove_overlapping_assignments(self, adjusted_subgenome_dict):
        """
        Purpose: This function removes overlapping subgenome assignments.
        Args: adjusted_subgenome_dict
        Returns: cleaned_subgenome_dict
        """
        cleaned_subgenome_dict = {}
        
        for key, df in adjusted_subgenome_dict.items():
            # Copy the DataFrame to avoid changing the original data
            new_df = df.copy()

            # Dynamically identify all subgenome columns in the DataFrame
            subgenomes = [col for col in df.columns if 'subgenome' in col]

            # Iterate through each row of the DataFrame
            for i in range(len(new_df)):
                current_range = (new_df.at[i, 'start'], new_df.at[i, 'end'])

                # Initialize dictionaries to hold the previous non-'x' value for each subgenome
                previous_values = {sg: None for sg in subgenomes}
                
                # Update previous_values by scanning backwards from current position
                for j in range(i - 1, -1, -1):
                    for sg in subgenomes:
                        if new_df.at[j, sg] != 'x' and previous_values[sg] is None:
                            previous_values[sg] = new_df.at[j, sg]

                # Assign previous values to the current row where applicable
                for sg in subgenomes:
                    if new_df.at[i, sg] == 'x' and previous_values[sg] is not None:
                        # Ensure the previous value is not assigned to another subgenome in the same row
                        if all(new_df.at[i, other_sg] != previous_values[sg] for other_sg in subgenomes if other_sg != sg):
                            new_df.at[i, sg] = previous_values[sg]

                # Remove duplicates within the same row
                for sg1 in subgenomes:
                    for sg2 in subgenomes:
                        if sg1 != sg2 and new_df.at[i, sg1] == new_df.at[i, sg2] and new_df.at[i, sg1] != 'x':
                            # If there is a duplication, revert the current subgenome value to 'x'
                            new_df.at[i, sg2] = 'x'

            new_df.to_excel(f'cleaned_next_{key}.xlsx', index=False)
            # Store the cleaned dataframe
            cleaned_subgenome_dict[key] = new_df

        return cleaned_subgenome_dict
    
    def swap_x_with_values_in_dict(self, adjusted_subgenome_dict):
        """
        Purpose: This function swaps 'x' values with the next available value in the same row.
        Args: adjusted_subgenome_dict
        Returns: cleaned_subgenome_dict
        """
        cleaned_subgenome_dict = {}  # To store the processed dataframes

        for key, df in adjusted_subgenome_dict.items():
            # Copy the dataframe to avoid changing the original data
            new_df = df.copy()

            # Dynamically identify all subgenome columns in the DataFrame
            subgenomes = [col for col in df.columns if 'subgenome' in col]

            # Store initial 'x' presence for each subgenome to avoid swapping at the beginning
            initial_x_presence = {sg: df.at[0, sg] == 'x' for sg in subgenomes}

            # Iterate through each row of the DataFrame
            for index in range(len(new_df)):
                # Iterate through the set of subgenomes to identify if a swap is needed
                for i in range(len(subgenomes)-1):  # -1 because we'll compare each subgenome with the next one
                    if index == 0:
                        continue  # Skip the first row to avoid swapping at the beginning
                    current_subgenome = subgenomes[i]
                    next_subgenome = subgenomes[i + 1]

                    current_value = new_df.at[index, current_subgenome]
                    next_value = new_df.at[index, next_subgenome]

                    # Check conditions for swapping
                    if current_value != 'x' and next_value == 'x'and not initial_x_presence[current_subgenome]:
                        # Perform the swap
                        new_df.at[index, current_subgenome] = 'x'
                        new_df.at[index, next_subgenome] = current_value

            new_df.to_excel(f'cleaned_{key}.xlsx', index=False)
            # After processing all rows in the current dataframe, save it back to the cleaned_subgenome_dict
            cleaned_subgenome_dict[key] = new_df

        return cleaned_subgenome_dict
    
    def adjust_subgenome_for_continuity(self, adjusted_subgenome_dict, window_size=10):
        """
        Purpose: This function adjusts subgenome assignments to improve continuity.
        Args: adjusted_subgenome_dict, window_size
        Returns: adjusted_dict
        """
        adjusted_dict = {}  # To store the modified dataframes

        for key, df in adjusted_subgenome_dict.items():
            # Copy the dataframe to avoid changing the original data
            new_df = df.copy()

            # Dynamically identify all subgenome columns in the DataFrame
            subgenomes = [col for col in df.columns if 'subgenome' in col]

            # Iterate through each row of the DataFrame
            for index in range(len(new_df)):
                # Skip early rows to ensure we have a full window to look back on
                if index > window_size:
                    # Create a dictionary to hold counts of each chromosome name for each subgenome in the window
                    window_counts = {sg: {} for sg in subgenomes}

                    # Populate window_counts with chromosome name frequencies within the window for each subgenome
                    for i in range(index - window_size, index):
                        for sg in subgenomes:
                            chromo = new_df.at[i, sg]
                            if chromo != 'x':
                                window_counts[sg][chromo] = window_counts[sg].get(chromo, 0) + 1

                    # Iterate through subgenomes to see if a swap would improve continuity
                    for current_sg in subgenomes:
                        current_chromo = new_df.at[index, current_sg]
                        # Skip 'x' values and already dominant chromo
                        if current_chromo != 'x' and current_chromo not in window_counts[current_sg]:
                            # Look for a better fitting subgenome for the current chromosome based on window counts
                            for other_sg in subgenomes:
                                if current_sg != other_sg:
                                    # Check if the current chromosome is more common in the other subgenome within the window
                                    if window_counts[other_sg].get(current_chromo, 0) > window_counts[current_sg].get(current_chromo, 0):
                                        # Perform swap if the other subgenome's current value is less common in its own history
                                        other_chromo = new_df.at[index, other_sg]
                                        if other_chromo == 'x' or window_counts[other_sg].get(other_chromo, 0) < window_counts[current_sg].get(other_chromo, 0):
                                            # Swap the values
                                            new_df.at[index, current_sg], new_df.at[index, other_sg] = other_chromo, current_chromo
                                            break  # Only swap once for each row
            
            new_df.to_excel(f'final_{key}.xlsx', index=False)
            # After processing all rows in the current dataframe, update the dictionary with the adjusted DataFrame
            adjusted_dict[key] = new_df

        return adjusted_dict
    
    def reduce_redundant_continuous_assignments(self, df):
        """ 
        Purpose: This function reduces redundant continuous assignments in subgenome columns.
        Args: df
        Returns: new_df
        """
        # Copy the dataframe to avoid changing the original data
        new_df = df.copy()

        # Define the subgenome columns
        subgenome_columns = [col for col in new_df.columns if 'subgenome' in col]

        # Iterate through each subgenome column
        for col in subgenome_columns:
            # Initialize the last seen unique value and its position
            last_unique_value = None
            last_unique_position = None

            # Iterate through each row of the DataFrame
            for index, row in new_df.iterrows():
                current_value = row[col]

                # If the current value is different from the last unique value or is 'x',
                # then reset the last unique value and position
                if current_value != last_unique_value or current_value == 'x':
                    # If entering a new sequence, mark the last sequence as redundant except the first one
                    if last_unique_position is not None and index - last_unique_position > 1:
                        for redundant_index in range(last_unique_position + 1, index):
                            new_df.at[redundant_index, col] = 'x'
                    last_unique_value = current_value
                    last_unique_position = index
                # No else needed, since we handle marking redundant directly when we encounter a change

            # Clean up for sequences ending at the last element
            if last_unique_position is not None and len(new_df) - last_unique_position > 1:
                for redundant_index in range(last_unique_position + 1, len(new_df)):
                    new_df.at[redundant_index, col] = 'x'

        return new_df

    def adjust_all_subgenome_placements(self, adjusted_subgenome_dict):
        """
        Purpose: This function adjusts all subgenome placements to maintain continuity.
        Args: adjusted_subgenome_dict
        Returns: cleaned_subgenome_dict
        """
        cleaned_subgenome_dict = {}

        for key, df in adjusted_subgenome_dict.items():
            cleaned_df = self.reduce_redundant_continuous_assignments(df)
            cleaned_df.to_excel(f'finalist_{key}.xlsx', index=False)
            cleaned_subgenome_dict[key] = cleaned_df

        return cleaned_subgenome_dict


    def fill_x_with_previous_or_next_chr(self, df):
        """
        Purpose: This function fills 'x' values in subgenome columns with the previous or next chromosome name.
        Args: df
        Returns: df
        """
        # Iterate through each subgenome column
        for col in [c for c in df.columns if 'subgenome' in c]:
            # Iterate through each row in the column
            for i in range(len(df)):
                # If the current value is 'x'
                if df.loc[i, col] == 'x':
                    # Initialize a variable to hold the replacement value
                    replacement_chr = None
                    # Look backward for the previous non-'x' value that doesn't conflict with other subgenomes in the same row
                    for j in range(i - 1, -1, -1):  # Start from the previous row and move back
                        prev_chr = df.loc[j, col]
                        if prev_chr != 'x' and all(df.loc[i, other_col] != prev_chr for other_col in df.columns if 'subgenome' in other_col and other_col != col):
                            replacement_chr = prev_chr
                            break  # Stop after the first non-conflicting replacement
                    # If no previous non-'x' value was found, look ahead for the next non-'x' value
                    if replacement_chr is None:
                        for j in range(i + 1, len(df)):
                            next_chr = df.loc[j, col]
                            if next_chr != 'x' and all(df.loc[i, other_col] != next_chr for other_col in df.columns if 'subgenome' in other_col and other_col != col):
                                replacement_chr = next_chr
                                break  # Stop after the first non-conflicting replacement
                    # If a replacement value was found, use it to replace 'x'
                    if replacement_chr:
                        df.loc[i, col] = replacement_chr
        return df

    def combine_common_values_in_subgenomes(self, df):
        """
        Purpose: This function combines adjacent rows with common values in subgenomes.
        Args: df
        Returns: combined_df
        """
        combined_list = []
        previous_row = None

        for index, row in df.iterrows():
            if previous_row is not None and all(row[subgenome] == previous_row[subgenome] for subgenome in df.columns if 'subgenome' in subgenome):
                # Update the 'end' of the previous row in the combined list
                combined_list[-1]['end'] = row['end']
            else:
                combined_list.append(row.to_dict())
            previous_row = row

        return pd.DataFrame(combined_list)

    def merge_subgenome_blocks(self, adjusted_subgenome_dict):
        """
        Purpose: This function furge merges chromosome blocks based on the adjusted subgenome assignments.
        Args: adjusted_subgenome_dict
        Returns: merged_dfs
        """
        merged_dfs = {}  # Initialize a dictionary for the merged DataFrames

        for key, df in adjusted_subgenome_dict.items():
            normalized_df = pd.DataFrame(columns=['start', 'end'] + [col for col in df.columns if 'subgenome' in col])
            df.sort_values(by=['start', 'end'], inplace=True)
            
            unique_starts_ends = sorted(set(df['start'].tolist() + df['end'].tolist()))
            
            previous_end = None  # Keep track of the end of the previous block

            for i in range(len(unique_starts_ends) - 1):
                start = unique_starts_ends[i]
                end = unique_starts_ends[i + 1]

                # Adjust start if it is the same as the previous end
                if start == previous_end:
                    start += 1  # Increment start to avoid overlap

                # Skip if the adjustment makes start greater than end (no valid range)
                if start > end:
                    continue

                new_row_index = len(normalized_df)
                new_row = {'start': start, 'end': end}
                for subgenome in [col for col in df.columns if 'subgenome' in col]:
                    overlaps = df[(df['start'] <= start) & (df['end'] >= end)]
                    values = overlaps[subgenome].replace('x', pd.NA).dropna().unique().tolist()
                    new_row[subgenome] = ' '.join(values) if values else 'x'
                
                normalized_df.loc[new_row_index] = new_row
                previous_end = end  # Update previous_end for the next iteration
            
            subgenome_columns = [col for col in normalized_df.columns if 'subgenome' in col]
            normalized_df.dropna(axis=0, how='all', subset=subgenome_columns, inplace=True)
            normalized_df.replace({col: {'': 'x'} for col in subgenome_columns}, inplace=True)

            # Fill 'x' at the start of subgenomes with the next chr
            filled_normalized_df = self.fill_x_with_previous_or_next_chr(normalized_df)

            # Combine adjacent rows with common values in subgenomes
            combined_df = self.combine_common_values_in_subgenomes(filled_normalized_df)

            combined_df.to_excel(f'finalist_{key}.xlsx', index=False)
            merged_dfs[key] = combined_df

        return merged_dfs
    
    def resolve_subgenome_conflicts(self, df):
        """
        Purpose: This function resolves conflicts in subgenome assignments.
        Args: df
        Returns: df
        """
        # Assuming subgenomes are labeled as 'subgenome_1', 'subgenome_2', and 'subgenome_3'
        subgenome_cols = [col for col in df.columns if 'subgenome' in col]

        # Iterate through each row of the DataFrame
        for idx, row in df.iterrows():
            for col in subgenome_cols:
                # Split the values if there is more than one and remove extra spaces
                values = [val.strip() for val in str(row[col]).split() if val.strip() != 'x']

                if len(values) > 1:
                    # If there's more than one value, determine which one to keep
                    other_cols = [other_col for other_col in subgenome_cols if other_col != col]
                    other_values = [row[other_col] for other_col in other_cols if row[other_col] != 'x']

                    # Check if either of the values is present in another subgenome; if so, keep the other one
                    if values[0] in other_values and values[1] not in other_values:
                        df.at[idx, col] = values[1]  # Keep the second one
                    elif values[1] in other_values and values[0] not in other_values:
                        df.at[idx, col] = values[0]  # Keep the first one
                    else:
                        # If neither value appears in other subgenomes, default to the second value to maintain some continuity
                        df.at[idx, col] = values[1]  # Default to keep second
                elif len(values) == 1:
                    # If there's only one value (or none, but represented as a list from split), just keep it directly
                    df.at[idx, col] = values[0] if values else 'x'

        return df
    
    def resolve_subgenome_conflicts_driver(self, adjusted_subgenome_dict):
        """
        Purpose: This function resolves conflicts in subgenome assignments for all DataFrames.
        Args: adjusted_subgenome_dict
        Returns: adjusted_subgenome_dict
        """
        for key, original_df in adjusted_subgenome_dict.items():
            adjusted_df = self.resolve_subgenome_conflicts(original_df)
            # Now you can do further processing with adjusted_df or store it for later
            adjusted_subgenome_dict[key] = adjusted_df  # Replace the old dataframe with the adjusted one
            # Optionally, save to Excel or other formats as required
            adjusted_df.to_excel(f'adjusted_{key}.xlsx', index=False)

        return adjusted_subgenome_dict


    def remove_rows_with_duplicate_subgenomes_and_update_end(self, dfs):
        cleaned_dfs = {}  # Dictionary to store the cleaned DataFrames

        for key, original_df in dfs.items():
            print(f'Processing {key}:')

            new_df = original_df.copy()
            subgenome_columns = [col for col in new_df.columns if 'subgenome' in col]

            # Function to check for duplicates
            def has_duplicates(df):
                for _, row in df.iterrows():
                    if len(set(row[subgenome_columns])) < len(subgenome_columns):
                        return True
                return False

            # Repeat until no duplicates are found
            while has_duplicates(new_df):
                for index, row in new_df.iterrows():
                    duplicates = [col for col in subgenome_columns if list(row[subgenome_columns]).count(row[col]) > 1]

                    if duplicates:
                        for next_index in range(index + 1, len(new_df)):
                            next_row = new_df.iloc[next_index]
                            next_row_values = set(next_row[subgenome_columns])
                            if len(next_row_values) == len(subgenome_columns):
                                for col in duplicates:
                                    new_df.at[index, col] = next_row[col]
                                break

            new_df.to_excel(f'no_duplicate_{key}.xlsx', index=False)
            cleaned_dfs[key] = new_df

        return cleaned_dfs


    def check_duplicates_in_subgenomes(self, dfs):
        """
        This function checks for duplicates in subgenome assignments within each row.
        Args: 
            dfs: A dictionary of DataFrames with subgenome columns.
        Returns: 
            None
        """
        for key, original_df in dfs.items():
            print(f'Checking duplicates in {key}:')

            # Assuming subgenomes are labeled as 'subgenome_1', 'subgenome_2', 'subgenome_3', etc.
            subgenome_columns = [col for col in original_df.columns if 'subgenome' in col]

            # Initialize a counter for duplicate instances
            duplicate_counter = 0

            # Iterate through each row of the DataFrame
            for index, row in original_df.iterrows():
                values_seen = set()  # Track values seen in this row for different subgenomes
                duplicates = []  # Store pairs of subgenomes with duplicate values in this row

                # Check each subgenome's value in this row
                for sg in subgenome_columns:
                    value = row[sg]
                    if value != 'x' and value in values_seen:
                        # If the value (not 'x') has been seen before in this row, it's a duplicate
                        duplicates.append((sg, value))
                    values_seen.add(value)

                # If duplicates were found, print the row index, subgenome, and value
                if duplicates:
                    duplicate_counter += 1
                    duplicate_details = ', '.join([f'{sg}: {val}' for sg, val in duplicates])
                    print(f'Row {index} has duplicates: {duplicate_details}')

            # If no duplicates found, print that information
            if duplicate_counter == 0:
                print(f'No duplicates found in {key}.')
            print()  # Add a newline for better readability between different keys

    def get_densities_subgenome_placements(self, adjusted_sub, ploidy_numpy, chrsets, ploidy_numpy_overlap):
        """
        Purpose: This function calculates the densities of subgenome placements.
        Args: adjusted_sub, ploidy_numpy, chrsets
        Returns: density_dfs
        """
        density_dfs = {}

        # Determine total columns needed for each key based on chrsets
        total_columns_needed = {key: (int(value) + 1) * 2 for key, value in chrsets.items()}
        total_col = 0
        print("total_columns_needed", total_columns_needed)

        for key, value in chrsets.items():
            total_col += (int(value)) * 2
            print("total_col", total_col)

        for key, df in adjusted_sub.items():
            np_array = ploidy_numpy[key].astype(float)  # Ensure numpy array is of float type for calculations
            print(np_array)
            np_array_overlap = ploidy_numpy_overlap[key].astype(float)
            densities = []

            for index, row in df.iterrows():
                start, end = int(row['start']), int(row['end'])
                block_length = end - start + 1

                if block_length > 0:
                    block_densities = []

                    subgenomes = [col for col in df.columns if 'subgenome' in col]

                    # Calculate density only for relevant chromosome columns as indicated by subgenome data
                    for subgenome_col in subgenomes:  # Modify if there are more subgenomes
                        chromo_info = row[subgenome_col]

                        if chromo_info != 'x':
                            # Convert chromosome info to column index
                            if '.r' in chromo_info:  # Reverse strand
                                print(chromo_info)
                                if chromo_info != 'Un.r':
                                    col_idx = (int(chromo_info[3:-2]) - 1) * 2 + 1  # Convert to 0-based index and get reverse strand
                                    print("!Un.r",col_idx)
                                if chromo_info == 'Un.r':
                                    col_idx = total_col + 1
                                    print("Un.r",col_idx)

                            else:  # Forward strand
                                if chromo_info != 'Un' and chromo_info != 'overlap':
                                    col_idx = (int(chromo_info[3:]) - 1) * 2  # Convert to 0-based index
                                    print("forward",col_idx)
                                if chromo_info == 'Un':
                                    col_idx = total_col
                                    print("Un",col_idx)
                                if chromo_info == 'overlap':
                                    col_idx = total_col + 2
                                    print("overlap",col_idx)
                                    

                            # Ensure the column index is within the range for this ploidy level
                            if col_idx < total_columns_needed[key]:
                                if np_array.shape[1] > col_idx:
                                    sum_values = np.sum(np_array[start:end + 1, col_idx], axis=0)
                                    density = sum_values / block_length
                                else:
                                    sum_values = np.sum(np_array_overlap[start:end + 1, col_idx], axis=0)
                                    density = sum_values / block_length
                            else:
                                density = 0  # Default to zero if out of range
                        else:
                            density = 0  # Assign zero for 'x'

                        block_densities.append(density)

                    densities.append([start, end] + block_densities)

            # Convert the list of densities into a DataFrame
            column_names = ['start', 'end'] + [f'subgenome_{i+1}' for i in range(len(subgenomes))]  # Adjust according to actual subgenome count
            density_df = pd.DataFrame(densities, columns=column_names)
            density_dfs[key] = density_df

        # Save each density DataFrame as a separate sheet in an Excel file
        with pd.ExcelWriter('densities_subgenome_placements.xlsx') as writer:
            for key, df in density_dfs.items():
                df.to_excel(writer, sheet_name=key)

        return density_dfs
    

    def swap_max_density_to_subgenome1_for_all(self, adjusted_sub, density_dict, n_sub):
        """
        This function reassigns subgenomes based on their density rankings:
        the most dense one in subgenome_1, the second most dense one in subgenome_2,
        and the least dense one in subgenome_3.
        """
        updated_subgenome_dict = {}
        updated_density_dict = {}

        for key in adjusted_sub:
            df = adjusted_sub[key]
            densities = density_dict[key]

            assert len(df) == len(densities), f"DataFrames lengths do not match for {key}"

            updated_df = df.copy()
            updated_densities = densities.copy()

            #generalize the subgenome columns using n_sub
            subgenome_cols = [f'subgenome_{i}' for i in range(1, n_sub + 1)]

            for index, row in updated_df.iterrows():
                # Extract the corresponding densities row
                density_values = updated_densities.loc[index, subgenome_cols].values
                
                # Skip rows where the highest density is 1 or greater
                if np.abs(density_values).max() >= 1:
                    continue

                # Rank subgenomes based on their densities
                ranked_indices = np.argsort(-np.abs(density_values))  # argsort by descending absolute density
                ranked_subgenomes = [updated_df.at[index, subgenome_cols[i]] for i in ranked_indices]
                ranked_densities = [density_values[i] for i in ranked_indices]

                # Reassign the ranked subgenome and density values to their new order
                for new_idx, col in enumerate(subgenome_cols):
                    updated_df.at[index, col] = ranked_subgenomes[new_idx]
                    updated_densities.at[index, col] = ranked_densities[new_idx]

            # Store the updated dataframes back in the dictionaries
            updated_df.to_excel(f'new_adjusted_{key}.xlsx', index=False)
            updated_densities.to_excel(f'new_densities_{key}.xlsx', index=False)
            updated_subgenome_dict[key] = updated_df
            updated_density_dict[key] = updated_densities

        return updated_subgenome_dict, updated_density_dict

    
    
    def get_subgenomes(self, adjusted_subgenome, n_sub, chrsets, gene_id_matrix_df_, header, last_two_columns, last_column):
        """
        Purpose: This function assigns gene IDs to subgenomes based on the adjusted subgenome assignments.
        Args: adjusted_subgenome, n_sub, chrsets, gene_id_matrix_df_, header
        Returns: updated_subgenomes
        """
        updated_subgenomes = {}
        # Calculate the total number of columns needed for each chrset based on its value.
        total_columns_needed = {key: (int(value)) * 2 for key, value in chrsets.items()}

        # Create a cumulative total to keep track of the column ranges for each chrset.
        cumulative_total_columns = {}
        cumulative_sum = 0
        for key, value in total_columns_needed.items():
            cumulative_total_columns[key] = (cumulative_sum, cumulative_sum + value)
            cumulative_sum += value

        # Iterate over each key in adjusted_subgenome to map each subgenome to its chromosomes.
        for key, df in adjusted_subgenome.items():
            gene_id_matrix_df = gene_id_matrix_df_[key].reset_index(drop=True)
            gene_id_matrix_df_extended = pd.concat([gene_id_matrix_df, last_two_columns, last_column], axis=1)

            # Determine the chrset key corresponding to the current adjusted_subgenome key
            chrset_key = None
            for k in chrsets:
                if k in key:  # or any other logic that matches the keys
                    chrset_key = k
                    break

            if not chrset_key:
                continue  # Skip if no corresponding chrset is found

            # Apply the relevant column names based on chrset.
            start_col, end_col = cumulative_total_columns[chrset_key]
            relevant_columns = header[start_col:end_col]
            relevant_columns = relevant_columns + header[-3:]  # Add the last three columns to the relevant columns list
            gene_id_matrix_df_extended.columns = relevant_columns

            subgenome_cols = [col for col in df.columns if 'subgenome' in col]

            # Initialize with 'x' for all entries
            new_df = pd.DataFrame(index=range(0,len(gene_id_matrix_df_extended)-1), columns=subgenome_cols)
            new_df.fillna('x', inplace=True)
            
        #    Iterate through each row of the subgenome dataframe
            for _, row in df.iterrows():
                # For each subgenome, find the gene ID and update
                for i in range(1, n_sub + 1):
                    subgenome = f'subgenome_{i}'              

                    gene_rows = gene_id_matrix_df_extended.iloc[int(row['start']):int(row['end'])]

                    for gene_row_index in range(int(row['start']), int(row['end']) + 1):
                        if gene_row_index < len(gene_id_matrix_df_extended):
                            gene_row = gene_id_matrix_df_extended.iloc[gene_row_index]

                            new_df.at[gene_row_index, f'subgenome_{i}'] = gene_row[row[subgenome]]

            new_df.to_excel(f'new_subgenomes_{key}.xlsx', index=False)
            updated_subgenomes[key] = new_df

        return updated_subgenomes
    
    
    def adjust_subgenome_segmental_dup(self, updated_subgenome, adjusted_subgenome_results, collinear, n_sub, noisy_breakpoints, ploidy_status, chrsets, file, gene_id_matrix_df_):
        """
        Purpose: This function adds missing pieces of gene IDs to the subgenomes.
        Args: updated_subgenome, adjusted_subgenome_results, collinear, n_sub, noisy_breakpoints, ploidy_status, chrsets, file
        Returns: adjusted_subgenomes
        """
        # slice the collinear matrix based on the ploidy status
        # look for ploidy status: either diploid, tetraploid, hexaploid etc...
        
        adjusted_subgenomes = {}

        for key, df in updated_subgenome.items():
            gene_id_matrix_df = gene_id_matrix_df_[key].reset_index(drop=True)
            gene_id_matrix_df.columns = [f'chr{num}' for num in range(len(gene_id_matrix_df.columns))]

            valid_indices = noisy_breakpoints[key][noisy_breakpoints[key].iloc[:, -2]].index.tolist()  # Assuming last column indicates validity

            for i in valid_indices:
                # Ensure indices are within bounds
                if 0 <= i < len(df):
                    new_df = adjusted_subgenome_results[key].copy()
                    row = new_df.iloc[i]  # This is fine, 'row' is a Series now
                    density_columns = new_df.columns[2:-3]  # Access columns from the DataFrame, not from 'row'

                    # Iterate through each valid index                 
                    start_bp, end_bp = adjusted_subgenome_results[key].loc[i, ['start', 'end']]
                    used_genes = set()  # Initialize the set of used genes for each block
                    
                    # Iterate through each subgenome column
                    for sub_index in range(1, n_sub + 1):
                        # Identify candidate columns for replacement based on density
                        candidate_cols = []
                        for col in density_columns:
                            # Check if this column is not part of the current subgenomes and has a value greater than zero
                            if col not in row['subgenome_1':'subgenome_' + str(n_sub)].values and row[col] > 0:

                                candidate_cols.append(col)
                            
                        subgenome = f'subgenome_{sub_index}'
                        continuous_x_blocks = self.find_continuous_x_blocks(df[subgenome].iloc[int(start_bp):int(end_bp) + 1], 10)
                        
                        # Check and fill continuous blocks with 'x'
                        for block_start, block_end in continuous_x_blocks:
                            self.fill_in_genes_from_other_columns(block_start, block_end, df, gene_id_matrix_df, sub_index, n_sub, candidate_cols, used_genes)  # Pass the set of used genes
                            # Save the updated DataFrame for this key to an Excel file
            
            excel_filename = f'{key}_subgenomes.xlsx'
            df.to_excel(excel_filename, index=False)

            adjusted_subgenomes[key] = df  # Store adjusted subgenome data

        return adjusted_subgenomes
    
    