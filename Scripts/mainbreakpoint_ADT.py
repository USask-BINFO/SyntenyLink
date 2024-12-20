import Syntenylink_ADT as S
import numpy as np
import warnings
import pandas as pd
import sys



    # gap_threshold is the threshold to determine if there is a gap or not
    # min_block_length is the minimum length of the signal blocks

class BP(object):
    def __init__(self):
        """
        Purpose: Initialize the Main Breakpoint class object
        Args: None
        Returns: None
        """
        self.synteny = S.SyntenyLink()
    
    def main(self):
        """
        Purpose: This is the main function for the main breakpoint analysis. It will call the main function from the mainbreakpoint_ADT.py file and run the main function.
        Args: None
        Returns: ploidy_numpy, gap_threshold, minimum_block_length
        """

        # Get the input file as an argument (collinear file)
        input_file = sys.argv[sys.argv.index('-i') + 1]

        # Get the ploidy status as an argument
        ploidy_status = int(sys.argv[sys.argv.index('-s') + 1]) #if diploid 2, tetraploid 4, hexaploid 6, octaploid 8
        n_subgenomes = int(sys.argv[sys.argv.index('-n') + 1]) # 2, 3 ....
        gene_prefix = sys.argv[sys.argv.index('-p') + 1]

        #get the GT file default is None
        GT = None
        if '-gt' in sys.argv:
            GT = sys.argv[sys.argv.index('-gt') + 1]

        bed_file = sys.argv[sys.argv.index('-bed') + 1]

        # Get the chromosome sets as an argument
        if ploidy_status > 2:
            chrsets = {}  # Initialize an empty dictionary to hold the chromosome sets
            for i in range(1, (int(ploidy_status/2)) + 1):
                try:
                    chrsets[f'_ploidy_collinear_{i}'] = sys.argv[sys.argv.index('-chr' + str(i)) + 1]
                except ValueError as e:
                    print(f"Error: Argument '-chr{i}' not found.")
                    sys.exit(1)
                except IndexError as e:
                    print(f"Error: No value provided for '-chr{i}'.")
                    sys.exit(1)
        else:
            chrsets = {'_ploidy_collinear_1': sys.argv[sys.argv.index('-chr1') + 1]}

            # print(chrsets)

        if ploidy_status > 2:
            #get subgenome letters for chromosomes
            subgenome_letters = {}
            for i in range(1, (int(ploidy_status/2)) + 1):
                try:
                    subgenome_letters[f'_subgenome_{i}'] = sys.argv[sys.argv.index('-sub' + str(i)) + 1]
                except ValueError as e:
                    print(f"Error: Argument '-sub{i}' not found.")
                    sys.exit(1)
                except IndexError as e:
                    print(f"Error: No value provided for '-sub{i}'.")
                    sys.exit(1)
        else:
            subgenome_letters = {}
            # have two N s inside the dictionary
            subgenome_letters = {'_subgenome_1': 'N', '_subgenome_2': 'N'}

        
        if n_subgenomes < 2:
            print("Error: Invalid value for '-n'.")
            sys.exit(1)

        # Open the input file
        file = open(input_file, "r")

        # Convert the collinear file to a dataframe
        matrix = self.synteny.create_matrix(file)

        ref_gene_id = self.synteny.extract_reference_gene_id(matrix)

        # Create a collinear matrix
        collinear_all = self.synteny.create_collinear_all(matrix)

        # Update the collinear matrix with overlap genes
        collinear_all_mod = self.synteny.update_overlap_genes(collinear_all, bed_file,gene_prefix)

        #get overlap genes not added to collinear
        empty_df = self.synteny.overlap_genes_get(collinear_all, chrsets, bed_file, ploidy_status, subgenome_letters)

        # Split the collinear matrix without last three columns
        collinear = self.synteny.split_collinear(collinear_all_mod)

        # Get the header from the input matrix
        header_filtered, header = self.synteny.get_header(matrix)

        # convert 'x' into 0
        collinear = self.synteny.replace_x(collinear)


        # convert genes in forward strand into 1 and genes in reverse strand into -1
        collinear = self.synteny.replace_forward_reverse(header_filtered, collinear)

        # slice the collinear matrix based on the ploidy status
        # look for ploidy status: either diploid, tetraploid, hexaploid etc...
        ploidy_df = self.synteny.slice_collinear(ploidy_status, chrsets, collinear)
        ploidy_df_overlap = self.synteny.slice_collinear_merge_overlap(ploidy_status, chrsets, collinear)

        # convert the dataframe to a numpy array
        ploidy_numpy = self.synteny.convert_matrix_to_numpy(ploidy_df)
        ploidy_numpy_overlap = self.synteny.convert_matrix_to_numpy(ploidy_df_overlap)
     
        # get the gap statistics and the gap threshold  
        gap_stats, gap_threshold = self.synteny.gapthresh(ploidy_numpy)


        # get the block statistics and the minimum block length
        block_stats, minimum_block_length = self.synteny.min_blocklength(ploidy_numpy)


        # plot the gap statistics and the block statistics
        # self.synteny.plot_thresholds(gap_stats, block_stats)
        
        return n_subgenomes, ploidy_numpy, matrix, gap_threshold, minimum_block_length, file, ploidy_status, chrsets, header, header_filtered, GT, empty_df, ref_gene_id, bed_file, ploidy_numpy_overlap, gene_prefix, subgenome_letters
    

    def find_gaps(self, col, rows, gap_threshold):
        """
        Purpose: This function finds the gaps in each column of the collinear matrix
        Args: col, rows
        Returns: gaps
        """

        indices = np.where(np.abs(col) == 1)[0]

        if len(indices) == 0:
            return []  # No gaps if no '1's are found
        d = np.diff(indices)-1


        # if the difference is bigger than the gap threshold then store them
        initial_gaps = np.where(d > gap_threshold)[0]

        # number of gaps
        r = len(initial_gaps)
        # calculate gaps
        gaps = np.zeros((r, 2))
        # find the indices of beginning and end of gaps
        gaps[:, 0] = indices[initial_gaps] + 1
        gaps[:, 1] = indices[initial_gaps + 1] - 1
        # fix the end of gap index for the last gap
        if r > 0 and indices[-1] + 1 < rows:
            gaps = np.append(gaps, np.array([[indices[-1] + 1, rows]]), axis=0)
        return gaps
    

    def gap_calculation(self, ploidy_numpy, gap_threshold, minimum_block_length, n_sub):
        """
        Purpose: This function calculates the gaps in the collinear matrix.
        Returns: Dictionary of breakpoints for each key in ploidy_numpy.
        """
        j = 0  # Index for the gap_threshold and minimum_block_length arrays
        breakpoints = {}  # To store breakpoints for each ploidy set

        for key, np_array in ploidy_numpy.items():
            row, column = np_array.shape
            all_gaps = []  # Collect all gaps for each column in the current array

            # Find gaps in each column and store as 2D arrays in all_gaps
            for i in range(column):
                gaps = self.find_gaps(np_array[:, i], row, gap_threshold[j])
                all_gaps.append(np.array(gaps).reshape(-1, 2))  # Ensure 2D structure even if empty

            break_point_indices = np.zeros(5000)

            # Find the best columns to start with
            col_ind = np.argpartition(np.sum(np_array[:minimum_block_length[j], :], axis=0), -(n_sub-1))[-(n_sub-1):]

            # Find the first gap
            if n_sub == 2:
                i = int(all_gaps[col_ind[0]][0, 0])
            else:
                i = min(int(all_gaps[col_ind[k]][0, 0]) for k in range(n_sub - 1))

            break_point_indices[0] = i  # Initialize the first breakpoint
            count = 1  # Start counting breakpoints

            # Loop to find the rest of the breakpoints
            while i <= row - 1:
                i = int(i)  # Ensure i is integer for indexing
                col_ind = np.argpartition(
                    np.sum(np_array[i: min(row, i + minimum_block_length[j]), :], axis=0), -n_sub
                )[-n_sub:]

                i_vals = []
                for k in range(n_sub):
                    current_gaps = all_gaps[col_ind[k]]
                    if current_gaps.size == 0:
                        i_vals.append(row)
                    else:
                        i_vals.append(current_gaps[current_gaps[:, 0] > i][0, 0] if (current_gaps[:, 0] > i).any() else row)
                
                i = min(i_vals)
                
                # Add breakpoint if it meets gap threshold condition
                if abs(i - break_point_indices[count - 1]) > gap_threshold[j]:
                    break_point_indices[count] = i
                    count += 1

            # Resize array to include only found breakpoints
            break_point_indices = np.resize(break_point_indices, count)

            # Add the last row as the final breakpoint if it’s not already included
            if break_point_indices[-1] != row:
                break_point_indices = np.append(break_point_indices, row)

            breakpoints[key] = break_point_indices
            j += 1  # Move to the next set of parameters

        return breakpoints

    

    def small_gap_calculation(self, ploidy_numpy, gap_threshold, minimum_block_length, n_sub):
        """
        Purpose: This function calculates the smaller gaps in the main breakpoints
        Args: ploidy_numpy, gap_threshold, minimum_block_length, n_sub
        Returns: break_point_indices

        """
        breakpoints = []

        for key, value in ploidy_numpy.items():
            np_array = value  # Assuming it's already a numpy array
            row, column = np_array.shape

            all_gaps = []  # Collect all gap lengths for the current key

            for i in range(column):
                if len(gap_threshold) == 0 or len(minimum_block_length) == 0:
                    gap_threshold = []
                    gap_threshold.append(5)
                    minimum_block_length = []
                    minimum_block_length.append(5)
                gaps = self.find_gaps(np_array[:, i], row, gap_threshold[0])
                if isinstance(gaps, list):  # Convert list to numpy array if necessary
                    gaps = np.array(gaps)
                if gaps.size > 0:  # Add only non-empty gaps
                    all_gaps.append(gaps)

            break_point_indices = np.zeros(5000)
            col_ind = np.argpartition(np.sum(np_array[:minimum_block_length[0], :], axis=0), -(n_sub-1))[-(n_sub-1):]

            # Check for the first gap based on the best columns
            i = row  # Default if no gaps are found
            for k in range(min(len(col_ind), n_sub)):  # Ensure k is within bounds
                if k < len(all_gaps) and col_ind[k] < len(all_gaps) and len(all_gaps[col_ind[k]]) > 0:
                    i = min(i, int(all_gaps[col_ind[k]][0][0]))  # Update i to the minimum starting gap

            break_point_indices[0] = i
            count = 1

            # Find the rest of the breakpoints
            while i < row:
                i = int(i)
                col_ind = np.argpartition(np.sum(np_array[i:min(row, i + minimum_block_length[0]), :], axis=0), -n_sub)[-n_sub:]
                next_gap_values = [row] * n_sub  # Default to row if no further gaps

                for k in range(min(len(col_ind), n_sub)):
                    if k < len(all_gaps) and col_ind[k] < len(all_gaps) and len(all_gaps[col_ind[k]]) > 0:
                        filtered_gaps = all_gaps[col_ind[k]][all_gaps[col_ind[k]][:, 0] > i]
                        if len(filtered_gaps) > 0:
                            next_gap_values[k] = np.min(filtered_gaps[:, 0])

                new_i = min(next_gap_values)
                if abs(new_i - break_point_indices[count - 1]) > gap_threshold[0]:
                    break_point_indices[count] = new_i
                    count += 1

                i = new_i  # Update i for the next iteration

            break_point_indices = np.resize(break_point_indices, count)

            # Update the last breakpoint if necessary
            if break_point_indices[-1] != row:
                break_point_indices = np.append(break_point_indices, row)

            breakpoints = break_point_indices

        return breakpoints

    def small_bp(self, ploidy_numpy, key_bp, breakpoints, small_breakpoints, start_bp_index):
        """
        Purpose: This function calculates the smaller breakpoints
        Args: ploidy_numpy, key_bp, breakpoints, small_breakpoints, start_bp_index
        Returns: new_breakpoints

        """
        new_breakpoints = {}  # Dictionary to store new breakpoints

        for key, np_array in ploidy_numpy.items():
            if key == key_bp:
                if small_breakpoints is not None:
                    # get the index of value equal to start_bp_index in breakpoints list
                    i = breakpoints.index(start_bp_index)

                    if 0 <= i < len(breakpoints):
                        # Prepare a list to collect all new breakpoint positions
                        new_bp_positions = []

                        # Track the original breakpoint to be replaced or where to start inserting
                        original_bp_position = breakpoints[i]


                        for j, small_bp in enumerate(small_breakpoints):
                            # Calculate new breakpoint position relative to start_bp_index
                            new_bp = start_bp_index + small_bp
                            if i != len(breakpoints)-2 and i < len(breakpoints) - 1:
                                if len(small_breakpoints) == 1:
                                    # For the first small breakpoint, we replace the original breakpoint
                                    print(breakpoints, i)
                                    breakpoints[i+1] = new_bp
                                else:
                                    if j == len(small_breakpoints) - 1:
                                        # For the last small breakpoint, we replace the original breakpoint
                                        breakpoints[i+1] = new_bp
                                    else:
                                        # For additional small breakpoints, just prepare to insert them
                                        new_bp_positions.append(new_bp)
                            elif i == len(breakpoints)-2:
                                if new_bp+5 < int(breakpoints[i+1]):
                                    
                                    if len(small_breakpoints) == 1:
                                        # For the first small breakpoint, we replace the original breakpoint
                                        breakpoints[i+1] = new_bp
                                    else:
                                        if j == len(small_breakpoints) - 1:
                                            # For the last small breakpoint, we replace the original breakpoint
                                            breakpoints[i+1] = new_bp
                                        else:
                                            # For additional small breakpoints, just prepare to insert them
                                            new_bp_positions.append(new_bp)

                        # Insert new breakpoints while maintaining the order
                        for new_bp in sorted(new_bp_positions):  # Ensure new breakpoints are processed in order
                            # Find the correct position to insert new breakpoints
                            # This ensures that even after replacing the first breakpoint, others are inserted correctly
                            insert_position = next((idx for idx, val in enumerate(breakpoints) if val > new_bp), len(breakpoints))
                            breakpoints.insert(insert_position, new_bp)

                        #sort the breakpoints
                        breakpoints.sort()
                        # Convert the list to a numpy array                       
                        new_breakpoints[key] = np.array(breakpoints)  # Convert back to a numpy array if 
                        return new_breakpoints


    def get_densities_bp(self, break_point_indices, ploidy_numpy):
        """
        Purpose: This function calculates the densities of each block and saves them as an Excel file.
        Args: break_point_indices
        Returns: density of each block
        """
        # Get the ploidy status, ploidy_numpy, gap_threshold, and minimum_block_length
        # n_sub, ploidy_numpy, gap_threshold, minimum_block_length = self.main()
        # Get the break point indices
        density = {}

        # Initialize an empty dictionary to store DataFrames for each key
        density_dfs = {}

        # for each key in the ploidy_numpy dictionary
        for key, value in ploidy_numpy.items():
            np_array = value  # Assuming it's already a numpy array
            row, column = np_array.shape

            # calculate densities of each block
            b = len(break_point_indices[key])  # number of breakpoints
            densities = np.zeros((b, column+2))  # densities for each blocks
            for i in range(b):
                if i == 0:  # for block one
                    #add start and end of the block start as 0 and end as break_point_indices[i]
                    densities[i,0] = 0
                    densities[i,1] = break_point_indices[key][i]

                    densities[i, 2:] = (
                        np.sum(np_array[: int(break_point_indices[key][i]), :], axis=0)
                        / int(break_point_indices[key][i])
                    ).astype(float)

                else:  # for all else
                    #add start and end of the block start as break_point_indices[i-1] and end as break_point_indices[i]
                    densities[i,0] = break_point_indices[key][i-1]+1
                    densities[i,1] = break_point_indices[key][i]

                    densities[i, 2:] = (
                        np.sum(
                            np_array[int(break_point_indices[key][i - 1])+1: int(break_point_indices[key][i]), :],
                            axis=0,
                        )
                        / (int(break_point_indices[key][i]) - int(break_point_indices[key][i - 1]))
                    ).astype(float)
                    
            density[key] = densities
            # Convert the numpy array to a pandas DataFrame
            density_dfs[key] = pd.DataFrame(densities, columns=[f'Column_{i+1}' for i in range(column+2)])

        # Save each density DataFrame as a separate sheet in an Excel file
        with pd.ExcelWriter('densities.xlsx') as writer:
            for key, df in density_dfs.items():
                df.to_excel(writer, sheet_name=str(key))
        
        return density


    def get_noisy_breakpoints(self, densities, n_sub):
        """
        Purpose: This function calculates the noisy breakpoints. A noisy breakpoint is a breakpoint that has columns with densities that exceed the number of subgenomes (n_sub). 
        Args: densities
        Returns: noisy breakpoints
        """
        # n_sub, ploidy_numpy, gap_threshold, minimum_block_length = self.main()
        subgenomes_adjusted = {}

        for key, value in densities.items():
            # Ensure value is treated as a DataFrame
            if isinstance(value, pd.DataFrame):
                df_temp = value[:,2:]
            else:  # Assuming value is a numpy array if not already a DataFrame
                df_temp = pd.DataFrame(value[:,2:])

            # Initialize a DataFrame to hold pairwise aggregated data
            col_names = pd.DataFrame(index=df_temp.index)

            # Pairwise aggregate forward and reverse strands
            for col in range(0, df_temp.shape[1] - 1, 2):  # Adjusting for DataFrame column count
                forward_strand = df_temp.iloc[:, col]
                reverse_strand = df_temp.iloc[:, col + 1].abs()  # Assuming reverse strand values are negative
                # Aggregate by summing the absolute values of forward and reverse strands
                col_names[f'pair{col // 2}'] = forward_strand.abs() + reverse_strand

            # Count non-zero elements in each aggregated column
            col_names["Non_zero"] = col_names.apply(lambda x: (x != 0).sum(), axis=1)

            # Check rows with non-zero counts exceeding the number of subgenomes
            col_names["Exceeds_n_sub"] = col_names["Non_zero"] > n_sub
            #print the TRUE values in the Exceeds_n_sub column
            # print(col_names[col_names["Exceeds_n_sub"] == True])
            col_names["Less_than_n_sub"] = col_names["Non_zero"] < n_sub

            # Now col_names contains a boolean column "Exceeds_n_sub" indicating if the non-zero count exceeds n_sub
            subgenomes_adjusted[key] = col_names
        
        return subgenomes_adjusted


    def smaller_breakpoints(self, breakpoints, noisy_breakpoints, n_sub, ploidy_numpy):
        """
        Purpose: This function calculates the smaller breakpoints. A smaller breakpoint is a breakpoint that has columns with densities that exceed the number of subgenomes (n_sub).
        Args: breakpoints, noisy_breakpoints
        Returns: smaller breakpoints
        """
        # n_sub, ploidy_numpy, gap_threshold, minimum_block_length = self.main()
        flattened_dict = {}

        for key, np_array in ploidy_numpy.items():
            if key in noisy_breakpoints and key in breakpoints:
                valid_indices = noisy_breakpoints[key][noisy_breakpoints[key].iloc[:, -2] == True].index.tolist()
            
                updated_breakpoints = list(breakpoints[key])  # Work on a copy to maintain original just in case

                for i, index in enumerate(valid_indices):
                    if i > 0:
                        updated_breakpoints = list(updated_breakpoints[key])  # Work on a copy to maintain original just in case
                    if index < len(updated_breakpoints):
                        start_bp_index = int(breakpoints[key][index - 1] + 1 if index > 0 else 0)
                        end_bp_index = int(breakpoints[key][index])
                        # Extract and process the segment
                        segment = np_array[start_bp_index:end_bp_index + 1, :]
                        if segment.size > 0:
                            df_segments = pd.DataFrame(np.vstack(segment))
                            df_numpy = df_segments.to_numpy()
                            df_numpy_dict = {key: df_numpy}
                            # Recalculate stats and breakpoints based on this segment
                            gap_stats, gap_threshold = self.synteny.gapthresh_minblock(df_numpy_dict)
                            block_stats, minimum_block_length = self.synteny.min_blocklength_minblock(df_numpy_dict)
                            new_breakpoints = self.small_gap_calculation(df_numpy_dict, gap_threshold, minimum_block_length, n_sub)
                            
                            # Update breakpoints for the next iteration
                            
                            updated_breakpoints = self.small_bp(df_numpy_dict, key, updated_breakpoints, new_breakpoints, int(breakpoints[key][index - 1]))

                # Once all indices are processed, update the original breakpoints
                breakpoints[key] = updated_breakpoints

        for outer_key, inner_dict in breakpoints.items():
            for inner_key, value in inner_dict.items():
                # Assign the value to the outer key
                flattened_dict[outer_key] = value

        return flattened_dict


    def get_subgenomes_unstranded(self, noisy_sbp, densities, n_sub, ):
        """
        Purpose: This function assigns subgenomes to each block based on the densities of each block.
        Args: densities, n_sub
        Returns: subgenomes_result
        """
        warnings.filterwarnings("ignore")
        subgenomes_result = {}

        for key, density_matrix in noisy_sbp.items():
            # Copy the original density matrix except for the last two columns
            subgenome_data = density_matrix.iloc[:, :-3].copy()
            #add first and second column of densities[key] to the subgenome_data
            subgenome_data.insert(0, 'start', densities[key][:, 0])
            subgenome_data.insert(1, 'end', densities[key][:, 1])


            # Prepare a DataFrame to hold subgenome assignments, using the same index as the original
            for sub in range(1, n_sub+1):
                subgenome_data[f'subgenome_{sub}'] = pd.NA  # Initialize subgenome columns with NA

            for i in range(len(density_matrix)):
                block_data = density_matrix.iloc[i, :-3]  # Work with data excluding the last two columns

                current_subgenomes = []
                for _ in range(n_sub):
                    max_col = block_data.idxmax()
                    current_subgenomes.append(max_col)
                    block_data = block_data.drop(max_col)

                num_nonzero = density_matrix.iloc[i, -3]
                if i > 0 and num_nonzero < n_sub:
                    # Fill missing subgenomes with those from the previous block
                    for j in range(n_sub - num_nonzero):
                        current_subgenomes[n_sub - 1 - j] = subgenome_data.iloc[i - 1][f'subgenome_{n_sub - j}']

                # Assign current subgenomes to this block in the new DataFrame
                for idx, sg in enumerate(current_subgenomes, start=1):
                    subgenome_data.at[i, f'subgenome_{idx}'] = sg

            # Once all subgenome assignments are completed, check for initial rows needing updates
            first_valid_row = next((idx for idx, row in density_matrix.iterrows() if row[-2] >= n_sub), None)
            if first_valid_row is not None:
                for j in range(first_valid_row):
                    if density_matrix.iloc[j, -3] < n_sub:
                        # Update these initial rows with the first valid subgenomes
                        for sub in range(1, n_sub+1):
                            subgenome_data.at[j, f'subgenome_{sub}'] = subgenome_data.at[first_valid_row, f'subgenome_{sub}']

            # Store the final DataFrame for this key in the results
            subgenomes_result[key] = subgenome_data

        return subgenomes_result
    
    def adjust_subgenome_unstranded(self, subgenome_results, n_sub):
        """
        Purpose: This function adjusts the subgenome assignments to maintain continuity.
        Args: subgenome_results, n_sub
        Returns: adjusted_subgenome_results
        """
        adjusted_subgenome_results = {}

        for key, df in subgenome_results.items():
            for i in range(1, len(df)):
                # Extract current and previous subgenome assignments as lists
                current_subgenomes = [df.at[i, f'subgenome_{j}'] for j in range(1, n_sub + 1)]
                prev_subgenomes = [df.at[i - 1, f'subgenome_{j}'] for j in range(1, n_sub + 1)]

                # Continuity correction logic
                # Compare each current subgenome to the previous ones to find the best match
                for j, curr_sg in enumerate(current_subgenomes):
                    best_match_idx = j  # Assume the best match is the current position
                    min_distance = n_sub  # Max possible distance is n_sub-1

                    for k, prev_sg in enumerate(prev_subgenomes):
                        distance = abs(j - k)  # Calculate "distance" between positions
                        if curr_sg == prev_sg and distance < min_distance:
                            best_match_idx = k
                            min_distance = distance

                    # If a better match is found, swap subgenomes to maintain continuity
                    if best_match_idx != j:
                        # Swap in current_subgenomes list
                        current_subgenomes[j], current_subgenomes[best_match_idx] = current_subgenomes[best_match_idx], current_subgenomes[j]
                        # Update the dataframe
                        df.at[i, f'subgenome_{j+1}'] = current_subgenomes[j]
                        df.at[i, f'subgenome_{best_match_idx+1}'] = current_subgenomes[best_match_idx]

            for i in range(1, len(df)):
                # Ensure no duplicate subgenome assignments in the same row
                assigned_subgenomes = set()
                for p in range(1, n_sub + 1):
                    current_sg = df.at[i, f'subgenome_{p}']
                    if current_sg in assigned_subgenomes:
                        # If duplicate, replace with the value from the previous row
                        df.at[i, f'subgenome_{p}'] = df.at[i-1, f'subgenome_{p}']
                        # Check if this replacement causes another duplicate
                        if df.at[i, f'subgenome_{p}'] in assigned_subgenomes:
                            # Find the previous subgenome value that is not in the current assigned set
                            for q in range(1, n_sub + 1):
                                if df.at[i-1, f'subgenome_{q}'] not in assigned_subgenomes:
                                    df.at[i, f'subgenome_{p}'] = df.at[i-1, f'subgenome_{q}']
                                    break
                    assigned_subgenomes.add(df.at[i, f'subgenome_{p}'])

            excel_filename = f'{key}_subgenome_assignment.xlsx'
            df.to_excel(excel_filename, index=False)
            adjusted_subgenome_results[key] = df

        return adjusted_subgenome_results



    def get_subgenomes_stranded(self, adjusted_subgenome_results, file, n_sub, ploidy_status, chrsets, bed_file, gene_prefix):
        """
        Purpose: update the subgenome assignments based on the strand information and regenerate the updated subgenome assignments by extracting gene ids from the collinear matrix using block information.
        Args: adjusted_subgenome_results, file, n_sub, ploidy_status, chrsets, bed_file
        Returns: updated_subgenomes
        """

        # slice the collinear matrix based on the ploidy status
        # look for ploidy status: either diploid, tetraploid, hexaploid etc...
        gene_id_matrix_df_ = self.synteny.slice_collinear(ploidy_status, chrsets, pd.DataFrame(self.synteny.split_collinear(self.synteny.update_overlap_genes(self.synteny.create_collinear_all(self.synteny.create_matrix(file)), bed_file, gene_prefix))))
        # Dictionary to store the updated subgenomes for each key
        updated_subgenomes = {}

        for key, subgenome_df in adjusted_subgenome_results.items():
            gene_id_matrix_df = gene_id_matrix_df_[key].reset_index(drop=True)
            gene_id_matrix_df.columns = [f'chr{num}' for num in range(len(gene_id_matrix_df.columns))]
            
            # Create a new DataFrame to store subgenome gene IDs for this key
            new_df = pd.DataFrame(index=range(0,len(gene_id_matrix_df)-1), columns=[f'subgenome_{i}' for i in range(1, n_sub + 1)])
            new_df.fillna('x', inplace=True)

            # Last selected strand with a gene ID for any subgenome
            last_selected_strand = None

            # Iterate through each row of the subgenome dataframe
            for _, row in subgenome_df.iterrows():
                # For each subgenome, find the gene ID and update
                for i in range(1, n_sub + 1):
                    subgenome = f'subgenome_{i}'
                    pair_index = int(row[subgenome].replace('pair', ''))  # Extract index from pair name
                    matrix_col_idx_forward = ((pair_index) * 2)
                    matrix_col_idx_reverse = matrix_col_idx_forward + 1

                    gene_rows = gene_id_matrix_df.iloc[int(row['start']):int(row['end'])]

                    for gene_row_index in range(int(row['start']), int(row['end']) + 1):
                        if gene_row_index < len(gene_id_matrix_df):
                            gene_row = gene_id_matrix_df.iloc[gene_row_index]

                            forward_gene = gene_row[matrix_col_idx_forward] != 'x'
                            reverse_gene = gene_row[matrix_col_idx_reverse] != 'x'

                            # If only one of the strands has a gene, select it and update the last selected strand
                            if forward_gene and not reverse_gene:
                                new_df.at[gene_row_index, f'subgenome_{i}'] = gene_row[matrix_col_idx_forward]
                                last_selected_strand = 'forward'
                            elif reverse_gene and not forward_gene:
                                new_df.at[gene_row_index, f'subgenome_{i}'] = gene_row[matrix_col_idx_reverse]
                                last_selected_strand = 'reverse'
                            # If both strands have genes, use the last selected strand for prioritization
                            elif forward_gene and reverse_gene:
                                if last_selected_strand == 'forward':
                                    new_df.at[gene_row_index, f'subgenome_{i}'] = gene_row[matrix_col_idx_forward]
                                else:
                                    new_df.at[gene_row_index, f'subgenome_{i}'] = gene_row[matrix_col_idx_reverse]
                            else:
                                new_df.at[gene_row_index, f'subgenome_{i}'] = 'x'
            
            updated_subgenomes[key] = new_df
            # Save the updated subgenome DataFrame to an Excel file
            excel_filename = f'{key}_subgenome_assignment_updated.xlsx'
            new_df.to_excel(excel_filename, index=False)

        return updated_subgenomes


    def find_continuous_x_blocks(self, series, min_size):
        """
        Purpose: This function finds continuous 'x' blocks in a series.
        Args: series, min_size
        Returns: blocks
        """
        is_x = (series == 'x')  # Reset index for continuous, zero-based indexing
        blocks = []
        start_index = None
        continuous = False  # This flag will keep track of continuity

        for i, x in enumerate(is_x):
            if x == True:
                if start_index is None:
                    start_index = i # Mark the start of a new 'x' block
                continuous = True  # We are inside a continuous 'x' block
            
            else:
                # We exit an 'x' block, so check if it was continuous and long enough
                if continuous and start_index is not None and (i - start_index) >= min_size:
                    blocks.append((start_index, i-1))  # Add the continuous block to the list
                # Reset for the next block
                start_index = None
                continuous = False
        
        # Handle case where the last characters in the series are 'x' and form a valid block
        if continuous and start_index is not None and (series.index[-1] - start_index) >= min_size:
            blocks.append((start_index, series.index[-1] - 1))
            
        return blocks

    def fill_in_genes_from_other_columns(self, block_start, block_end, updated_subgenome, gene_id_matrix, sub_index, n_sub, potential_cols, used_genes):
        """
        Purpose: This function fills in missing gene IDs in a block from other columns.
        Args: block_start, block_end, updated_subgenome, gene_id_matrix, sub_index, n_sub, potential_cols, used_genes
        Returns: None
        """

        for idx in range(block_start, block_end + 1):
            for col in potential_cols:
                pair_idx = int(col.replace('pair', ''))
                forward_idx, reverse_idx = pair_idx * 2, pair_idx * 2 + 1
                forward_gene, reverse_gene = gene_id_matrix.iloc[idx, forward_idx], gene_id_matrix.iloc[idx, reverse_idx]
                overlap_gene_forward = gene_id_matrix.iloc[idx, -2]
                overlap_gene_reverse = gene_id_matrix.iloc[idx, -1]

                if forward_gene != 'x' and forward_gene not in used_genes:
                    updated_subgenome.iloc[idx, sub_index - 1] = forward_gene
                    used_genes.add(forward_gene)  # Update the set of used genes
                    break  # Ensure only the first unique gene ID is used for the block
                elif reverse_gene != 'x' and reverse_gene not in used_genes:
                    updated_subgenome.iloc[idx, sub_index - 1] = reverse_gene
                    used_genes.add(reverse_gene)  # Update the set of used genes
                    break
                elif overlap_gene_forward != 'x' and overlap_gene_forward not in used_genes:
                    updated_subgenome.iloc[idx, sub_index - 1] = overlap_gene_forward
                    used_genes.add(overlap_gene_forward)
                    break
                elif overlap_gene_reverse != 'x' and overlap_gene_reverse not in used_genes:
                    updated_subgenome.iloc[idx, sub_index - 1] = overlap_gene_reverse
                    used_genes.add(overlap_gene_reverse)
                    break

        return updated_subgenome   


    def adjust_subgenome_stranded(self, updated_subgenome, adjusted_subgenome_results, collinear, n_sub, noisy_breakpoints, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns):
        """
        Purpose: This function adds missing pieces of gene IDs to the subgenomes.
        Args: updated_subgenome, adjusted_subgenome_results, collinear, n_sub, noisy_breakpoints, ploidy_status, chrsets, file
        Returns: adjusted_subgenomes
        """
        # slice the collinear matrix based on the ploidy status
        # look for ploidy status: either diploid, tetraploid, hexaploid etc...
        
        adjusted_subgenomes_stranded = {}
        used_genes = set()  # Initialize the set of used genes for each block

        for key, df in updated_subgenome.items():
            gene_id_matrix_df = gene_id_matrix_df_[key].reset_index(drop=True)
            gene_id_matrix_df_extended = pd.concat([gene_id_matrix_df, last_two_columns], axis=1)
            gene_id_matrix_df_extended.columns = [f'chr{num}' for num in range(len(gene_id_matrix_df_extended.columns))]

            valid_indices = noisy_breakpoints[key][noisy_breakpoints[key].iloc[:, -2]].index.tolist()  # Assuming last column indicates validity
            
            for i in valid_indices:
                # Ensure indices are within bounds
                if 0 <= i < len(df):
                    new_df = adjusted_subgenome_results[key].copy()
                    row = new_df.iloc[i]  # This is fine, 'row' is a Series now
                    density_columns = new_df.columns[2:-3]  # Access columns from the DataFrame, not from 'row'

                    # Iterate through each valid index                 
                    start_bp, end_bp = adjusted_subgenome_results[key].loc[i, ['start', 'end']]
                    used_genes = set(df.iloc[:, :n_sub].values.flatten())
                    # Iterate through each subgenome column
                    for sub_index in range(1, n_sub + 1):
                        # Identify candidate columns for replacement based on density
                        candidate_cols = []
                        for col in density_columns:
                            # Check if this column is not part of the current subgenomes and has a value greater than zero
                            if col not in row['subgenome_1':'subgenome_' + str(n_sub)].values and row[col] > 0:

                                candidate_cols.append(col)
                            
                        subgenome = f'subgenome_{sub_index}'
                        continuous_x_blocks = self.find_continuous_x_blocks(df[subgenome].iloc[int(start_bp):int(end_bp) + 1], 2)
                        
                        # Check and fill continuous blocks with 'x'
                        for block_start, block_end in continuous_x_blocks:
                            df = self.fill_in_genes_from_other_columns(block_start, block_end, df, gene_id_matrix_df_extended, sub_index, n_sub, candidate_cols, used_genes)  # Pass the set of used genes
                            # Save the updated DataFrame for this key to an Excel file
            
            excel_filename = f'{key}_subgenomes_dup.xlsx'
            df.to_excel(excel_filename, index=False)

            adjusted_subgenomes_stranded[key] = df  # Store adjusted subgenome data

        return adjusted_subgenomes_stranded
    
    def fill_in_genes_from_other_columns_segmental(self, block_start, block_end, updated_subgenome, gene_id_matrix, sub_index, n_sub, potential_cols, used_genes, empty_df):
        """
        Purpose: This function fills in missing gene IDs in a block from other columns.
        Args: block_start, block_end, updated_subgenome, gene_id_matrix, sub_index, n_sub, potential_cols, used_genes
        Returns: None
        """

        for idx in range(block_start, block_end + 1):
            genes_in_other_subgenomes = set()  # Initialize the set of genes in other subgenomes
            for col in potential_cols:
                pair_idx = int(col.replace('pair', ''))
                forward_idx, reverse_idx = pair_idx * 2, pair_idx * 2 + 1
                forward_gene, reverse_gene = gene_id_matrix.iloc[idx, forward_idx], gene_id_matrix.iloc[idx, reverse_idx]
                overlap_gene_forward = gene_id_matrix.iloc[idx, -2]
                overlap_gene_reverse = gene_id_matrix.iloc[idx, -1]
                overlap_idx = empty_df.iloc[idx, pair_idx]

                for sub in range(0, n_sub):
                    if forward_gene != 'x' and updated_subgenome.iloc[idx, sub] == forward_gene and forward_gene not in genes_in_other_subgenomes:
                        genes_in_other_subgenomes.add(forward_gene)
                    if reverse_gene != 'x' and updated_subgenome.iloc[idx, sub] == reverse_gene and reverse_gene not in genes_in_other_subgenomes:
                        genes_in_other_subgenomes.add(reverse_gene)
                    if overlap_gene_forward != 'x' and updated_subgenome.iloc[idx, sub] == overlap_gene_forward and overlap_gene_forward not in genes_in_other_subgenomes:
                        genes_in_other_subgenomes.add(overlap_gene_forward)
                    if overlap_gene_reverse != 'x' and updated_subgenome.iloc[idx, sub] == overlap_gene_reverse and overlap_gene_reverse not in genes_in_other_subgenomes:
                        genes_in_other_subgenomes.add(overlap_gene_reverse)
                    if overlap_idx != 'x' and updated_subgenome.iloc[idx, sub] == overlap_idx and overlap_idx not in genes_in_other_subgenomes:
                        genes_in_other_subgenomes.add(overlap_idx)


                if overlap_idx != 'x' and overlap_idx not in used_genes and overlap_idx not in genes_in_other_subgenomes:
                    updated_subgenome.iloc[idx, sub_index - 1] = overlap_idx
                    used_genes.add(overlap_idx)
                    break
                if forward_gene != 'x' and forward_gene not in used_genes and forward_gene not in genes_in_other_subgenomes:
                    updated_subgenome.iloc[idx, sub_index - 1] = forward_gene
                    used_genes.add(forward_gene)  # Update the set of used genes
                    break  # Ensure only the first unique gene ID is used for the block
                if reverse_gene != 'x' and reverse_gene not in used_genes and reverse_gene not in genes_in_other_subgenomes:
                    updated_subgenome.iloc[idx, sub_index - 1] = reverse_gene
                    used_genes.add(reverse_gene)  # Update the set of used genes
                    break
                if overlap_gene_forward != 'x' and overlap_gene_forward not in used_genes and overlap_gene_forward not in genes_in_other_subgenomes:
                    updated_subgenome.iloc[idx, sub_index - 1] = overlap_gene_forward
                    used_genes.add(overlap_gene_forward)
                    break
                if overlap_gene_reverse != 'x' and overlap_gene_reverse not in used_genes and overlap_gene_reverse not in genes_in_other_subgenomes:
                    updated_subgenome.iloc[idx, sub_index - 1] = overlap_gene_reverse
                    used_genes.add(overlap_gene_reverse)
                    break
        return updated_subgenome
                
    def adjust_subgenome_stranded_segmental(self, updated_subgenome, adjusted_subgenome_results, collinear, n_sub, noisy_breakpoints, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns, empty_df, ref_gene_id):
        """
        Purpose: This function adds missing pieces of gene IDs to the subgenomes.
        Args: updated_subgenome, adjusted_subgenome_results, collinear, n_sub, noisy_breakpoints, ploidy_status, chrsets, file
        Returns: adjusted_subgenomes
        """
        # slice the collinear matrix based on the ploidy status
        # look for ploidy status: either diploid, tetraploid, hexaploid etc...
        total_columns_needed = {key: int(value) for key, value in chrsets.items()}

         # Create a cumulative total to keep track of the column ranges for each chrset.
        cumulative_total_columns = {}
        cumulative_sum = 0
        if ploidy_status > 2:
            for key, value in total_columns_needed.items():

                cumulative_total_columns[key] = (cumulative_sum, cumulative_sum + value)
                cumulative_sum += value
            adjusted_subgenomes = {}
        
        elif ploidy_status == 2:
            for key, value in total_columns_needed.items():
                cumulative_total_columns[key] = (cumulative_sum, cumulative_sum + value)
            adjusted_subgenomes = {}


        for key, df in updated_subgenome.items():
            # get the value of the key from total_columns_needed
            total_columns = cumulative_total_columns[key]
            # reset column index starting from 0
            empty_df_filtered = empty_df.iloc[:, total_columns[0]:total_columns[1]]
            # Reset column index by assigning a new list of column names
            empty_df_filtered.columns = range(empty_df_filtered.shape[1])

            gene_id_matrix_df = gene_id_matrix_df_[key].reset_index(drop=True)
            gene_id_matrix_df_extended = pd.concat([gene_id_matrix_df, last_two_columns], axis=1)
            gene_id_matrix_df_extended.columns = [f'chr{num}' for num in range(len(gene_id_matrix_df_extended.columns))]

            less_than_n_sub_indices = noisy_breakpoints[key][noisy_breakpoints[key].iloc[:, -1]].index.tolist()  # Indices with less than n_sub non-zero values

            for i in less_than_n_sub_indices:
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
                            candidate_cols.append(col)
                            
                        subgenome = f'subgenome_{sub_index}'
                        continuous_x_blocks = self.find_continuous_x_blocks(df[subgenome].iloc[int(start_bp):int(end_bp) + 1], 2)
                        
                        # Check and fill continuous blocks with 'x'
                        for block_start, block_end in continuous_x_blocks:
                            df = self.fill_in_genes_from_other_columns_segmental(block_start, block_end, df, gene_id_matrix_df_extended, sub_index, n_sub, candidate_cols, used_genes, empty_df_filtered)
                            # Save the updated DataFrame for this key to an Excel file
            #add ref_gene_ids to the updated_subgenome as the first column no key in ref_gene_id. it's just a list of gene_ids
            df.insert(0, 'ref_gene_id', ref_gene_id)
            
            excel_filename = f'{key}_subgenomes.xlsx'
            df.to_excel(excel_filename, index=False)

            adjusted_subgenomes[key] = df  # Store adjusted subgenome data

        return adjusted_subgenomes
        
    

    # def accuracy_gene_id(self, GT_file_path, n_sub, df_subgenome_dict):
    #     """
    #     Purpose: This function calculates the accuracy of subgenome assignments compared to the ground truth.
    #     Args: GT_file_path, n_sub, df_subgenome_dict, new_chr_placement
    #     Returns: results
    #     """
    #     # Read the ground truth DataFrame from the provided file path
    #     df_groundtruth = pd.read_excel(GT_file_path)

    #     # Initialize a dictionary to store the results for each subgenome
    #     results = {}

    #     # If df_subgenome_dict is a dictionary, extract a specific DataFrame; otherwise, use it directly
    #     # This assumes that you need to iterate over each key in the dictionary
    #     key_num = 0
    #     for key, df_subgenome in df_subgenome_dict.items():  # Now iterating through dictionary items
    #         key_num += 1
    #         # Initialize accuracy data for this particular subgenome DataFrame
    #         results[key] = []

    #         for j in range(n_sub):
    #             # Initialize metrics for accuracy calculation
    #             exact_matches = 0
    #             total_genes = 0
    #             missing_genes = 0

    #             # Compare gene placements between ground truth and the subgenome DataFrame
    #             for i, row_groundtruth in df_groundtruth.iterrows():
    #                 if key_num == 1:
    #                     first_letter = 'A'
    #                 elif key_num == 2:
    #                     first_letter = 'C'
    #                 if i < len(df_subgenome):  # Ensure index is within the range of df_subgenome
    #                     gene_id_groundtruth = row_groundtruth[f"{first_letter}_subgenome{j+1}"]
    #                     if gene_id_groundtruth != 'x':
    #                         total_genes += 1
    #                         row_subgenome = df_subgenome.iloc[i]  # df_subgenome used directly
    #                         gene_id_subgenome = row_subgenome[f"subgenome_{j+1}"]
    #                         if gene_id_groundtruth == gene_id_subgenome:
    #                             exact_matches += 1
    #                         else:
    #                             missing_genes += 1

    #             # Calculate and store the overlap percentage for this subgenome
    #             overlap_percentage = exact_matches / total_genes if total_genes > 0 else 0
    #             results[key].append({
    #                 'subgenome': j+1,
    #                 'exact_matches': exact_matches,
    #                 'total_genes': total_genes,
    #                 'missing_genes': missing_genes,
    #                 'overlap_percentage': overlap_percentage
    #             })
    #             print(results[key])

    #     # Now, results contain the accuracy data for each subgenome in each key of df_subgenome_dict
    #     return results
    
    # def accuracy_gene_id_new(self, GT_file_path, n_sub, df_subgenome, key_num):
    #     """
    #     Purpose: This function calculates the accuracy of subgenome assignments compared to the ground truth.
    #     Args: GT_file_path, n_sub, df_subgenome_dict, new_chr_placement
    #     Returns: results
    #     """
    #     # Read the ground truth DataFrame from the provided file path
    #     df_groundtruth = pd.read_excel(GT_file_path)

    #     # Initialize a dictionary to store the results for each subgenome
    #     results = {}

    #     # If df_subgenome_dict is a dictionary, extract a specific DataFrame; otherwise, use it directly
    #     # This assumes that you need to iterate over each key in the dictionary
    #     # Initialize accuracy data for this particular subgenome DataFrame
    #     results[key_num] = []

    #     for j in range(n_sub):
    #         # Initialize metrics for accuracy calculation
    #         exact_matches = 0
    #         total_genes = 0
    #         missing_genes = 0

    #         # Compare gene placements between ground truth and the subgenome DataFrame
    #         for i, row_groundtruth in df_groundtruth.iterrows():
    #             if key_num == '_ploidy_collinear_1':
    #                 first_letter = 'A'
    #             elif key_num == '_ploidy_collinear_2':
    #                 first_letter = 'C'
    #             if i < len(df_subgenome):  # Ensure index is within the range of df_subgenome
    #                 gene_id_groundtruth = row_groundtruth[f"{first_letter}_subgenome{j+1}"]
    #                 if gene_id_groundtruth != 'x':
    #                     total_genes += 1
    #                     row_subgenome = df_subgenome.iloc[i]  # df_subgenome used directly
    #                     gene_id_subgenome = row_subgenome[f"subgenome_{j+1}"]
    #                     if gene_id_groundtruth == gene_id_subgenome:
    #                         exact_matches += 1
    #                     else:
    #                         missing_genes += 1

    #         # Calculate and store the overlap percentage for this subgenome
    #         overlap_percentage = exact_matches / total_genes if total_genes > 0 else 0
    #         results[key_num].append({
    #             'subgenome': j+1,
    #             'exact_matches': exact_matches,
    #             'total_genes': total_genes,
    #             'missing_genes': missing_genes,
    #             'overlap_percentage': overlap_percentage
    #         })
    #         print(results[key_num])

    #     # Now, results contain the accuracy data for each subgenome in each key of df_subgenome_dict
    #     return results
    
    
    # def accuracy(self, GT_file_path, n_sub, df_subgenome_dict, new_chr_placement):
    #     """
    #     Purpose: This function calculates the accuracy of subgenome assignments compared to the ground truth.
    #     Args: GT_file_path, n_sub, df_subgenome_dict, new_chr_placement
    #     Returns: results
    #     """
    #     # Read the ground truth DataFrame from the provided file path
    #     df_groundtruth = pd.read_excel(GT_file_path)

    #     # Initialize a dictionary to store the results for each subgenome
    #     results = {}

    #     # If df_subgenome_dict is a dictionary, extract a specific DataFrame; otherwise, use it directly
    #     # This assumes that you need to iterate over each key in the dictionary
    #     for key, df_subgenome in df_subgenome_dict.items():  # Now iterating through dictionary items
    #         # Extract the first value which is not 'x' from the first subgenome of new_chr_placement
    #         first_chr_name = next((val for val in new_chr_placement[key]['subgenome_1'] if val != 'x'), 'x')
    #         first_letter = first_chr_name[0] if first_chr_name != 'x' and isinstance(first_chr_name, str) else 'x'

    #         # Initialize accuracy data for this particular subgenome DataFrame
    #         results[key] = []

    #         for j in range(n_sub):
    #             # Initialize metrics for accuracy calculation
    #             exact_matches = 0
    #             total_genes = 0
    #             missing_genes = 0

    #             # Compare gene placements between ground truth and the subgenome DataFrame
    #             for i, row_groundtruth in df_groundtruth.iterrows():
    #                 if i < len(df_subgenome):  # Ensure index is within the range of df_subgenome
    #                     gene_id_groundtruth = row_groundtruth[f"{first_letter}_subgenome{j+1}"]
    #                     if gene_id_groundtruth != 'x':
    #                         total_genes += 1
    #                         row_subgenome = df_subgenome.iloc[i]  # df_subgenome used directly
    #                         gene_id_subgenome = row_subgenome[f"subgenome_{j+1}"]
    #                         if gene_id_groundtruth == gene_id_subgenome:
    #                             exact_matches += 1
    #                         else:
    #                             missing_genes += 1

    #             # Calculate and store the overlap percentage for this subgenome
    #             overlap_percentage = exact_matches / total_genes if total_genes > 0 else 0
    #             results[key].append({
    #                 'subgenome': j+1,
    #                 'exact_matches': exact_matches,
    #                 'total_genes': total_genes,
    #                 'missing_genes': missing_genes,
    #                 'overlap_percentage': overlap_percentage
    #             })
    #             print(results[key])

    #     # Now, results contain the accuracy data for each subgenome in each key of df_subgenome_dict
    #     return results
    
    
    # def find_missing_genes(self, GT_file_path, df_subgenome_dict):
    #     """
    #     Purpose: This function finds missing genes in the subgenome assignments compared to the ground truth.
    #     Args: GT_file_path, df_subgenome_dict
    #     Returns: missing_genes_dict
    #     """
    #     # Read the ground truth DataFrame from the provided file path
    #     df_groundtruth = pd.read_excel(GT_file_path)
        
    #     # Initialize a dictionary to store the missing genes for each key
    #     missing_genes_dict = {}

    #     # Split the ground truth DataFrame into parts corresponding to each key
    #     gt_keys = list(df_subgenome_dict.keys())
    #     # gt_split = {gt_keys[0]: df_groundtruth.iloc[:, :3], gt_keys[1]: df_groundtruth.iloc[:, 3:6]}
    #     gt_split = {gt_keys[0]: df_groundtruth.iloc[:, :3]}
    #     # Iterate over each key in the df_subgenome_dict
    #     for key in gt_keys:
    #         # Convert the respective part of the ground truth DataFrame to a set of unique genes, excluding 'x'
    #         gt_genes = set(gt_split[key].replace('x', pd.NA).stack().unique())

    #         # Check if the key exists in df_subgenome_dict and convert that DataFrame to a set of unique genes, excluding 'x'
    #         if key in df_subgenome_dict:
    #             df_genes = set(df_subgenome_dict[key].replace('x', pd.NA).stack().unique())
    #         else:
    #             df_genes = set()  # If key doesn't exist in df_subgenome_dict, use an empty set for comparison

    #         # Find the difference between ground truth genes and DataFrame genes for the current key
    #         missing_genes = gt_genes - df_genes

    #         print(len(missing_genes))
    #         # Store the missing genes for this key
    #         missing_genes_dict[key] = missing_genes

    #     # Return the dictionary containing the missing genes for each key
    #     return missing_genes_dict