import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

class SyntenyLink(object):
    def __init__(self):
        """
        Purpose:
            Initialize a SyntenyLink object instance
        Args:
            None
        Returns:
            None
        """
        self.__gaps = []
        self.__blocklength = []
        self.__matrix = []
        self.__collinear=[]
        self.__collinear_all = []
        self.__header_filtered = []
        self.__header = []
        self.__ploidy_data = {}
        self.__ploidy_data_overlap = {}
    
    def create_matrix(self, file):
        """
        Purpose:
            Create a matrix from the input file
        Args:
            file: The input file
        Returns:
            A matrix
                    
        """
        self.__matrix = []  # Ensure matrix is cleared for each call
        for name in file:
            name = name.rstrip()
            rows = re.split(r'\t+', name)
            self.__matrix.append(rows)  # Directly append the split line
        return self.__matrix
    
    def extract_reference_gene_id(self, matrix):
        """
        Purpose:
            Extract the reference gene ID from the input matrix
        Args:
            matrix: The input matrix
        Returns:
            The reference gene ID
        """
        #extract all rows first column from the matrix except the first row
        reference_gene_id = [row[0] for row in matrix[1:]]
        return reference_gene_id

    
    def split_collinear(self,matrix):
        """
        Purpose:
            Filter the collinear matrix to remove the last three columns.
        Args:
            matrix: The input matrix, a list of lists.
        Returns:
            A new matrix with the last three columns removed from each row.
        """
        
        if isinstance(matrix, pd.DataFrame):
            return matrix.iloc[:, :-3]
        elif isinstance(matrix, list):
            # Ensure all elements of matrix are lists
            if all(isinstance(row, list) for row in matrix):
                return [row[:-3] for row in matrix]
            else:
                raise TypeError("All elements of matrix must be lists.")
        else:
            raise TypeError("Matrix must be a list of lists or a DataFrame.")
    
    def split_collinear_all(self,matrix):
        """
        Purpose:
            Filter the collinear matrix to remove the last three columns.
        Args:
            matrix: The input matrix, a list of lists.
        Returns:
            A new matrix with the last three columns removed from each row.
        """
        
        if isinstance(matrix, pd.DataFrame):
            return matrix.iloc[:, :]
        elif isinstance(matrix, list):
            # Ensure all elements of matrix are lists
            if all(isinstance(row, list) for row in matrix):
                return [row[:] for row in matrix]
            else:
                raise TypeError("All elements of matrix must be lists.")
        else:
            raise TypeError("Matrix must be a list of lists or a DataFrame.")



    def create_collinear_all(self, matrix):
        """
        Purpose:
            Create a collinear matrix from the input matrix
        Args:
            matrix: The input matrix
        Returns:
            A collinear matrix
        """
        for index in range (1,len(matrix)):
            self.__collinear_all.append(matrix[index][1:])
            
        return self.__collinear_all
    

    def update_overlap_genes(self, data, bedfile, gene_prefix):
        """
        Purpose:
            Update the columns replacing 'x' with overlap genes in the collinear matrix
        Args:
            data: The collinear matrix (can be a list of lists or a DataFrame)
            bedfile: Path to the BED file
            gfffile: Path to the GFF file
        Returns:
            The collinear matrix with updated overlap genes
        """
        # Convert data to DataFrame if it's a list of lists
        if isinstance(data, list):
            data = pd.DataFrame(data)

        # Read the BED file and store gene info in a dictionary
        gene_info_bed = {}
        with open(bedfile, 'r') as bf:
            for line in bf:
                line = line.rstrip().split("\t")
                gene_info_bed[line[3]] = line[0]  # Store chromosome number and position by gene ID

        for row_index, row in data.iterrows():
            overlap_gene = row.iloc[-1]  # The last element is the overlap gene

            if overlap_gene == 'x':
                continue

            overlap_gene_str = str(overlap_gene)
            overlap_prefix = overlap_gene_str.split(gene_prefix)[0]

            found_gene_id = False
            target_col = None

            # Check in the BED file info
            if overlap_gene in gene_info_bed:
                target_col = gene_info_bed[overlap_gene]

            if target_col:
                for col_index, gene_id in enumerate(row[:-1]):  # Exclude the overlap gene
                    gene_id_str = str(gene_id)
                    if gene_id_str.startswith(overlap_prefix):
                        target_col = col_index
                        found_gene_id = True
                        break  # Break the inner loop once a matching gene ID is found

                # Replace 'x' in the identified target column with the gene ID from the overlap column, if applicable
                if found_gene_id and row.iloc[target_col] == 'x':
                    row.iloc[target_col] = overlap_gene_str

        return data

    
    def overlap_genes_get(self, data, chrsets, bedfile, ploidy_status, subgenome_letters):
        """
        Purpose:
            Get the genes in overlap column which do not have a match in the collinear matrix
        Args:
            data: The collinear matrix (can be a list of lists or a DataFrame)
            chrsets: A dictionary containing chromosome sets information
            bedfile: Path to the bed file
            subgenome_letters: Dictionary mapping subgenome indices to their letter codes
            ploidy_status: Ploidy level to determine column adjustments
        Returns:
            A DataFrame with updated overlap genes and an Excel sheet for genes with no match.
        """
        total_columns_needed = {key: int(value) for key, value in chrsets.items()}

        # Convert data to DataFrame if it's a list of lists
        if isinstance(data, list):
            data = pd.DataFrame(data)

        num_rows = len(data)
        num_cols = sum(total_columns_needed.values())
        empty_df = pd.DataFrame('x', index=range(num_rows), columns=range(num_cols))

        # Read the bed file and store gene info in a dictionary
        gene_info = {}
        with open(bedfile, 'r') as bf:
            for line in bf:
                line = line.rstrip().split("\t")
                gene_info[line[3]] = line[0]  # Store chromosome by gene ID

        # Generalized: get the first subgenome letter dynamically and exclude it
        subgenome_values = list(subgenome_letters.values())
        first_subgenome_letter = subgenome_values[0]  # Use the first subgenome letter for adjustments

        for row_index, row in data.iterrows():
            overlap_gene = row.iloc[-1]  # The last element is the overlap gene

            if overlap_gene == 'x':
                continue

            overlap_gene_str = str(overlap_gene)

            if overlap_gene in gene_info:
                target_col = gene_info[overlap_gene]
                overlap_letter = target_col[0]
                
                # Handle "Chr" and replace it with the relevant subgenome letter
                if "Chr" in target_col:
                    target_col = target_col.replace("Chr", subgenome_letters.get(f'_subgenome_{overlap_letter}', overlap_letter))

                try:
                    # Extract chromosome number and adjust column based on subgenome
                    overlap_number = int(target_col[1:])
                    if overlap_letter != first_subgenome_letter and ploidy_status > 2:
                        cumulative_sum = 0
                        for key, value in total_columns_needed.items():
                            if key == overlap_letter:
                                break
                            cumulative_sum += value
                        cumulative_sum += overlap_number
                        col = cumulative_sum
                    else:
                        col = overlap_number - 1

                    # Update or log the overlap gene
                    if row.iloc[col] == 'x':
                        row.iloc[col] = overlap_gene_str
                    elif row.iloc[col] != overlap_gene_str:
                        empty_df.at[row_index, col] = overlap_gene_str

                except ValueError:
                    print(f"Unexpected chromosome format or non-numeric part in target_col: {target_col}")

        # Save unmatched genes to Excel
        empty_df.to_excel("overlap_genes.xlsx", index=False)
        return empty_df

    def replace_x(self, collinear):
        """
        Purpose:
            Replace 'x' with 0 in the collinear matrix
        Args:
            collinear: The collinear matrix
        Returns:
            The collinear matrix with 'x' replaced with 0
        """
        collinear = pd.DataFrame(collinear)
        collinear = collinear.replace(r'^x$', 0, regex=True)
        return collinear

    def get_header(self, matrix):
        """
        Purpose:
            Get the header from the input matrix
        Args:
            matrix: The input matrix
        Returns:
            The header
        """
        self.__header = matrix[0][1:]
        self.__header_filtered = matrix[0][1:-3]
        return self.__header_filtered, self.__header

    def replace_forward_reverse(self, header, collinear):
        """
        Purpose:
            Replace genes in the forward strand with 1 and genes in the reverse strand with -1
        Args:
            header: The header
            collinear: The collinear matrix 
        Returns:
            The collinear matrix with genes in the forward strand replaced with 1 and genes in the reverse strand replaced with -1
        """
        for row in range(len(collinear)):
            for column in range(len(header)):
                cell_value = collinear.iloc[row, column]
                if cell_value != 0 and re.match(r'^[a-zA-Z]', str(cell_value)):
                    if len(header[column].split(".")) > 1:
                        if header[column].split(".")[len(header[column].split(".")) - 1] == 'r':
                            collinear.iloc[row, column] = -1
                        else:
                            collinear.iloc[row, column] = 1
                    else:
                        collinear.iloc[row, column] = 1
        return collinear
    
    def slice_collinear(self, ploidy_status, chrsets, collinear):
        """
        Purpose:
            Slice the collinear matrix based on the ploidy status
        Args:
            ploidy_status: The ploidy status (diploid = 2, tetraploid = 4, hexaploid = 6, octaploid = 8, etc.)
            chrsets: The number of chromosome sets
            collinear: The collinear matrix
        Returns:
            The sliced collinear matrix in a dictionary
        """
        if ploidy_status == 2:
            chr_key = f'_ploidy_collinear_1'

            start_bound = 0  # Initialize starting bound

            upper_bound = start_bound + (int(chrsets[f'_ploidy_collinear_1'])) * 2

            # Slice the DataFrame from start_bound to upper_bound
            self.__ploidy_data[chr_key] = collinear.iloc[:, start_bound:upper_bound]

        elif ploidy_status > 2:

              # Initialize an empty dictionary to store the slices
            start_bound = 0  # Initialize starting bound

            for chr in range(1, (int(ploidy_status/2)) + 1):
                chr_key = f'_ploidy_collinear_{chr}'
                # Dynamically determine the current upper bound for the slice
                upper_bound = start_bound + (int(chrsets[f'_ploidy_collinear_{chr}'])) * 2
                
                # Slice the DataFrame from start_bound to upper_bound
                self.__ploidy_data[chr_key] = collinear.iloc[:, start_bound:upper_bound]
                
                
                # Update start_bound for the next slice to start right after the current upper_bound
                start_bound = upper_bound
        return self.__ploidy_data
    
    def slice_collinear_merge_overlap(self, ploidy_status, chrsets, collinear):
        """
        Purpose:
            Slice the collinear matrix and add last three columns including overlap column based on the ploidy status
        Args:
            ploidy_status: The ploidy status (diploid = 2, tetraploid = 4, hexaploid = 6, octaploid = 8, etc.)
            chrsets: The number of chromosome sets
            collinear: The collinear matrix
        Returns:
            The sliced collinear matrix in a dictionary
        """
        if ploidy_status == 2:
            chr_key = f'_ploidy_collinear_1'

            start_bound = 0  # Initialize starting bound

            upper_bound = start_bound + int(chrsets[f'_ploidy_collinear_1']) * 2

            # Slice the DataFrame from start_bound to upper_bound
            self.__ploidy_data_overlap[chr_key] = collinear.iloc[:, start_bound:upper_bound]
            #add last three columns at the end
            self.__ploidy_data_overlap[chr_key] = pd.concat([self.__ploidy_data_overlap[chr_key], collinear.iloc[:, -3:]], axis=1)
            
        elif ploidy_status > 2:

              # Initialize an empty dictionary to store the slices
            start_bound = 0  # Initialize starting bound

            for chr in range(1, (int(ploidy_status/2)) + 1):
                chr_key = f'_ploidy_collinear_{chr}'
                # Dynamically determine the current upper bound for the slice
                upper_bound = start_bound + int(chrsets[f'_ploidy_collinear_{chr}']) * 2
                
                # Slice the DataFrame from start_bound to upper_bound
                self.__ploidy_data_overlap[chr_key] = collinear.iloc[:, start_bound:upper_bound]
                #add last three columns at the end
                self.__ploidy_data_overlap[chr_key] = pd.concat([self.__ploidy_data_overlap[chr_key], collinear.iloc[:, -3:]], axis=1)
                
                # Update start_bound for the next slice to start right after the current upper_bound
                start_bound = upper_bound

        return self.__ploidy_data_overlap
    
    def convert_matrix_to_numpy(self, ploidy_data):
        """
        Purpose:
            Convert the ploidy data to a numpy array if it is not already
        Args:
            ploidy_data: The ploidy data
        Returns:
            The ploidy data as a numpy array
        """
        for key, value in ploidy_data.items():
            if not isinstance(value, np.ndarray):  # Check if the value is already a NumPy array
                ploidy_data[key] = value.to_numpy()  # Convert to NumPy array if it is a DataFrame
            row, column = ploidy_data[key].shape
        return ploidy_data


    def find_gaps(self, col):
        """
        Purpose:
            Find the gaps in a column of the numpy array and return gap lengths
        Args:
            col: A column from a numpy array
        Returns:
            A list of gap lengths in the column
        """
        indices = np.where(np.abs(col) == 1)[0]
        if len(indices) == 0:
            return []  # No gaps if no '1's are found
        gaps = np.diff(indices) - 1  # Subtract 1 to get actual gap lengths between '1's
        return gaps[gaps > 0].tolist()  # Filter out non-positive gaps and convert to list
    
    def gapthresh(self, ploidy_data):
        """
        Purpose:
            Calculate gap threshold for each key in the ploidy_data dictionary
        Args:
            None
        Returns:
            A list of gap thresholds for each key in the ploidy_data dictionary
        """
        gap_stats = {}  # Dictionary to store gap statistics for each key

        for key, value in ploidy_data.items():
            np_array = value  # Assuming it's already a numpy array
            row, column = np_array.shape
            all_gaps = []  # Collect all gap lengths for the current key

            for i in range(column):
                gaps = self.find_gaps(np_array[:, i])
                all_gaps.extend(gaps)  # Collect gaps across all columns for the current key

            # Filter out gaps below a certain threshold, e.g., 20
            filtered_gaps = [gap for gap in all_gaps if gap > 20]

            # If there are no gaps meeting the criteria for the current key
            if not filtered_gaps:
                gap_stats[key] = {"No significant gaps": None}
                continue  # Skip to the next key

            unique_gaps, gap_counts = np.unique(filtered_gaps, return_counts=True)
            
            # Calculate frequency statistics for the current key
            max_gap_length = unique_gaps[np.argmax(gap_counts)]
            mean_gap_length = np.mean(unique_gaps)
            median_gap_length = np.median(unique_gaps)

            max_gap_freq = np.max(gap_counts)
            mean_gap_freq = np.mean(gap_counts)
            median_gap_freq = np.median(gap_counts)

            # Store the statistics in the gap_stats dictionary under the current key
            gap_stats[key] = {
                "unique_gaps": unique_gaps,
                "gap_counts": gap_counts,
                "max_gap_length": max_gap_length,
                "mean_gap_length": mean_gap_length,
                "median_gap_length": median_gap_length,
                "max_gap_freq": max_gap_freq,
                "mean_gap_freq": mean_gap_freq,
                "median_gap_freq": median_gap_freq
            }

            self.__gaps.append(gap_stats[key]["max_gap_length"])

        return gap_stats, self.__gaps
    
    def gapthresh_minblock(self, ploidy_data):
        """
        Purpose:
            Calculate gap threshold for each key in the ploidy_data dictionary
        Args:
            None
        Returns:
            A list of gap thresholds for each key in the ploidy_data dictionary
        """
        gap_stats = {}  # Dictionary to store gap statistics for each key
        self.__gaps = []

        for key, value in ploidy_data.items():
            np_array = value  # Assuming it's already a numpy array
            row, column = np_array.shape
            all_gaps = []  # Collect all gap lengths for the current key

            for i in range(column):
                gaps = self.find_gaps(np_array[:, i])
                all_gaps.extend(gaps)  # Collect gaps across all columns for the current key

            # Filter out gaps below a certain threshold, e.g., 20
            filtered_gaps = [gap for gap in all_gaps if gap > 10]

            # If there are no gaps meeting the criteria for the current key
            if not filtered_gaps:
                gap_stats[key] = {"No significant gaps": None}
                continue  # Skip to the next key

            unique_gaps, gap_counts = np.unique(filtered_gaps, return_counts=True)
            
            # Calculate frequency statistics for the current key
            max_gap_length = unique_gaps[np.argmax(gap_counts)]
            mean_gap_length = np.mean(unique_gaps)
            median_gap_length = np.median(unique_gaps)

            max_gap_freq = np.max(gap_counts)
            mean_gap_freq = np.mean(gap_counts)
            median_gap_freq = np.median(gap_counts)

            # Store the statistics in the gap_stats dictionary under the current key
            gap_stats[key] = {
                "unique_gaps": unique_gaps,
                "gap_counts": gap_counts,
                "max_gap_length": max_gap_length,
                "mean_gap_length": mean_gap_length,
                "median_gap_length": median_gap_length,
                "max_gap_freq": max_gap_freq,
                "mean_gap_freq": mean_gap_freq,
                "median_gap_freq": median_gap_freq
            }

            self.__gaps.append(gap_stats[key]["max_gap_length"])

        return gap_stats, self.__gaps
    
    def find_block_lengths(self, col):
        """
        Purpose:
            Find the block lengths in a column of the numpy array and return block lengths
        Args:
            col: A column from a numpy array
        Returns:
            A list of block lengths in the column
        """
        indices = np.where(np.abs(col) == 1)[0]
        if len(indices) < 2:  # Need at least two '1's to form a block
            return [] if len(indices) == 0 else [1]  # No blocks or a single '1' considered as a block of length 1
        block_starts_ends = np.diff(indices, prepend=-1, append=col.size)  # Identify gaps to determine blocks
        block_lengths = block_starts_ends[block_starts_ends > 1]  # Block lengths are where the gap is more than 1
        return block_lengths.tolist()
    
    def min_blocklength_minblock(self, ploidy_data):
        """
        Purpose:
            Calculate minimum block length statistics for each key in the ploidy_data dictionary
        Args:
            None
        Returns:
            A list of minimum block lengths for each key in the ploidy_data dictionary
        """
        block_stats = {}  # Dictionary to store block length statistics for each key
        self.__blocklength = []

        for key, np_array in ploidy_data.items():
            all_blocks = []  # Collect all block lengths for the current key

            row, column = np_array.shape
            for i in range(column):
                blocks = self.find_block_lengths(np_array[:, i])
                all_blocks.extend(blocks)  # Collect blocks across all columns for the current key

            # Filter out block lengths below a certain threshold, e.g., 50
            filtered_blocks = [block for block in all_blocks if block > 20]

            if not filtered_blocks:
                block_stats[key] = {"No significant blocks": None}
                continue  # Skip to the next key

            unique_blocks, block_counts = np.unique(filtered_blocks, return_counts=True)
            
            # Calculate statistics for the current key
            min_block_length = unique_blocks[np.argmax(block_counts)]
            mean_block_length = np.mean(unique_blocks)
            median_block_length = np.median(unique_blocks)

            max_block_freq = np.max(block_counts)
            mean_block_freq = np.mean(block_counts)
            median_block_freq = np.median(block_counts)

            # Store the statistics in the block_stats dictionary under the current key
            block_stats[key] = {
                "unique_blocks": unique_blocks,
                "block_counts": block_counts,
                "min_block_length": min_block_length,
                "mean_block_length": mean_block_length,
                "median_block_length": median_block_length,
                "max_block_freq": max_block_freq,
                "mean_block_freq": mean_block_freq,
                "median_block_freq": median_block_freq
            }

            self.__blocklength.append(block_stats[key]["min_block_length"])

        return block_stats, self.__blocklength


    def min_blocklength(self, ploidy_data):
        """
        Purpose:
            Calculate minimum block length statistics for each key in the ploidy_data dictionary
        Args:
            None
        Returns:
            A list of minimum block lengths for each key in the ploidy_data dictionary
        """
        block_stats = {}  # Dictionary to store block length statistics for each key

        for key, np_array in ploidy_data.items():
            all_blocks = []  # Collect all block lengths for the current key

            row, column = np_array.shape
            for i in range(column):
                blocks = self.find_block_lengths(np_array[:, i])
                all_blocks.extend(blocks)  # Collect blocks across all columns for the current key

            # Filter out block lengths below a certain threshold, e.g., 50
            filtered_blocks = [block for block in all_blocks if block > 50]

            if not filtered_blocks:
                block_stats[key] = {"No significant blocks": None}
                continue  # Skip to the next key

            unique_blocks, block_counts = np.unique(filtered_blocks, return_counts=True)
            
            # Calculate statistics for the current key
            min_block_length = unique_blocks[np.argmax(block_counts)]
            mean_block_length = np.mean(unique_blocks)
            median_block_length = np.median(unique_blocks)

            max_block_freq = np.max(block_counts)
            mean_block_freq = np.mean(block_counts)
            median_block_freq = np.median(block_counts)

            # Store the statistics in the block_stats dictionary under the current key
            block_stats[key] = {
                "unique_blocks": unique_blocks,
                "block_counts": block_counts,
                "min_block_length": min_block_length,
                "mean_block_length": mean_block_length,
                "median_block_length": median_block_length,
                "max_block_freq": max_block_freq,
                "mean_block_freq": mean_block_freq,
                "median_block_freq": median_block_freq
            }

            self.__blocklength.append(block_stats[key]["min_block_length"])

        return block_stats, self.__blocklength

    def plot_thresholds(self, gap_stats, block_stats):
        """
        Purpose:
            Plot the gap and block length distributions across the ploidy data
        Args:
            gap_stats: The gap statistics
            block_stats: The block length statistics
        Returns:
            None
        """
        
        for key, value in gap_stats.items():
            unique_gaps = value["unique_gaps"]
            gap_counts = value["gap_counts"]

            max_gap_length = value["max_gap_length"]

            # Find the index of the gap length closest to the gap length with maximum frequency
            index_max_gap = np.argmax(gap_counts)

            # Calculate the start and end indices for the range
            start_index = max(0, index_max_gap - 10)
            end_index = min(len(unique_gaps), index_max_gap + 20)

            # Extend the x-axis range
            extended_start = max(0, start_index - 1)
            extended_end = min(len(unique_gaps), end_index + 1)

            # Extract the range of values around the maximum gap length
            plot_gaps = unique_gaps[extended_start:extended_end]
            plot_counts = gap_counts[extended_start:extended_end]

            # Plotting frequency of gaps within the extended range
            fig, ax = plt.subplots(figsize=(10, 6))

            # Connect the points to form a curve
            ax.plot(plot_gaps, plot_counts, marker='o', linestyle='-', linewidth=2)


            # Add a red dashed line indicating the peak value
            peak_gap = unique_gaps[index_max_gap]
            ax.axvline(peak_gap, color='red', linestyle='--', linewidth=3)

            # Add text annotation for the peak value
            peak_count = gap_counts[index_max_gap]
            ax.text(peak_gap, peak_count, f'Peak: {max_gap_length}', color='red', ha='center', va='bottom', fontsize=15)

            ax.set_xlim(plot_gaps.min(), plot_gaps.max())  # Set x-axis limits to the extended range
            ax.set_xlabel('Gap Length', fontsize=12)
            ax.set_ylabel('Frequency', fontsize=12)
            ax.set_title('Frequency of Gap Lengths', fontsize=16)
            ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

            plt.show()

            # Save plot as figure
            # fig.savefig(f'{key}_gap_frequency.png', dpi=300, bbox_inches='tight')
        
        for key, value in block_stats.items():
            unique_lengths = value["unique_blocks"]
            length_counts = value["block_counts"]

            # Calculate the maximum frequency block length
            max_length = unique_lengths[np.argmax(length_counts)]
            max_frequency = np.max(length_counts)

            # Plot frequency of filtered block lengths with smooth curves
            fig, ax = plt.subplots(figsize=(10, 6))

            # Find the index of the maximum length
            max_length_index = np.argmax(length_counts)

            # Get the values 10 before and after the maximum length index
            start_index = max(0, max_length_index - 10)
            end_index = min(len(unique_lengths), max_length_index + 20)

            # Extend the x-axis range
            extended_start = max(0, start_index - 1)
            extended_end = min(len(unique_lengths), end_index + 1)

            x = unique_lengths[extended_start:extended_end]
            y = length_counts[extended_start:extended_end]

            # Use spline interpolation for smooth curves
            x_smooth = np.linspace(x.min(), x.max(), 300)
            y_smooth = make_interp_spline(x, y)(x_smooth)

            ax.plot(x_smooth, y_smooth, color='blue', linewidth=2)
            ax.scatter(x, y, color='red', marker='o', label='Data Points')

            # Add a red dashed line indicating the peak value
            ax.axvline(max_length, color='red', linestyle='--')

            # Add text annotation for the peak value
            ax.text(max_length, max_frequency, f'Peak: {max_length}', color='red', ha='center', va='bottom')

            ax.set_xlim(x.min(), x.max())  # Set x-axis limits to the extended range
            ax.set_xlabel('Block Length', fontsize=12)
            ax.set_ylabel('Frequency', fontsize=12)
            ax.set_title('Frequency of Block Lengths', fontsize=16)
            ax.legend()
            ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

            plt.show()

            # Save plot as figure
            # fig.savefig(f'{key}_block_length_frequency.png', dpi=300, bbox_inches='tight')



                
