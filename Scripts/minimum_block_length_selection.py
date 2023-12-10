import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
import sys

# Function for finding gaps in each column
def find_gaps(col, gap_threshold, m):
    """
    Finds the gaps in each column of the matrix
    col: column of the matrix
    gap_threshold: threshold for the gap
    m: number of rows in the matrix

    returns: lengths of the blocks
    """
    # Indices for ones in the given column
    indices_one = np.where(col == 1)[0]
    # Calculate difference between consecutive indices
    d = np.diff(indices_one)

    # If the difference is bigger than the gap threshold, store the lengths of the blocks
    initial_gaps = np.where(d > gap_threshold)[0]
    block_lengths = indices_one[initial_gaps + 1] - indices_one[initial_gaps]

    return block_lengths

# Get the input file as an argument (collinear file)
input_file = sys.argv[sys.argv.index('-i') + 1]

# Convert the collinear file to a dataframe
C_df_csv = pd.read_csv(input_file, sep='\t')
# Make a copy of the dataframe
C_df = C_df_csv.copy()

C_df_updated = C_df.iloc[1:, 1:-3]
# Convert all the 'x' to 0 and all the other entries which are not equal to 'x' to 1 omitting first column and last three columns and first row
C_df_updated = C_df_updated.replace(r'^x$', 0, regex=True)
# Convert all the other entries that start with a letter to 1
C_df_updated = C_df_updated.replace(r'^[a-zA-Z]', 1, regex=True)
# Set first row index starts from 0
C_df_updated.index = C_df_updated.index - 1

def split_df(C_df_updated, gap_thresholds, k= 0):
    # Split the dataframe based on the first letter of column names
    dfs = {}
    for column in C_df_updated.columns:
        print(column)
        first_letter = column[0]
        if first_letter not in dfs:
            dfs[first_letter] = pd.DataFrame()
        dfs[first_letter][column] = C_df_updated[column]
    
    # Print the resulting dataframes
    for key, value in dfs.items():
        first_letter_get = key
        print(f"Dataframe with columns starting with '{key}':")
        print(value)
        print()
        # Convert the dataframe to a numpy array
        C = value.to_numpy()
        m, n = C.shape
        print(m, n)

        # Gap_threshold is the threshold to determine if there is a gap or not
        # Min_block_length is the minimum length of the signal blocks

        gap_threshold = gap_thresholds[k]
        print(f"Gap threshold for dataframe {k+1}: {gap_threshold}")
        k += 1

        #21

        # Get all gaps in all columns
        block_lengths = [find_gaps(C[:, i], gap_threshold, m) for i in range(n)]
        print(block_lengths)

        # Filter block lengths greater than 0
        filtered_lengths = np.concatenate(block_lengths)
        filtered_lengths = filtered_lengths[filtered_lengths > 50]

        # Calculate frequency of filtered block lengths
        unique_lengths, length_counts = np.unique(filtered_lengths, return_counts=True)

        # Find the maximum, mean, and median block lengths and their frequencies
        max_length = unique_lengths[np.argmax(length_counts)]
        mean_length = np.mean(unique_lengths)
        median_length = np.median(unique_lengths)

        max_frequency = np.max(length_counts)
        mean_frequency = np.mean(length_counts)
        median_frequency = np.median(length_counts)

        print("Block Length for Maximum Frequency:", max_length)
        print("Block Length for Mean Frequency:", mean_length)
        print("Block Length for Median Frequency:", median_length)
        print("Max Frequency:", max_frequency)
        print("Mean Frequency:", mean_frequency)
        print("Median Frequency:", median_frequency)

        # New condition: Consider block lengths only greater than 100 if mean frequency is less than 600
        if median_length < 600:
            filtered_lengths = filtered_lengths[filtered_lengths > 100]
            unique_lengths, length_counts = np.unique(filtered_lengths, return_counts=True)
            # Recalculate maximum, mean, and median block lengths and their frequencies
            max_length = unique_lengths[np.argmax(length_counts)]
            mean_length = np.mean(unique_lengths)
            median_length = np.median(unique_lengths)

            max_frequency = np.max(length_counts)
            mean_frequency = np.mean(length_counts)
            median_frequency = np.median(length_counts)

            print("Updated Block Length for Maximum Frequency:", max_length)
            print("Updated Block Length for Mean Frequency:", mean_length)
            print("Updated Block Length for Median Frequency:", median_length)
            print("Updated Max Frequency:", max_frequency)
            print("Updated Mean Frequency:", mean_frequency)
            print("Updated Median Frequency:", median_frequency)

        # Calculate the maximum frequency block length
        max_length = unique_lengths[np.argmax(length_counts)]
        max_frequency = np.max(length_counts)

        # Plot frequency of filtered block lengths with smooth curves
        fig, ax = plt.subplots(figsize=(10, 6))

        # Find the index of the maximum length
        max_length_index = np.argmax(length_counts)

        # Get the values 10 before and after the maximum length index
        start_index = max(0, max_length_index - 10)
        end_index = min(len(unique_lengths), max_length_index + 11)

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
        fig.savefig(f'{first_letter_get}_block_length_frequency.png', dpi=300, bbox_inches='tight')

# Get the gap thresholds as command line arguments
gap_thresholds = []
if '-g' in sys.argv:
    index = sys.argv.index('-g')
    for i in range(index+1, len(sys.argv)):
        if sys.argv[i].startswith('-'):
            break
        gap_thresholds.append(float(sys.argv[i]))

split_df(C_df_updated, gap_thresholds)