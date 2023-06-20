import numpy as np
import pandas as pd
import re
import warnings
import sys
import pickle
import csv
import matplotlib.pyplot as plt
from scipy import stats

# Function for finding gaps in each column
def find_gaps(col, m):
    """
    Finds the gaps in each column of the matrix
    col: column of the matrix
    m: number of rows in the matrix

    returns: indices of the beginning and end of the gaps
    """
    gaps = []
    indices_one = np.where(col == 1)[0]
    # print(indices_one)
    # calculate difference between consecutive indices
    d = np.diff(indices_one)
    print(d)
    gaps.extend(d)

    return gaps

# Get the input file as an argument (collinear file)
input_file = sys.argv[1]

# Convert the collinear file to a dataframe
C_df_csv = pd.read_csv(input_file, sep=' ', header=None)
# Make a copy of the dataframe
C_df = C_df_csv.copy()
C_df_updated = C_df.iloc[1:, 1:-3]

# Convert all the 'x' to 0 and all the other entries which are not equal to 'x' to 1
C_df_updated = C_df_updated.replace(r'^x$', 0, regex=True)
# Convert all the other entries starting with a letter to 1
C_df_updated = C_df_updated.replace(r'^[a-zA-Z]', 1, regex=True)

# Set first row index starts from 0
C_df_updated.index = C_df_updated.index - 1

# Convert the dataframe to a numpy array
C = C_df_updated.to_numpy()
m, n = C.shape
print(m, n)

gaps_cell = [find_gaps(C[:, i], m) for i in range(n)]

# Calculate frequency of gaps inside each block
freqs = []

for gaps in gaps_cell:
    # Ignore gaps equal to zero
    gaps = [gap for gap in gaps if gap > 20]
    freqs.extend(gaps)

unique_gaps, gap_counts = np.unique(freqs, return_counts=True)

# Find the gap length for maximum, mean, and median frequency
max_gap_length = unique_gaps[np.argmax(gap_counts)]
mean_gap_length = np.mean(unique_gaps)
median_gap_length = np.median(unique_gaps)

max_gap_freq = np.max(gap_counts)
mean_gap_freq = np.mean(gap_counts)
median_gap_freq = np.median(gap_counts)

print("Gap Length for Maximum Frequency:", max_gap_length)
print("Gap Length for Mean Frequency:", mean_gap_length)
print("Gap Length for Median Frequency:", median_gap_length)
print("Max Frequency:", max_gap_freq)
print("Mean Frequency:", mean_gap_freq)
print("Median Frequency:", median_gap_freq)

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
fig.savefig('gap_frequency.png', dpi=300, bbox_inches='tight')
