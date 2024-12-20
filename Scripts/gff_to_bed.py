import pandas as pd
import re

# Function to extract ID from the attributes column
def extract_id(attributes):
    match = re.search(r'ID=([^;]+)', attributes)
    return match.group(1) if match else None

# Read the input file
input_file = '/binfo-nas1/data/genomes/Descurainia/sophia/flixweed.gff'  # Replace with the actual file path
df = pd.read_csv(input_file, sep='\t', header=None, names=['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])

# Filter for mRNA rows
mrna_df = df[df['Type'] == 'mRNA']

# Extract necessary columns and process the Attributes column to get the ID
mrna_df['ID'] = mrna_df['Attributes'].apply(extract_id)

# Select and rename the relevant columns
result_df = mrna_df[['Chromosome', 'ID', 'Start', 'End']]

# Display the result
print(result_df)

# Optionally, save the result to a new file
output_file = 'Sophia/flixweed.bed'  # Replace with the desired output file path
result_df.to_csv(output_file, sep='\t', index=False)
