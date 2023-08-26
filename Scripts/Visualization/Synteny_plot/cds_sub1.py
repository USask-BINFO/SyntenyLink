import pandas as pd

# Define the paths to the subgenome and CDS files
subgenome_file_path = 'Final_result_A_Bju.xlsx'
cds_file_path = 'Bjuncea_3DH.cds_20211001.fasta'
output_file_path = 'Bju_sub2.cds'

# Read the subgenome file using pandas
subgenome_df = pd.read_excel(subgenome_file_path)

# Create a set of gene IDs from subgenome1 column
subgenome1_genes = set(subgenome_df['subgenome2'].loc[subgenome_df['subgenome2'] != 'x'])

# Open the output file to write the results
with open(output_file_path, 'w') as outfile:
    # Read the CDS file line by line
    with open(cds_file_path, 'r') as cds_file:
        write_sequence = False
        for line in cds_file:
            if line.startswith('>'):
                gene_id = line[1:].strip().split()[0]  # Split by whitespace and get the first part
                write_sequence = gene_id in subgenome1_genes
            if write_sequence:
                outfile.write(line)

print(f'Sequences belonging to subgenome1 have been extracted to {output_file_path}')
