import pandas as pd

# Define the paths to the subgenome and BED files
subgenome_file_path = 'Final_result_B_Bju.xlsx'
bed_file_path = 'Bju.bed'
output_file_path = 'Bju_sub1.bed'

# Read the subgenome file using pandas
subgenome_df = pd.read_excel(subgenome_file_path)

# Extract gene IDs belonging to subgenome1 and ignore 'x'
subgenome1_genes = [gene for gene in subgenome_df['subgenome1'] if gene != 'x']

# Open the output file to write the results
with open(output_file_path, 'w') as outfile:
    # Read the BED file line by line
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            columns = line.strip().split('\t')
            gene_id = columns[3]  # Get the gene ID from the fourth column (index 3)
            if gene_id in subgenome1_genes:
                outfile.write(line)
                subgenome1_genes.remove(gene_id)  # Remove the found gene ID to speed up future searches

print(f'Genes belonging to subgenome1 have been extracted to {output_file_path}')
