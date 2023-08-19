# Importing necessary libraries
import pandas as pd

# Loading the Excel file
file_path = "Final_result_bni.xlsx"
data = pd.read_excel(file_path)

# Displaying the first few rows to understand the structure of the data
data.head()
# Loading the BED files for Arabidopsis and Brassica genomes
arabidopsis_bed_file = "TAIR10.bed"
# brassica_bed_file = "/mnt/data/Brapa_genome_v3.0_genes.bed"

# Reading the BED files
arabidopsis_bed = pd.read_csv(arabidopsis_bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_id'])
# brassica_bed = pd.read_csv(brassica_bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_id'])

# Displaying the first few rows of both BED files
# arabidopsis_bed.head(), brassica_bed.head()
# Reloading the Arabidopsis BED file with the correct format
arabidopsis_bed = pd.read_csv(arabidopsis_bed_file, sep='\t', header=None,
                              names=['chromosome', 'start', 'end', 'gene_id', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])

# Displaying the first few rows to confirm the correct format
arabidopsis_bed.head()

# Matching the Arabidopsis gene IDs from the Excel file with the corresponding entries in the Arabidopsis BED file
arabidopsis_coordinates = arabidopsis_bed[arabidopsis_bed['gene_id'].isin(data['AT_Geneid'])]

# Merging the Arabidopsis coordinates with the original Excel file based on the gene IDs
combined_data_arabidopsis = pd.merge(data, arabidopsis_coordinates[['gene_id', 'chromosome', 'start', 'end']], left_on='AT_Geneid', right_on='gene_id', how='left')
combined_data_arabidopsis.drop(columns=['gene_id'], inplace=True)  # Dropping the duplicate gene_id column

# Displaying the first few rows of the combined data
combined_data_arabidopsis.head()

# Defining the path for the new file
output_file_path = "Final_result_with_Arabidopsis_coordinates_bni.xlsx"

# Saving the combined data with Arabidopsis coordinates to the new file
combined_data_arabidopsis.to_excel(output_file_path, index=False)

# Providing the link to download the file
output_file_path
