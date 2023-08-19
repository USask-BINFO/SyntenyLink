# Importing necessary libraries
import pandas as pd

# Reading the Arabidopsis Data with Subgenomes (Update the path accordingly)
final_result_path = "Final_result_bni.xlsx"
data = pd.read_excel(final_result_path)

# Reading the Brassica BED File
brassica_bed_path = "Bnigra_NI100.v2.genes.bed"
brassica_bed_data = pd.read_csv(brassica_bed_path, sep='\t', header=None, names=['chromosome','gene_id', 'start', 'end'])

# Splitting and Formatting the Brassica BED Data
# brassica_bed_data[['chromosome','gene_id', 'start', 'end']] = brassica_bed_data['chromosome'].str.split('\t', expand=True) # Assuming tab-separated values
# brassica_bed_data['start'] = brassica_bed_data['start'].astype(int)
# brassica_bed_data['end'] = brassica_bed_data['end'].astype(int)

# Merging Arabidopsis Data with Brassica Coordinates for Subgenome1
data_subgenome1 = pd.merge(data, brassica_bed_data, left_on='subgenome1', right_on='gene_id', how='left')

# The resulting DataFrame (e.g., data_subgenome1) contains the merged data for Subgenome1
data = data_subgenome1


# Function to merge Brassica coordinates with subgenome data and save to Excel file
def merge_and_save_subgenome(subgenome_column, brassica_coordinates, output_file_path):
    subgenome_data = data[['AT_Geneid', subgenome_column]].copy()
    subgenome_data = pd.merge(subgenome_data, brassica_coordinates, left_on=subgenome_column, right_on='gene_id', how='left')
    subgenome_data.drop(columns=['gene_id'], inplace=True)
    subgenome_data.to_excel(output_file_path, index=False)

    
# Paths for the Excel files for subgenomes
output_file_subgenome1_excel_path = "Final_result_with_Brassica_subgenome1.xlsx"
output_file_subgenome2_excel_path = "Final_result_with_Brassica_subgenome2.xlsx"
output_file_subgenome3_excel_path = "Final_result_with_Brassica_subgenome3.xlsx"

# Merging and saving the Excel files for each subgenome
merge_and_save_subgenome('subgenome1', brassica_bed_data, output_file_subgenome1_excel_path)
merge_and_save_subgenome('subgenome2', brassica_bed_data, output_file_subgenome2_excel_path)
merge_and_save_subgenome('subgenome3', brassica_bed_data, output_file_subgenome3_excel_path)
