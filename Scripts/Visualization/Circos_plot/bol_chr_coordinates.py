import pandas as pd
# Loading the Brassica BED file
brassica_bed_path = "Boleracea.v2.1.genes.bed"
brassica_bed_data = pd.read_csv(brassica_bed_path, sep='\t', header=None, names=['chromosome', 'gene_id', 'start', 'end'])

# Displaying the first few rows of the Brassica BED file to understand its structure
brassica_bed_data.head()

# # Splitting the combined string into separate columns
# brassica_bed_data[['chromosome', 'gene_id', 'start', 'end']] = brassica_bed_data['chromosome']

# # Converting the 'start' and 'end' columns to integers
# brassica_bed_data['start'] = brassica_bed_data['start'].astype(int)
# brassica_bed_data['end'] = brassica_bed_data['end'].astype(int)

# Calculating the start and end coordinates for each Brassica chromosome
brassica_chromosomes = brassica_bed_data.groupby('chromosome').agg(
    start=('start', 'min'),
    end=('end', 'max')
).reset_index()

# Saving the Brassica chromosomes data to an Excel file
brassica_chromosomes_excel_path = "brassica_chromosomes_coordinates.xlsx"
brassica_chromosomes.to_excel(brassica_chromosomes_excel_path, index=False)

# Displaying the first few rows of the Brassica chromosomes data and providing the download link
brassica_chromosomes.head(), brassica_chromosomes_excel_path
