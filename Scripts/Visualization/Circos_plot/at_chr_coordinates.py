import pandas as pd
# Loading the Arabidopsis BED file
arabidopsis_bed_path = "TAIR10.bed"
arabidopsis_bed_data = pd.read_csv(arabidopsis_bed_path, sep='\t', header=None,
                              names=['chromosome', 'start', 'end', 'gene_id', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])

# Saving the Arabidopsis chromosomes data to an Excel file
arabidopsis_chromosomes_excel_path = "arabidopsis_chromosomes_coordinates.xlsx"

# # Splitting the combined string into separate columns
# arabidopsis_bed_data[['chromosome', 'start', 'end', 'gene_id']] = arabidopsis_bed_data['chromosome'].str.split('\t', expand=True)

# # Converting the 'start' and 'end' columns to integers
# arabidopsis_bed_data['start'] = arabidopsis_bed_data['start'].astype(int)
# arabidopsis_bed_data['end'] = arabidopsis_bed_data['end'].astype(int)

# Calculating the start and end coordinates for each Arabidopsis chromosome
arabidopsis_chromosomes = arabidopsis_bed_data.groupby('chromosome').agg(
    start=('start', 'min'),
    end=('end', 'max')
).reset_index()

# Saving the Arabidopsis chromosomes data to an Excel file
arabidopsis_chromosomes.to_excel(arabidopsis_chromosomes_excel_path, index=False)

# Displaying the first few rows of the Arabidopsis chromosomes data and providing the download link
arabidopsis_chromosomes.head(), arabidopsis_chromosomes_excel_path

