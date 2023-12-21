import re

def extract_gene_ids(file_path):
    gene_ids = set()
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) > 2:
                    gene_ids.add(parts[1].strip())  # Extract gene ID 1
                    gene_ids.add(parts[2].strip())  # Extract gene ID 2
    return gene_ids

def filter_gtf_file(gtf_file_path, gene_ids, output_file_path):
    with open(gtf_file_path, 'r') as gtf_file, open(output_file_path, 'w') as output_file:
        for line in gtf_file:
            if any(gene_id in line for gene_id in gene_ids):
                output_file.write(line)

# File paths
collinear_file_path = 'modified_collinear_file_sub2.colinearity'
gtf_file_path = 'ce_ba.gff'
output_gtf_file_path = 'filtered_gtf_sub2.gff'

# Extract gene IDs
extracted_gene_ids = extract_gene_ids(collinear_file_path)

# Filter the GTF file
filter_gtf_file(gtf_file_path, extracted_gene_ids, output_gtf_file_path)

print(f"Filtered GTF file has been created and saved to {output_gtf_file_path}")
