import re

# Function to remove the last decimal from gene IDs if present
def remove_version(gene_id):
    return re.sub(r'\.\d+$', '', gene_id)

# Read the gene combination file and store the desired combinations
desired_combinations = set()
with open('sub1_list.txt', 'r') as gene_combination_file:
    for line in gene_combination_file:
        if not line.startswith('locus_id'):
            columns = line.strip().split('\t')
            if len(columns) == 2:
                gene1, gene2 = columns
                gene1 = remove_version(gene1)
                gene2 = remove_version(gene2)
                desired_combinations.add((gene1, gene2))

# Read the collinear file and retain only the synteny blocks with desired gene combinations
output_lines = []
block_lines = []
alignment_header = ''
with open('ce_ba.collinearity', 'r') as collinear_file:
    for line in collinear_file:
        if line.startswith('## Alignment'):
            if block_lines:
                filtered_block_lines = [
                    line for line in block_lines if
                    (remove_version(line.split('\t')[1]), remove_version(line.split('\t')[2])) in desired_combinations
                ]
                if filtered_block_lines:
                    output_lines.append(alignment_header)
                    output_lines.extend(filtered_block_lines)
            block_lines = []
            alignment_header = line.strip()
        elif not line.startswith('#'):
            block_lines.append(line.strip())

# Write the modified collinear file
with open('modified_collinear_file_sub1.colinearity', 'w') as output_file:
    for line in output_lines:
        output_file.write(line + '\n')

print("Modified collinear file created successfully.")
