#!/bin/bash
# This script is used to generate synteny links between two genomes
# Usage: ./SyntenyLink.sh ref_pep.fa query_cds.fa query_pep.fasta query_cds.fasta query.gff3 subject.gff3 ref_gene_list.txt ground_truth.txt bed_file.bed chr1 chr2

# # Check if the correct number of arguments are provided
# if [ "$#" -ne 11 ]; then
#     echo "Error: You must enter exactly 11 arguments"
#     echo "Usage: ./SyntenyLink.sh data/Bna/TAIR10_pep_20101214_updated.fa data/Bna/TAIR10_cds_20101214_updated.fa data/Bna/Bnapus_3DH.pep_20211001.fasta data/Bna/Bnapus_3DH.cds_20211001.fasta query.gff3 subject.gff3 data/ath_genes.txt groundtruth_bna.xlsx data/Bnapus_3DH.mrna_20211001.bed 10 9"
#     exit 1
# fi

# # Clone the repository
# # git clone https://git.cs.usask.ca/qnm481/syntenylink.git || { echo "Error: Failed to clone the repository"; exit 1; }
# cd syntenylink || { echo "Error: Failed to change directory to syntenylink"; exit 1; }

# # Install the required packages
# pip install -r requirements.txt || { echo "Error: Failed to install required packages"; exit 1; }

# # conda init
# # # Create a new environment
# conda create --name blast_env python=3.8

# # # Activate the environment
# conda activate blast_env


# # Install legacy blast (check if it's already installed)
# if ! command -v blastall &> /dev/null; then
#     conda install -c bioconda blast-legacy || { echo "Error: Failed to install blast-legacy"; exit 1; }
# fi

# Create the blast databases
# makeblastdb -in "$1" -dbtype prot -out ref_pep || { echo "Error: Failed to create reference protein BLAST database"; exit 1; }
# makeblastdb -in "$2" -dbtype nucl -out query_cds || { echo "Error: Failed to create query nucleotide BLAST database"; exit 1; }

# # Run the blast (increase -a value to use more CPUs if available)
# blastall -i "$3" -p blastp -d ref_pep -m 8 -e 1e-5 -F F -v 5 -b 5 -o abc.blast -a 4 || { echo "Error: Failed to run protein BLAST"; exit 1; }
# blastall -i "$4" -p blastn -d query_cds -m 8 -e 1e-5 -F F -v 5 -b 5 -o abc.blastn.blast -a 4 || { echo "Error: Failed to run nucleotide BLAST"; exit 1; }

# Filter the blast results
# perl SyntenyLink_bf.pl abc.blast > abc_blast_filtered_modified.txt || { echo "Error: Failed to filter protein BLAST results"; exit 1; }
# perl SyntenyLink_bf.pl abc.blastn.blast > abc_blastn_filtered_modified.txt || { echo "Error: Failed to filter nucleotide BLAST results"; exit 1; }

# Rearrange the blast results matching the input format of DAGchainer
# python3 transform_blast_to_dagchaine.py abc_blast_filtered_modified.txt "$5" "$6" || { echo "Error: Failed to transform protein BLAST results to DAGchainer format"; exit 1; }
# python3 transform_blast_to_dagchaine.py abc_blastn_filtered_modified.txt "$5" "$6" --output_file transformed_blast_output_with_selected_columns.blastn.blast || { echo "Error: Failed to transform nucleotide BLAST results to DAGchainer format"; exit 1; } 

# Run DAGchainer (ensure DAGchainer is installed before running)
# perl run_DAG_chainer.pl -i transformed_blast_output_with_selected_columns.blast -s -I || { echo "Error: Failed to run DAGchainer for protein"; exit 1; }
# perl run_DAG_chainer.pl -i transformed_blast_output_with_selected_columns.blastn.blast -s -I || { echo "Error: Failed to run DAGchainer for nucleotide"; exit 1; }

# Generate the colinear file (provide gene list file as input)
perl SyntenyLink_st.pl -d transformed_blast_output_with_selected_columns.blast.aligncoords -g "$7" || { echo "Error: Failed to generate colinear file"; exit 1; }

# Run SyntenyLink
python3 main.py -i transformed_blast_output_with_selected_columns.blast.aligncoords.success.colinear -s 2 -n 3 -p g -bed "$8" -chr1 "${9}" -sub1 chr -dag transformed_blast_output_with_selected_columns.blastn.blast.aligncoords || { echo "Error: Failed to run SyntenyLink"; exit 1; }
