import pandas as pd
import os
import re
import argparse

def extract_id_chromosome_start_end(df, feature_type, id_pattern, id_column=8):
    # Extract rows with the specified feature type
    df_feature = df[df[2] == feature_type].copy()
    # Extract the ID from the specified pattern
    df_feature['extracted_id'] = df_feature[id_column].str.extract(id_pattern, expand=False)

    # Group by extracted ID to get the start from the first CDS and end from the last CDS
    start_dict = df_feature.groupby('extracted_id')[3].min().to_dict()
    end_dict = df_feature.groupby('extracted_id')[4].max().to_dict()
    chromosome_dict = dict(zip(df_feature['extracted_id'], df_feature[0]))
    return start_dict, end_dict, chromosome_dict

def add_new_start_and_chromosome_to_blast(blast_file_path, query_file_path, subject_file_path):
    # Read the full BLAST file into a Pandas DataFrame
    df_blast = pd.read_csv(blast_file_path, sep='\t', header=None)

    # Read the query GFF file
    query_extension = os.path.splitext(query_file_path)[1]
    if query_extension == '.gff':
        df_query = pd.read_csv(query_file_path, sep='\t', header=None, comment='#')
        query_start_dict, query_end_dict, query_chromosome_dict = extract_id_chromosome_start_end(
            df_query, 'CDS', r'Parent=([^;]+)'
        )
    elif query_extension == '.bed':
        # Handle BED
        df_query = pd.read_csv(query_file_path, sep='\t', header=None)
        query_start_dict = dict(zip(df_query[3], df_query[1]))
        query_chromosome_dict = dict(zip(df_query[3], df_query[0]))
        query_end_dict = dict(zip(df_query[3], df_query[2]))

    # Read the subject GFF file
    subject_extension = os.path.splitext(subject_file_path)[1]
    if subject_extension == '.gff3':
        df_subject = pd.read_csv(subject_file_path, sep='\t', header=None, comment='#')
        subject_start_dict, subject_end_dict, subject_chromosome_dict = extract_id_chromosome_start_end(
            df_subject, 'CDS', r'Parent=([^;]+)'
        )
        # Extract the first attribute before the comma in the Parent field
        subject_start_dict = {k.split(',')[0]: v for k, v in subject_start_dict.items()}
        subject_end_dict = {k.split(',')[0]: v for k, v in subject_end_dict.items()}
        subject_chromosome_dict = {k.split(',')[0]: v for k, v in subject_chromosome_dict.items()}
    elif subject_extension == '.bed':
        # Handle BED
        df_subject = pd.read_csv(subject_file_path, sep='\t', header=None)
        subject_start_dict = dict(zip(df_subject[3], df_subject[1]))
        subject_chromosome_dict = dict(zip(df_subject[3], df_subject[0]))   
        subject_end_dict = dict(zip(df_subject[3], df_subject[2]))

    # Map the dictionaries to the BLAST DataFrame to create new columns
    df_blast['new_query_start'] = df_blast[0].map(query_start_dict)
    df_blast['new_query_end'] = df_blast[0].map(query_end_dict)
    df_blast['new_subject_start'] = df_blast[1].map(subject_start_dict)
    df_blast['new_subject_end'] = df_blast[1].map(subject_end_dict)
    df_blast['query_chromosome'] = df_blast[0].map(query_chromosome_dict)
    df_blast['subject_chromosome'] = df_blast[1].map(subject_chromosome_dict)

    # Update start and end columns for query and subject
    df_blast[6] = df_blast['new_query_start']
    df_blast[7] = df_blast['new_query_end']
    df_blast[8] = df_blast['new_subject_start']
    df_blast[9] = df_blast['new_subject_end']

    # Select and rearrange specific columns
    df_blast_final = df_blast[['query_chromosome', 0, 6, 7, 'subject_chromosome', 1, 8, 9, 10, 11]]

    # Rename columns
    df_blast_final.columns = [
        'chromosome_query', 'query_geneid', 'query_start', 'query_end',
        'chromosome_subject', 'subject_geneid', 'subject_start', 'subject_end',
        'evalue', 'bit_score'
    ]

    # Write the output file
    df_blast_final.to_csv('transformed_blast_output_with_selected_columns.blast', sep='\t', index=False, header=False)

    return df_blast_final

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BLAST output with GFF3 or BED files.')
    parser.add_argument('blast_file', type=str, help='Path to the BLAST output file')
    parser.add_argument('query_file', type=str, help='Path to the query GFF3 or BED file')
    parser.add_argument('subject_file', type=str, help='Path to the subject GFF3 or BED file')
    parser.add_argument('--output_file', type=str, default='transformed_blast_output_with_selected_columns.blast', help='Path to the output file')

    args = parser.parse_args()

    # Run the function
    output_df = add_new_start_and_chromosome_to_blast(args.blast_file, args.query_file, args.subject_file)

    # Save the output DataFrame to the specified output file
    output_df.to_csv(args.output_file, sep='\t', index=False, header=False)