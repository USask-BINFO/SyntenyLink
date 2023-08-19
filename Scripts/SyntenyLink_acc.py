
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import math
import random

def subgenome_overlap(GT_file_path, df_subgenome_density_file_path, df_synteny, num_subgenomes, first_letter_get):
    # Read the data
    df_subgenome_density = pd.read_excel(df_subgenome_density_file_path)
    # Create a list of subgenomes
    subgenomes = [[] for _ in range(num_subgenomes)]

    for i in range(len(df_subgenome_density["Row start #"])):
        for j in range(df_subgenome_density["Row start #"].values[i], df_subgenome_density["Row end #"].values[i] + 1):
            for k in range(num_subgenomes):
                subgenomes[k].append(df_synteny[df_subgenome_density[f"subgenome{k+1}"].values[i]][j])

    # Create a dataframe of subgenomes
    df_subgenome = pd.DataFrame(subgenomes).transpose()

    # Change the column names
    df_subgenome.columns = [f"subgenome{i+1}" for i in range(num_subgenomes)]
    # print(df_subgenome)
    #save df_subgenome as an excel file
    df_subgenome.to_excel("final_result.xlsx")

    # Load the dataframes
    df_groundtruth = pd.read_excel(GT_file_path) # dataframe for groundtruth

    for j in range(num_subgenomes):
        # Check for exact matches of gene ids and 'x' between subgenome1 of both dataframes
        exact_matches = 0
        total_genes = 0
        missing_genes = 0
        missing_genes_list = []

        for i, row_groundtruth in df_groundtruth.iterrows():
            if i < len(df_groundtruth) - 1:
                    gene_id_groundtruth = row_groundtruth[f"{first_letter_get}_subgenome{j+1}"]
                    if gene_id_groundtruth != 'x':
                        total_genes += 1
                        if i < len(df_subgenome):
                            row_subgenome = df_subgenome.iloc[i]
                            gene_id_subgenome = row_subgenome[f"subgenome{j+1}"]
                            if gene_id_groundtruth == gene_id_subgenome:
                                exact_matches += 1
                            else:
                                missing_genes += 1
                                missing_genes_list.append(gene_id_groundtruth)

        # Calculate the overlap percentage for each subgenome
        overlaps = exact_matches / total_genes

        print(f"Exact match percentage for subgenome{j+1}: {overlaps:.2%}")
        print(f"Exact match number for subgenome{j+1}: {exact_matches}")
        print(f"Missing genes for subgenome{j+1}: {missing_genes}")

def subgenome_overlap_separate(GT_file_path, df_subgenome_density_file_path, df_synteny, num_subgenomes, s, first_letter_get):
    # Read the data
    df_subgenome_density = pd.read_excel(df_subgenome_density_file_path)

    # Create a list of subgenomes
    subgenomes = [[] for _ in range(num_subgenomes)]

    for i in range(len(df_subgenome_density["Row start #"])):
        for j in range(df_subgenome_density["Row start #"].values[i], df_subgenome_density["Row end #"].values[i] + 1):
            for k in range(num_subgenomes):
                subgenomes[k].append(df_synteny[df_subgenome_density[f"subgenome{k+1}"].values[i]][j])

    # Create a dataframe of subgenomes
    df_subgenome = pd.DataFrame(subgenomes).transpose()

    # Change the column names
    df_subgenome.columns = [f"subgenome{i+1}" for i in range(num_subgenomes)]
    # print(df_subgenome)

    # Load the dataframes
    df_groundtruth = pd.read_excel(GT_file_path) # dataframe for groundtruth

    # Check for exact matches of gene ids and 'x' between subgenome1 of both dataframes
    exact_matches = 0
    total_genes = 0
    missing_genes = 0
    missing_genes_list = []

    for i, row_groundtruth in df_groundtruth.iterrows():
        if i < len(df_groundtruth) - 1:
                gene_id_groundtruth = row_groundtruth[f"{first_letter_get}_subgenome{s+1}"]
                if gene_id_groundtruth != 'x':
                    total_genes += 1
                    if i < len(df_subgenome):
                        row_subgenome = df_subgenome.iloc[i]
                        gene_id_subgenome = row_subgenome[f"subgenome{s+1}"]
                        if gene_id_groundtruth == gene_id_subgenome:
                            exact_matches += 1
                        else:
                            missing_genes += 1
                            missing_genes_list.append(gene_id_groundtruth)

    # Calculate the overlap percentage for each subgenome
    overlaps = exact_matches / total_genes

    print(f"Exact match percentage for subgenome{s+1}: {overlaps:.2%}")
    print(f"Exact match number for subgenome{s+1}: {exact_matches}")
    print(f"Missing genes for subgenome{s+1}: {missing_genes}")