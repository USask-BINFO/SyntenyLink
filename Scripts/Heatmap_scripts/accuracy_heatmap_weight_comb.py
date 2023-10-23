import pandas as pd
import re
import warnings
import sys
import pickle
import csv
import os
import seaborn as sns
import matplotlib.pyplot as plt


def subgenome_overlap(GT_file_path, df_subgenome_density_file_path, df_synteny, num_subgenomes, sub_no, first_letter_get):
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
    df_subgenome.columns = [f"{first_letter_get}_subgenome{i+1}" for i in range(num_subgenomes)]
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
                gene_id_groundtruth = row_groundtruth[f"{first_letter_get}_subgenome{sub_no+1}"]
                if gene_id_groundtruth != "x":
                    total_genes += 1
                    if i < len(df_subgenome):
                        row_subgenome = df_subgenome.iloc[i]
                        gene_id_subgenome = row_subgenome[f"{first_letter_get}_subgenome{sub_no+1}"]
                        if gene_id_groundtruth == gene_id_subgenome:
                            exact_matches += 1
                        else:
                            missing_genes += 1
                            missing_genes_list.append(gene_id_groundtruth)

    # Calculate the overlap percentage for each subgenome
    overlaps = exact_matches / total_genes
    print(overlaps, sub_no)

    return overlaps

def generate_heatmap(acc_matrix, subgenome_index, w1, w2):
    sub_acc_matrix = acc_matrix[subgenome_index]
    sns.set(font_scale=0.8)
    plt.figure(figsize=(8, 10))
    ax = sns.heatmap(sub_acc_matrix, cmap='YlGnBu', annot=True, fmt=".2f", cbar=False,
                     xticklabels=[str(round(i*0.05,2)) for i in range(1, 21)], 
                     yticklabels=[str(round(i*0.05,2)) for i in range(1, 21)])
    plt.title("Subgenome " + str(subgenome_index + 1) + " Accuracy Matrix")
    plt.xlabel("Weight2")
    plt.ylabel("Weight1")
    plt.savefig("subgenome_" + str(subgenome_index + 1) + "weight_heatmap.png", dpi=300, bbox_inches='tight')