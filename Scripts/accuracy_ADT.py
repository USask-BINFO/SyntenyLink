import Translocation as T
import pandas as pd
from itertools import permutations

class Acc(object):
    def __init__(self):
        """
        Purpose: Initialize the Main Breakpoint class object
        Args: None
        Returns: None
        """
        self.synteny = T.Translocation()
    
    def main(self):
        """
        Purpose: Run the main function
        Args: None
        Returns: None
        """
        
        # Run the main function
        formatted_dfs, header, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns, GT, n_sub, original_ploidy_numpy, last_column, ploidy_numpy_overlap, subgenome_letters = self.synteny.main()

        # # Convert gene ids to chromosome names they belong to
        chr = self.synteny.get_subgenome_chr(formatted_dfs, header, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns)
        x_updated = self.synteny.update_chr_x_with_previous_or_next_chr(chr)
        removed_errors = self.synteny.remove_inconsistencies_in_subgenome_assignments(x_updated)
        all_x_updated = self.synteny.update_chr_x_with_previous_or_next_chr(removed_errors)
        block_chr = self.synteny.chr_block_subgenome_assignment(all_x_updated)
        sub_pro = self.synteny.remove_overlapping_assignments(block_chr)
        sub_almost = self.synteny.swap_x_with_values_in_dict(sub_pro)
        adjusted_sub = self.synteny.adjust_subgenome_for_continuity(sub_almost, 20)
        sub_processed = self.synteny.adjust_all_subgenome_placements(adjusted_sub)
        sub_processed_final = self.synteny.merge_subgenome_blocks(sub_processed)
        final_sub = self.synteny.resolve_subgenome_conflicts_driver(sub_processed_final)
        final_sub = self.synteny.resolve_subgenome_conflicts_driver(sub_processed_final)
        sub_no_dup = self.synteny.remove_rows_with_duplicate_subgenomes_and_update_end(final_sub)

        # Check for duplicates in subgenomes
        self.synteny.check_duplicates_in_subgenomes(sub_no_dup)

        # Get the densities of the final subgenome placements
        densities_final = self.synteny.get_densities_subgenome_placements(sub_no_dup, original_ploidy_numpy, chrsets,ploidy_numpy_overlap)

        # Swap the subgenome assignments based on density
        swapped_sub, _ = self.synteny.swap_max_density_to_subgenome1_for_all(sub_no_dup, densities_final, n_sub)

        # Get the final subgenome placements
        subgenomes_final = self.synteny.get_subgenomes(swapped_sub, n_sub, chrsets, gene_id_matrix_df_, header, last_two_columns, last_column)

        return subgenomes_final, GT, n_sub, swapped_sub, subgenome_letters
    

    def accuracy(self, GT_file_path, n_sub, df_subgenome_dict, new_chr_placement, key_comb, subgenome_letters):
        """
        Purpose: Calculate accuracy of subgenome assignments compared to the ground truth.
        Args:
            GT_file_path: Path to ground truth file
            n_sub: Number of subgenomes
            df_subgenome_dict: Dictionary of subgenome DataFrames
            new_chr_placement: Dictionary with chromosome placement information
            key_comb: Key for combination results
            subgenome_letters: Dictionary of subgenome letter mappings
        Returns:
            Dictionary of accuracy results for each subgenome.
        """
        # Read the ground truth DataFrame from the provided file path
        df_groundtruth = pd.read_excel(GT_file_path)

        # Initialize a dictionary to store results for each subgenome
        results = {}

        # Iterate over each key in the chromosome placement dictionary
        for key, df_subgenome in new_chr_placement.items():
            # Extract the first non-'x' value in 'subgenome_1' and determine its letter
            first_chr_name = next((val for val in new_chr_placement[key]['subgenome_1'] if val != 'x'), 'x')
            first_letter = first_chr_name[0] if first_chr_name != 'x' and isinstance(first_chr_name, str) else 'x'

            # Set up accuracy results for this subgenome DataFrame
            results[key_comb] = []

            for j in range(n_sub):
                # Initialize metrics for accuracy calculation
                exact_matches = 0
                total_genes = 0
                missing_genes = 0

                # Retrieve subgenome-specific letter from subgenome_letters dictionary
                subgenome_letter = subgenome_letters.get(f'_subgenome_{j+1}', first_letter)

                # Compare gene placements between ground truth and subgenome DataFrame
                for i, row_groundtruth in df_groundtruth.iterrows():
                    if i < len(df_subgenome_dict):  # Ensure index is within range of df_subgenome
                        gene_id_groundtruth = row_groundtruth.get(f"{subgenome_letter}_subgenome{j+1}", 'x')
                        if gene_id_groundtruth != 'x':
                            total_genes += 1
                            row_subgenome = df_subgenome_dict.iloc[i]
                            gene_id_subgenome = row_subgenome.get(f"subgenome_{j+1}", 'x')
                            if gene_id_groundtruth == gene_id_subgenome:
                                exact_matches += 1
                            else:
                                missing_genes += 1

                # Calculate and store overlap percentage for this subgenome
                overlap_percentage = exact_matches / total_genes if total_genes > 0 else 0
                results[key_comb].append({
                    'subgenome': j + 1,
                    'exact_matches': exact_matches,
                    'total_genes': total_genes,
                    'missing_genes': missing_genes,
                    'overlap_percentage': overlap_percentage
                })

        return results


    def all_possible_combinations(self, updated_subgenomes, n_sub, swap_sub, GT, subgenome_letters):
        """
        Purpose: Generates two combinations of subgenome assignments for any number of subgenomes:
                1. The original DataFrame
                2. The DataFrame with subgenome_2 and subgenome_3 columns swapped
        Args: 
            updated_subgenomes: A dictionary with keys as identifiers and values as DataFrames representing subgenomes.
            n_sub: The total number of subgenomes.
            swap_sub: The dictionary that you pass to the accuracy function.
            GT: The ground truth file path.
            subgenome_letters: Dictionary mapping subgenome indices to their letter codes.
        Returns: 
            all_combinations: A dictionary with keys as identifiers and values as lists of DataFrames representing the two combinations of subgenome assignments.
            all_results: A dictionary containing accuracy results for each subgenome combination.
        """
        all_combinations = {}
        all_results = {}  # Dictionary to store accuracy results

        for key, df in updated_subgenomes.items():
            # Calculate accuracy for the original DataFrame
            if GT is not None:
                accuracy_before = self.accuracy(GT, n_sub, df, swap_sub, f"{key}_original", subgenome_letters)
                all_results[f"{key}_original"] = accuracy_before
                print(f"Accuracy for {key}_original:")

            if n_sub > 2:
                new_df = df.copy()
                # Swap subgenome_2 and subgenome_3 columns
                new_df['subgenome_2'], new_df['subgenome_3'] = new_df['subgenome_3'].copy(), new_df['subgenome_2'].copy()

                # Calculate accuracy for the swapped version
                if GT is not None:
                    accuracy_after = self.accuracy(GT, n_sub, new_df, swap_sub, f"{key}_after_swap", subgenome_letters)
                    all_results[f"{key}_after_swap"] = accuracy_after
                    print(f"Accuracy for {key}_after_swap:")

            # Save the original and swapped DataFrames to Excel
            try:
                with pd.ExcelWriter(f"{key}_combinations.xlsx", engine='openpyxl') as writer:
                    df.to_excel(writer, sheet_name=f"{key}_original", index=False)
                    if n_sub >= 3:
                        new_df.to_excel(writer, sheet_name=f"{key}_after_swap", index=False)
                print(f"{key}_combinations.xlsx saved successfully.")
            except Exception as e:
                print(f"Error saving Excel file for {key}: {e}")

            # Store the original and swapped DataFrames in a list for each key
            all_combinations[key] = [df]
            if n_sub > 2:
                all_combinations[key].append(new_df)

        print(all_results)

        return all_combinations, all_results




    def find_missing_genes(self, GT_file_path, df_subgenome_dict):
        """
        Purpose: This function finds missing genes in the subgenome assignments compared to the ground truth.
        Args: GT_file_path, df_subgenome_dict
        Returns: missing_genes_dict
        """
        # Read the ground truth DataFrame from the provided file path
        df_groundtruth = pd.read_excel(GT_file_path)
        
        # Initialize a dictionary to store the missing genes for each key
        missing_genes_dict = {}

        # Split the ground truth DataFrame into parts corresponding to each key
        gt_keys = list(df_subgenome_dict.keys())
        gt_split = {gt_keys[0]: df_groundtruth.iloc[:, :3], gt_keys[1]: df_groundtruth.iloc[:, 3:6]}
        # gt_split = {gt_keys[0]: df_groundtruth.iloc[:, :3]}
        # Iterate over each key in the df_subgenome_dict
        for key in gt_keys:
            # Convert the respective part of the ground truth DataFrame to a set of unique genes, excluding 'x'
            gt_genes = set(gt_split[key].replace('x', pd.NA).stack().unique())

            # Check if the key exists in df_subgenome_dict and convert that DataFrame to a set of unique genes, excluding 'x'
            if key in df_subgenome_dict:
                df_genes = set(df_subgenome_dict[key].replace('x', pd.NA).stack().unique())
            else:
                df_genes = set()  # If key doesn't exist in df_subgenome_dict, use an empty set for comparison

            # Find the difference between ground truth genes and DataFrame genes for the current key
            missing_genes = gt_genes - df_genes

            print(len(missing_genes))
            # Store the missing genes for this key
            missing_genes_dict[key] = missing_genes
        
        #save the missing genes in a file
        with open('missing_genes.txt', 'w') as f:
            for key, value in missing_genes_dict.items():
                f.write(f"{key}: {len(value)} missing genes\n")
                f.write(f"{value}\n\n")

        # Return the dictionary containing the missing genes for each key
        return missing_genes_dict
    
