import mainbreakpoint_ADT as MBP
import pandas as pd
import copy as cp
import sys
from itertools import combinations

class Graph(object):
    def __init__(self):
        """
        Purpose: Initialize the Graph class object
        Args: None
        Returns: None
        """
        self.synteny = MBP.BP()

    def main(self):
        """
        Purpose: This is the main function for the graph analysis. It will call the main function from the graph_ADT.py file and run the main function.
        Args: None
        Returns: 
        """
        # Get the input file as an argument (blastn dagchainer output file)
        dag = sys.argv[sys.argv.index('-dag') + 1]

        # Open and read the DAG file
        with open(dag, 'r') as f:
            dag_content = f.read()


        # Run the main function
        n_sub, ploidy_numpy, collinear, gap_threshold, minimum_block_length, file, ploidy_status, chrsets, header, header_filtered, GT, empty_df, ref_gene_id, bed_file, ploidy_numpy_overlap, gene_prefix, subgenome_letters = self.synteny.main()

        # Create a copy of the original ploidy numpy array
        original_ploidy_numpy = cp.deepcopy(ploidy_numpy)

        # Run the main breakpoint generation function
        main_breakpoints = self.synteny.gap_calculation(ploidy_numpy, gap_threshold, minimum_block_length, n_sub)

        # Get the densities of the main breakpoints
        density = self.synteny.get_densities_bp(main_breakpoints, ploidy_numpy)

        # Get the noisy breakpoints
        noisy = self.synteny.get_noisy_breakpoints(density, n_sub)

        # Run the smaller breakpoints function
        subgenomes = self.synteny.smaller_breakpoints(main_breakpoints, noisy, n_sub, ploidy_numpy)

        # Get the densities of the smaller breakpoints
        density_smallbp = self.synteny.get_densities_bp(subgenomes, ploidy_numpy)

        # Get the noisy breakpoints in the smaller breakpoints
        noisy_smallbp = self.synteny.get_noisy_breakpoints(density_smallbp, n_sub)

        # Get the subgenome assignment
        subgenome_assignment = self.synteny.get_subgenomes_unstranded(noisy_smallbp, density_smallbp, n_sub)

        # Adjust the subgenomes based on continuity
        adjusted_subgenomes = self.synteny.adjust_subgenome_unstranded(subgenome_assignment, n_sub)

        # Get the subgenomes stranded converting forward and reverse pair to stranded
        subgenomes_stranded = self.synteny.get_subgenomes_stranded(adjusted_subgenomes, file, n_sub, ploidy_status, chrsets, bed_file, gene_prefix)

        # Get the sliced collinear data for each key
        # Generate the DataFrame from update_overlap_genes
        updated_overlap_genes_df = pd.DataFrame(self.synteny.synteny.update_overlap_genes(self.synteny.synteny.create_collinear_all(self.synteny.synteny.create_matrix(file)), bed_file, gene_prefix))

        # Extract the last two columns
        last_two_columns = updated_overlap_genes_df.iloc[:, -3:-1]

        last_column = updated_overlap_genes_df.iloc[:, -1]

        # Assuming 'ploidy_status', 'chrsets', and 'file' are defined and accessible
        # Generate the 'gene_id_matrix_df_' DataFrame
        gene_id_matrix_df_ = self.synteny.synteny.slice_collinear(ploidy_status, chrsets, pd.DataFrame(self.synteny.synteny.update_overlap_genes(self.synteny.synteny.create_collinear_all(self.synteny.synteny.create_matrix(file)), bed_file, gene_prefix)))

        # Adjust the subgenomes stranded based on continuity
        adjusted_subgenomes_stranded =self.synteny.adjust_subgenome_stranded(subgenomes_stranded, adjusted_subgenomes, collinear, n_sub, noisy_smallbp, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns)

        # Create a copy of the adjusted subgenomes stranded
        original_adjusted_subgenomes_stranded = cp.deepcopy(adjusted_subgenomes_stranded)

        adjusted_subgenomes_stranded_segmental = self.synteny.adjust_subgenome_stranded_segmental(adjusted_subgenomes_stranded , adjusted_subgenomes, collinear, n_sub, noisy_smallbp, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns, empty_df, ref_gene_id)

        #copy the adjusted subgenomes stranded segmental
        adjusted_subgenomes_stranded_segmental_original = cp.deepcopy(adjusted_subgenomes_stranded_segmental)

        # Get the still missing genes
        # self.synteny.find_missing_genes(GT, adjusted_subgenomes_stranded_segmental)

        # print(adjusted_subgenomes_stranded_segmental)

        # self.synteny.accuracy_gene_id(GT, n_sub, adjusted_subgenomes_stranded_segmental)

        return adjusted_subgenomes_stranded_segmental_original, dag_content, GT, n_sub, header, ploidy_status, chrsets, file, gene_id_matrix_df_, last_two_columns, original_ploidy_numpy, last_column, ploidy_numpy_overlap, subgenome_letters
    
   # Parse synteny blocks
    def parse_synteny_blocks(self, synteny_blocks):
        """
        Purpose: Parse the synteny blocks from the DAG file
        Args:
            synteny_blocks (str): The synteny blocks from the DAG file
            Returns:
            block_dict (dict): A dictionary containing the synteny blocks
        """
        block_dict = {}
        current_block = None
        for line in synteny_blocks.strip().split('\n'):
            if line.startswith('## alignment'):
                current_block = line.strip()
                block_dict[current_block] = []
            else:
                block_dict[current_block].append(line.strip().split('\t'))
        return block_dict
    
    # Process the synteny blocks
    def process_synteny_blocks(self, block_dict):
        """
        Purpose: Process the synteny blocks to get the gene pairs
        Args:
            block_dict (dict): A dictionary containing the synteny blocks
            Returns:
            synteny_blocks (dict): A dictionary containing the gene pairs
        """
        synteny_blocks = {}
        block_num = 0
        for block, genes in block_dict.items():
            block_num += 1
            synteny_blocks[block_num] = []
            for gene in genes:
                ref_gene = gene[1].split('.')[0]
                query_gene = gene[5]
                synteny_blocks[block_num].append((ref_gene, query_gene))
        return synteny_blocks
                
    def weight_graph(self, adjusted_subgenomes_stranded_segmental):
        """
        Purpose: Weight the graph based on synteny blocks
        Args:
            adjusted_subgenomes_stranded_segmental (dict): The adjusted subgenomes with stranded segmental information
            Returns:
            final_gene_combinations (dict): The final gene combinations for each subgenome
        """
        final_gene_combinations = {}
        for key, df in adjusted_subgenomes_stranded_segmental.items():
            gene_combinations = {}
            for i, row in df.iterrows():
                for j in range(1, len(df.columns)):
                    if j not in gene_combinations:
                        gene_combinations[j] = []
                    gene_combinations[j].append((row[0], row[j]))
            final_gene_combinations[key] = gene_combinations
        return final_gene_combinations
    
    def create_graphs(self, final_gene_combinations, synteny_blocks):
        """
        Purpose: Create graphs for each subgenome
        Args:
            final_gene_combinations (dict): The final gene combinations for each subgenome
            synteny_blocks (dict): The synteny blocks
            Returns:
            graphs (dict): A dictionary containing the graphs for each subgenome
        """
        graphs = {}
        for key, gene_combinations in final_gene_combinations.items():
            # Convert gene_combinations to a DataFrame
            df = pd.DataFrame(gene_combinations).T
            num_layers = df.shape[1]
            num_subgenomes = df.shape[0]
            
            G = {}
            
            # Add all nodes to the graph
            for layer in range(num_layers):
                for subgenome in range(num_subgenomes):
                    node = df.iloc[subgenome, layer][1]
                    if node ==  'x':
                        # Modify the node name to add subgenome number and layer number
                        node = f'x_subgenome_{subgenome+1}_layer_{layer}'
                    G[node] = {'layer': layer, 'edges': []}
            
            # Save the graph as an excel file
            df_graph = pd.DataFrame(G).T
            df_graph.to_excel(f'graph_nodes_{key}.xlsx')
            
            # Add edges between nodes of consecutive layers
            for layer in range(1, num_layers):
                prev_layer_nodes = []
                current_layer_nodes = []
                
                for subgenome in range(num_subgenomes):
                    sub = subgenome + 1
                    if df.iloc[subgenome, layer - 1][1] == 'x':
                        prev_node = f'x_subgenome_{sub}_layer_{layer-1}'
                    else:
                        prev_node = df.iloc[subgenome, layer - 1][1]
                    prev_layer_nodes.append(prev_node)
                    
                for subgenome in range(num_subgenomes):
                    sub = subgenome + 1
                    if df.iloc[subgenome, layer][1] == 'x':
                        current_node = f'x_subgenome_{sub}_layer_{layer}'
                    else:
                        current_node = df.iloc[subgenome, layer][1]
                    current_layer_nodes.append(current_node)
                    
                prev_layer_refs = [df.iloc[subgenome, layer - 1][0] for subgenome in range(num_subgenomes)]
                current_layer_refs = [df.iloc[subgenome, layer][0] for subgenome in range(num_subgenomes)]
                
                for current_node, current_ref_gene in zip(current_layer_nodes, current_layer_refs):
                    for prev_node, prev_ref_gene in zip(prev_layer_nodes, prev_layer_refs):
                        weight = 1.0  # Default weight
                        for block_num, genes in synteny_blocks.items():
                            if (prev_ref_gene, prev_node) in genes and (current_ref_gene, current_node) in genes:
                                # Check continuity: ensure prev_ref_gene is followed by current_ref_gene
                                prev_index = genes.index((prev_ref_gene, prev_node))
                                current_index = genes.index((current_ref_gene, current_node))
                                if current_index == prev_index + 1:
                                    weight = 10.0  # Increase weight if there is a synteny block match and continuity

                        G[prev_node]['edges'].append((current_node, weight))
            
            graphs[key] = G
            
            # Save the graph as an excel file
            df_graph = pd.DataFrame(G).T
            df_graph.to_excel(f'graph_{key}.xlsx')

        return graphs


    def traverse_highest_weighted_path(self, graph_data, key, n_sub):
        """
        Purpose: Traverse the highest weighted path in the graph
        Args:
            graph_data (dict): The graph data containing the nodes and edges
            key (int): The key for the graph data
            n_sub (int): The number of subgenomes
            Returns:
            df (DataFrame): The DataFrame containing the subgenome paths
        """
        subgenome_paths = {i: [] for i in range(1, n_sub + 1)}
        used_nodes = set()

        # Get unique start nodes for each subgenome path
        start_nodes = list(graph_data.keys())[:n_sub]

        for i, start_node in enumerate(start_nodes, 1):
            if start_node in used_nodes:
                continue  # Skip if node is already used in another subgenome

            current_node = start_node
            subgenome_paths[i].append(current_node)
            used_nodes.add(current_node)

            # Build path for the current subgenome
            while True:
                next_node = None
                max_weight = -1

                if current_node is not None:
                    for edge, weight in graph_data[current_node]['edges']:
                        if edge not in used_nodes and (weight > max_weight or next_node is None):
                            max_weight = weight
                            next_node = edge

                if next_node is None or next_node in used_nodes:
                    break  # Stop if no valid next node or node already used

                subgenome_paths[i].append(next_node)
                used_nodes.add(next_node)
                current_node = next_node

        # Save to Excel
        df = pd.DataFrame({f"subgenome_{k}": pd.Series(v) for k, v in subgenome_paths.items()})
        df.replace(to_replace=r'^x_.*', value='x', regex=True, inplace=True)
                # Get the still missing genes
        # self.synteny.find_missing_genes(GT, df)

        # print(adjusted_subgenomes_stranded_segmental)

        # self.synteny.accuracy_gene_id_new(GT, n_sub, df, key)
        df.to_excel(f'paths_{key}.xlsx', index=True)

        print("Excel file with subgenome paths has been created.")

        return df
