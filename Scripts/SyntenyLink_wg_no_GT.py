import sys
import pandas as pd
import numpy as np
import pickle
import re
import SyntenyLink_acc as acc

def get_blocks(df_subgenome_density, df_synteny, num_blocks, n_subgenomes):
    '''
    This function takes in the subgenome density dataframe, the synteny dataframe, the number of blocks and the number of subgenomes
    and returns a dictionary with the blocks as the keys and the subgenomes as the values
    
    Parameters
    ----------
    df_subgenome_density: pandas dataframe
        The dataframe containing the subgenome density information
        df_synteny: pandas dataframe
        The dataframe containing the synteny information
        num_blocks: int
        The number of blocks
        num_subgenomes: int
        The number of subgenomes
        
    Returns
    -------
    blocks: dictionary
    The dictionary containing the blocks as the keys and the subgenomes as the values
    '''
    # store the blocks in a dictionary with the block number as the key
    blocks = {}
    for i in range(num_blocks):
        blocks[i] = {}
        for j in range(n_subgenomes):
            subgenome_col = "subgenome" + str(j+1)
            block_key = str(i) + "_" + df_subgenome_density[subgenome_col].values[i]
            start_row = df_subgenome_density.loc[i, "Row start #"]
            end_row = df_subgenome_density.loc[i, "Row end #"]
            block_val = df_synteny.loc[start_row:end_row, df_subgenome_density.loc[i, subgenome_col]]
            genes = df_synteny.loc[start_row:end_row, 'gene_id'].tolist()
            block_val_dict = {}
            for k in range(len(genes)):
                sub_key = genes[k]
                block_val_dict[sub_key] = block_val.iloc[k]
            blocks[i][block_key] = block_val_dict
    return blocks
    



def create_chains_dict(df_chains):
    """
    Creates a dictionary with keys as the zeroth column in df_chains and values as the list of values belonging to the same
    keys in the fourth column values in df_chains.
    
    Args:
    df_chains: pandas dataframe containing four columns
    
    Returns:
    chains: dictionary with keys as the zeroth column in df_chains and values as the list of values belonging to the same
    keys in the fourth column values in df_chains.
    """
    chains = {}
    for i in range(len(df_chains[0])):
        chains[df_chains[0][i]] = []
    for i in range(len(df_chains[0])):
        chains[df_chains[0][i]].append(df_chains[2][i])
    return chains

def calculate_density(blocks_dict, num_blocks):
    """
    Calculates the density of each block in each subgenome.
    
    Parameters:
        blocks_dict (dict): A dictionary containing all the blocks and their genes for each subgenome.
        num_blocks (int): The total number of blocks.
    
    Returns:
        dict: A dictionary containing the density of each block for each subgenome.
    """
    density_dict = {}
    for i in range(num_blocks):
        density_dict[i] = {}
        for j in range(len(blocks_dict[i])):
            subgenome_blocks = list(blocks_dict[i].values())[j]
            num_genes = len([gene for gene in subgenome_blocks.values() if gene != "x"])
            block_size = len(subgenome_blocks)
            density_dict[i][list(blocks_dict[i].keys())[j]] = (
                num_genes / block_size if block_size > 0 else 0
            )
    return density_dict

def calculate_pid(blocks, df, num_blocks, pickle_file):
    '''
    Calculates the PID of each block in each subgenome.
    
    Parameters:
        blocks (dict): A dictionary containing all the blocks and their genes for each subgenome.
        df (pandas dataframe): A dataframe containing the PID values for each gene pair.
        num_blocks (int): The total number of blocks.
        pickle_file (str): The path to the pickle file containing the PID values for each gene pair.
        
    Returns:
    dict: A dictionary containing the PID of each block for each subgenome.
    '''
    #Remove dot and number after dot in gene_id_AT
    df["gene_id_AT"] = [re.sub(r"\.\d+", "", gene_id) for gene_id in df["gene_id_AT"]]

    # Initialize the dictionary "pid"
    pid = {}
    for i in range(num_blocks):
        pid[i] = {}
        for j in range(len(blocks[i])):
            key = list(blocks[i].keys())[j]
            block = blocks[i][key]
            pid_values = []
            for sub_key in block.keys():
                if sub_key in df["gene_id_AT"].tolist():
                    val = block[sub_key]
                    # Create a dictionary to store gene_id as keys and PID as values
                    gene_id_to_pid = dict(zip(zip(df["gene_id_Brassica"], df["gene_id_AT"]), df["PID"]))
                    if val!="x":
                        brassica_id = df.loc[(df["gene_id_AT"] == sub_key) & (df["gene_id_Brassica"] == val), "gene_id_Brassica"].values
                        pid_val = gene_id_to_pid.get((brassica_id[0], sub_key), "x") if len(brassica_id) > 0 and val != "x" else "x"
                    else:
                        pid_val = "x"
                    pid_values.append(pid_val)
            pid[i][key] = pid_values
    #save in pickle file
    with open(pickle_file, "wb") as f:
        pickle.dump(pid, f)

    #load from pickle file
    # with open(pickle_file, "rb") as f:
    #     pid = pickle.load(f)
    return pid

def calculate_avg_pid(num_blocks, blocks, pid):
    avg_pid = {}
    for i in range(num_blocks):
        avg_pid[i] = {}
        for j in range(len(blocks[i])):
            avg_pid[i][list(blocks[i].keys())[j]] = 0
            if len([x for x in pid[i][list(blocks[i].keys())[j]] if x != "x"]) != 0:
                avg_pid[i][list(blocks[i].keys())[j]] = (
                    sum([x for x in pid[i][list(blocks[i].keys())[j]] if x != "x"])
                    / len([x for x in pid[i][list(blocks[i].keys())[j]]])
                ) / 100
    return avg_pid


def process_files_weight_cal(chains_file, blastn_file, df_subgenome_density_file_path, df_subgenomes_file_path, synteny_file, pickle_file, num_blocks, n_subgenomes):
    '''
    Processes the files and returns the required dataframes.
    
    Parameters:
        chains_file (str): Path to the chains file.
        blastn_file (str): Path to the blastn file.
        df_subgenome_density_file_path (str): Path to the subgenome density file.
        synteny_file (str): Path to the synteny file.
        pickle_file (str): Path to the pickle file.
        
    Returns:
    blocks_dict: A dictionary containing all the blocks and their genes for each subgenome.
    density_dict: A dictionary containing the density of each block for each subgenome.
    chains_dict: A dictionary containing the chains of each block for each subgenome.
    pid_dict: A dictionary containing the PID of each block for each subgenome.
    avg_pid_dict: A dictionary containing the average PID of each block for each subgenome.
    '''

    # read the blastn file
    blastn = pd.read_csv(blastn_file, sep="\t", header=None)
    df_blastn = blastn.iloc[:, 0:3]
    #have column names for the blastn file as gene_id_Brassica, gene_id_AT, PID
    df_blastn.columns = ["gene_id_Brassica", "gene_id_AT", "PID"]

    # read the chains file
    df_chains = pd.read_csv(chains_file, sep="\t", header=None)

    # get the column names from df_subgenome_density that start with N followed by a number
    df_subgenome = pd.read_excel(df_subgenome_density_file_path)
    column_names = ['gene_id'] + df_subgenome.filter(regex=r'^N\d+').columns.tolist()
    df_subgenome_density = pd.read_excel(df_subgenomes_file_path)

    df_synteny = synteny_file.rename(columns={synteny_file.columns[i]: column_names[i] if i > 0 else "gene_id" for i in range(len(column_names))})

    # Start the index from 0
    df_synteny.index = df_synteny.index - 1

    # number of blocks

    blocks= get_blocks(df_subgenome_density, df_synteny, num_blocks, n_subgenomes)
    density = calculate_density(blocks, num_blocks)
    chains_dict = create_chains_dict(df_chains)
    pid = calculate_pid(blocks, df_blastn, num_blocks, pickle_file)
    avg_pid = calculate_avg_pid(num_blocks, blocks, pid)
    
    return (blocks, density, chains_dict, pid, avg_pid)

def create_graph_input(chains_file, blastn_file, C_df_csv, n_subgenomes, num_blocks, first_letter_get):
    '''
    Creates a list of tuples containing the required dataframes for each subgenome.

    Returns:
        graph_inputs (list): A list of tuples containing the required dataframes for each subgenome.
    '''
    # subgenome_density_files = [f"subgenome_density_bra_sub{i+1}.xlsx" for i in range(n_subgenomes)]
    # pid_genes_removed_files = [f"pid_genes_removed_blastn_sub{i+1}.pickle" for i in range(n_subgenomes)]
    graph_input = process_files_weight_cal(chains_file, blastn_file, "Super_synteny_bl_sub_placement_density.xlsx", "Super_synteny_bl_sub_placement_density.xlsx", C_df_csv.iloc[:,:-3], f"{first_letter_get}_pid_blastn_bna.pickle", num_blocks, n_subgenomes)
    return graph_input

def create_block_graph(blocks, density, chains, avg_pid, num_subgenomes):
    '''
    Creates a graph of blocks.
    
    Parameters:
        blocks (dict): A dictionary containing all the blocks and their genes for each subgenome.
        density (dict): A dictionary containing the density of each block for each subgenome.
        chains (dict): A dictionary containing the chains of each block for each subgenome.
        avg_pid (dict): A dictionary containing the average PID of each block for each subgenome.
        num_subgenomes (int): Number of subgenomes.
    
    Returns:
    graph (dict): A dictionary containing the graph of blocks.
    '''
    
    graph = {}
    num_blocks = len(blocks)
    for i in range(num_blocks - 1):
        for k in range(num_subgenomes):
            graph[list(blocks[i].keys())[k]] = {}
        for s in range(num_subgenomes):
            graph[list(blocks[i + 1].keys())[s]] = {}
        for j in range(len(chains)):
            for l in range(num_subgenomes):
                for m in range(num_subgenomes):
                    if (
                        list(blocks[i].keys())[l].split("_")[1]
                        == list(blocks[i + 1].keys())[m].split("_")[1]
                        or list(blocks[i].keys())[l].split("_")[1].split(".")[0]
                        == list(blocks[i + 1].keys())[l].split("_")[1].split(".")[0]
                    ):
                        graph[list(blocks[i].keys())[l]][
                            list(blocks[i + 1].keys())[m]
                        ] = (
                            list(density[i].values())[l]
                            - (
                                list(density[i].values())[l]
                                - list(density[i + 1].values())[m]
                            )
                        ) + (
                            list(avg_pid[i].values())[l]
                            - (
                                list(avg_pid[i].values())[l]
                                - list(avg_pid[i + 1].values())[m]
                            )
                        ) 
                        new_list = (
                            list(blocks[i][list(blocks[i].keys())[l]].values())
                            + list(blocks[i + 1][list(blocks[i + 1].keys())[m]].values())
                        )
                        new_list = list(filter(lambda a: a != "x", new_list))
                        # concatenate the two blocks

                        if (
                            list(blocks[i].keys())[l].split("_")[1].split(".")[0]
                            == list(chains.keys())[j].split("_")[0]
                            or list(blocks[i].keys())[l].split("_")[1]
                            == list(chains.keys())[j].split("_")[0]
                        ):
                            if (
                                len(set(new_list) & set(
                                    chains[list(chains.keys())[j]]))
                                / float(
                                    len(set(new_list) | set(
                                        chains[list(chains.keys())[j]]))
                                )
                                * 100
                                > 80
                            ):
                                graph[list(blocks[i].keys())[l]][
                                    list(blocks[i + 1].keys())[m]
                                ] = 2
                            elif len(graph[list(blocks[i].keys())[l]]) != num_subgenomes:
                                graph[list(blocks[i].keys())[l]][
                                    list(blocks[i + 1].keys())[m]
                                ] = (
                                    list(density[i].values())[l]
                                    - (
                                        list(density[i].values())[l]
                                        - list(density[i + 1].values())[m]
                                    )
                                ) + (
                                    list(avg_pid[i].values())[l]
                                    - (
                                        list(avg_pid[i].values())[l]
                                        - list(avg_pid[i + 1].values())[m]
                                    )
                                )

                    elif (
                        list(blocks[i].keys())[l].split("_")[1]
                        == list(blocks[i + 1].keys())[m].split("_")[1].split(".")[0]
                        or list(blocks[i + 1].keys())[l].split("_")[1]
                        == list(blocks[i].keys())[m].split("_")[1].split(".")[0]
                    ):
                        graph[list(blocks[i].keys())[l]][
                            list(blocks[i + 1].keys())[m]
                        ] = (
                            list(density[i].values())[l]
                            - (
                                list(density[i].values())[l]
                                - list(density[i + 1].values())[m]
                            )
                        ) + (
                            list(avg_pid[i].values())[l]
                            - (
                                list(avg_pid[i].values())[l]
                                - list(avg_pid[i + 1].values())[m]
                            )
                        ) 
                    else:
                        graph[list(blocks[i].keys())[l]][list(blocks[i + 1].keys())[m]] = (
                            list(density[i].values())[l]
                            - (
                                list(density[i].values())[l]
                                - list(density[i + 1].values())[m]
                            )
                        ) + (
                            list(avg_pid[i].values())[l]
                            - (
                                list(avg_pid[i].values())[l]
                                - list(avg_pid[i + 1].values())[m]
                            )
                        )
                                
    return graph


# function for finding gaps in each columns

def add_edge_weights(blocks, graph, weight, num_blocks):
    '''
    This function adds edge weights to the graph
    
    Parameters
    ----------
    blocks : list
        list of blocks
    graph : dict
        graph of blocks
    weight : int
        weight of the edge
    
    Returns
    -------
    graph : dict
        graph of blocks with edge weights
    '''
    for i in range(num_blocks - 1):
        for node in blocks[i].keys():
            max_value = -1
            max_next_node = None
            for next_node in blocks[i + 1].keys():
                if (
                    node.split("_")[1] == next_node.split("_")[1]
                    or node.split("_")[1] == next_node.split("_")[1].split(".")[0]
                    or next_node.split("_")[1] == node.split("_")[1].split(".")[0]
                ):
                    if graph[node][next_node] > max_value:
                        max_value = graph[node][next_node]
                        max_next_node = next_node
            if max_next_node is not None:
                current_edge_value = max_value
                new_edge_value = current_edge_value + weight
                graph[node][max_next_node] = new_edge_value

    return graph


def get_nodes(start_node, end_node, num_blocks, graph):
    '''
    This function returns the nodes between the start and end node
    
    Parameters
    ----------
    start_node : str
        start node
    end_node : str
        end node
    num_blocks : int
        number of blocks
    graph : dict
        graph of blocks
    
    Returns
    -------
    nodes_sub : list
        list of nodes between the start and end node
    '''
    nodes_sub = []
    for i in range(num_blocks - 1):
        while start_node != end_node:
            nodes_sub.append(start_node)
            start_node = list(list(graph.values())[list(graph.keys()).index(start_node)])[
                list(
                    list(graph.values())[
                        list(graph.keys()).index(start_node)].values()
                ).index(
                    max(
                        list(
                            list(graph.values())[
                                list(graph.keys()).index(start_node)
                            ].values()
                        )
                    )
                )
            ]
            break
    nodes_sub.append(end_node)
    return nodes_sub


def remove_nodes_from_graph(nodes_sub, graph, blocks):
    '''
    This function removes the nodes from the graph
    
    Parameters
    ----------
    nodes_sub : list
        list of nodes between the start and end node
    graph : dict
        graph of blocks
    blocks : list
        list of blocks
    
    Returns
    -------
    graph : dict
        graph of blocks with nodes removed
    blocks : list
        list of blocks with nodes removed
    '''
    for i in range(len(nodes_sub)):
        graph.pop(nodes_sub[i], None)
        for j in range(len(graph)):
            graph[list(graph.keys())[j]].pop(nodes_sub[i], None)
    # remove blocks in blocks[] matching the nodes in nodes[]
    for i in range(len(nodes_sub)):
        for j in range(len(blocks)):
            blocks[j].pop(nodes_sub[i], None)
    return graph, blocks


def update_graph_edges(num_blocks, blocks, graph, weight):
    '''
    This function updates the graph edges
    
    Parameters
    ----------
    num_blocks : int
        number of blocks
    blocks : list
        list of blocks
    graph : dict
        graph of blocks
    weight : int
        weight of the edge
    
    Returns
    -------
    graph : dict
        graph of blocks with updated edges
    '''
    for i in range(num_blocks - 1):
        for node in blocks[i].keys():
            max_value = -1
            max_next_node = None
            for next_node in blocks[i + 1].keys():
                if (
                    node.split("_")[1] == next_node.split("_")[1]
                    or node.split("_")[1] == next_node.split("_")[1].split(".")[0]
                    or next_node.split("_")[1] == node.split("_")[1].split(".")[0]
                ):
                    if graph[node][next_node] > max_value:
                        max_value = graph[node][next_node]
                        max_next_node = next_node
            if max_next_node is not None:
                current_edge_value = max_value
                new_edge_value = current_edge_value - weight
                graph[node][max_next_node] = new_edge_value
    return graph

def node_traverse_most_weighted_path(n_subgenomes, df_synteny, chains_file, blastn_file, C_df_csv, num_blocks, first_letter_get):
    '''
    This function traverses the graph to find the most weighted path
    
    Parameters
    ----------
    n_subgenomes : int
        number of subgenomes
        
    Returns
    -------
    nodes_df : pandas dataframe
        dataframe with the nodes_sub lists as columns
    '''
    subgenome_density_files = "Super_synteny_bl_sub_placement_density.xlsx"
    # Create graph with edge weights and get nodes for first subgenome
    # Get input from command line
    for k in range(n_subgenomes):
        # if k == 0:
        #     weight_1 = float(sys.argv[sys.argv.index('-w1s1') + 1])
        #     weight_2 = float(sys.argv[sys.argv.index('-w2s1') + 1])
        # if k == 1:
        #     weight_1 = float(sys.argv[sys.argv.index('-w1s2') + 1])
        #     weight_2 = float(sys.argv[sys.argv.index('-w2s2') + 1])
        # if k == 2:
        #     weight_1 = float(sys.argv[sys.argv.index('-w1s3') + 1])
        #     weight_2 = float(sys.argv[sys.argv.index('-w2s3') + 1])
        # Get input from files
        subgenome_density = pd.read_excel(subgenome_density_files)
        graph_input = create_graph_input(chains_file, blastn_file, C_df_csv, n_subgenomes, num_blocks, first_letter_get)
        graph = create_block_graph(graph_input[0], graph_input[1], graph_input[2], graph_input[4], n_subgenomes)
        #get the first value in subgenome1 column (first row) in subgenome_density file as start node
        start_node = f"0_{subgenome_density.iloc[0]['subgenome1']}"
        #get the last value in subgenome1 column (last row) in subgenome_density file as end node
        end_node = f"{len(subgenome_density)-1}_{subgenome_density.iloc[-1]['subgenome1']}"
        nodes_sub1 = get_nodes(start_node, end_node, len(graph_input[0]), graph)
        # Remove nodes from graph and get nodes for remaining subgenomes
        blocks = graph_input[0]
        nodes_df = pd.DataFrame({"nodes_sub1": nodes_sub1}) # Define the nodes_df variable here
        for i in range(n_subgenomes - 1):
            graph, blocks = remove_nodes_from_graph(nodes_sub1, graph, blocks)
            # graph = update_graph_edges(len(blocks), blocks, graph, weight_2)
            start_node = f"0_{subgenome_density.iloc[0][f'subgenome{i+2}']}"
            end_node = f"{len(subgenome_density)-1}_{subgenome_density.iloc[-1][f'subgenome{i+2}']}"
            nodes_sub = get_nodes(start_node, end_node, len(blocks), graph)
            nodes_sub_name = f"nodes_sub{i+2}"
            nodes_df[nodes_sub_name] = nodes_sub
            remove_nodes_from_graph(nodes_sub, graph, blocks)

        # rename the columns as subgenomes
        nodes_df.columns = [f"subgenome{i}" for i in range(1, n_subgenomes+1)]

        #Append Row start # and Row end # columns from subgenome_density_bra as first two columns of the dataframe
        nodes_df.insert(0, "Row start #", subgenome_density["Row start #"])
        nodes_df.insert(1, "Row end #", subgenome_density["Row end #"])

        #remove everything before _ in the nodes_df subgenome columns
        for i in range(1, n_subgenomes+1):
            column_name = f"subgenome{i}"
            nodes_df[column_name] = nodes_df[column_name].str.split("_", n=1, expand=True)[1]

        nodes_df.to_excel(f"Super_synteny_graph_nodes_sub{k+1}.xlsx", index=False)
        # acc.subgenome_overlap_separate(GT,f"Super_synteny_graph_nodes_sub{k+1}.xlsx", df_synteny, 3, k, first_letter_get)

