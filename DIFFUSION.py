# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 18:50:13 2024

@author: atara
"""

# import operator
# import copy
import networkx as nx
from numpy import array
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import  expm_multiply
import matplotlib.pyplot as plt



def diffuse(matrix, heat_array, time):
    return expm_multiply(-matrix, heat_array,
                         start=0, stop = time, endpoint=True)[-1]

def sparse_laplacian(network, normalize = False):
    if normalize:
        return csc_matrix(nx.normalized_laplacian_matrix(network))
    else:
        return csc_matrix(nx.laplacian_matrix(network))

def create_heat_array(network, seed_genes, heat_value=1.0):
    heat_list = []
    seed_genes = list(seed_genes)
    for node in network.nodes:
        if node in seed_genes:
            heat_list.append(heat_value)
        else:
            heat_list.append(0.0)

    return array(heat_list)

def filter_node_list(node_list, nodes_to_remove):
    nodes_to_remove = list(nodes_to_remove)
    filtered_nodes = []
    for item in node_list:
        if item[0] not in nodes_to_remove:
            filtered_nodes.append(item[0])

    return filtered_nodes

def run_heat_diffusion(network, seed_genes, diffusion_time = 0.005, n_positions=None):
    # Get the sparse Lapliacian matrix of the network.
    matrix = sparse_laplacian(network, normalize=True)

    # Create heat array.
    heat_array = create_heat_array(network, seed_genes)

    # Diffuse heat.
    diffused_heat_array = diffuse(matrix, heat_array, diffusion_time)

    # Get the heat of each node.
    node_heat = {node_name: diffused_heat_array[i] 
                 for i, node_name in enumerate(network.nodes())}

    # Sort the nodes by heat value.
    sorted_nodes = sorted(node_heat.items(), key=lambda x:x[1], reverse=True)

    # Remove from the sorted nodes the seed genes.
    predicted_genes = filter_node_list(sorted_nodes, seed_genes)

    if n_positions:
        # Return the top `n_positions` genes
        # wrt their heat value.
        return predicted_genes[:n_positions]

    return predicted_genes

def mapping_the_diffusion(network, seed_genes, predicted_genes):
    # plt.show()
    pos = nx.kamada_kawai_layout(network)
    # nx.draw_networkx_nodes(network, pos, node_size = 10, alpha=0.7)    
    # nx.draw_networkx_edges(network, pos, width = 0.5, alpha = 0.5)
    # plt.tight_layout()
    
    # plt.show()
    # color_map = ['blue' if node in seed_genes
    #              else 'orange' if node in predicted_genes
    #              else 'blue' for node in seed_genes]
    # nx.draw_networkx_nodes(network, pos, node_color = color_map, 
    #                    node_size = 10, alpha=0.7)    
    # nx.draw_networkx_edges(network, pos, width = 0.5, alpha = 0.5)
    # plt.tight_layout()    
    plt.show()
    nx.draw_networkx_nodes(network, pos, nodelist = list(network.nodes), 
                           node_color = "blue", node_size = 25, 
                           edgecolors = "blue")
    nx.draw_networkx_nodes(network, pos, nodelist = list(seed_genes), 
                           node_color = "red", node_size = 25, 
                           edgecolors = "red")
    nx.draw_networkx_nodes(network, pos, nodelist = list(predicted_genes), 
                           node_color = "orange", node_size = 20, 
                           edgecolors = "orange")
    nx.draw_networkx_edges(network, pos, width = 0.5, alpha = 0.5)
    plt.tight_layout()    


#%%
# # to call the algorithm
# predicted_genes = run_heat_diffusion(
#     network=PPI_LCC,
#     seed_genes=seeds,
#     diffusion_time=0.05,
#     n_positions=50)