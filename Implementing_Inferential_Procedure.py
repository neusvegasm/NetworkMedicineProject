# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 21:59:57 2024

@author: atara
"""
# import os
# import csv
import pandas as pd
import pickle
import networkx as nx
from DIAMOnD import DIAMOnD
from DiaBLE import DiaBLE


Graph_of_PPI_file_name = "Graph_of_PPI.pickle"
Graph_of_Disease_Interactome_file_name = "Graph_of_Disease_Interactome.pickle"
Graph_of_LCC_of_disease_interactome_file_name = "Graph_of_LCC_of_disease_interactome.pickle"
Seed_genes_Symbol_file_name = "Seed_genes_Symbol.pickle"
PPI_list_variable_corrected_updated = "PPI_of_interest_corrected updated"
number_of_nodes_inferred = 100
node_inferred_by_DIAMOnD_file = 'Node Inferred by DIAMOnD.pickle'
node_inferred_by_DiaBLE_file = 'Node Inferred by DiaBLE'
First_gene_symbol_indicator = "Official Symbol Interactor A"
Second_gene_symbol_indicator = "Official Symbol Interactor B"
Top_positions = [25, 50, 65, 130, 252] 




with open(Graph_of_PPI_file_name, "rb") as file:
    Graph_of_PPI = pickle.load(file)
with open(Graph_of_Disease_Interactome_file_name, "rb") as file:
    Graph_of_Disease_Interactome = pickle.load(file)
with open(Graph_of_LCC_of_disease_interactome_file_name, "rb") as file:
    Graph_of_LCC_of_disease_interactome = pickle.load(file)
Seed_genes_Symbols = pd.read_pickle(Seed_genes_Symbol_file_name)
with open(PPI_list_variable_corrected_updated, 'rb') as file:
    PPI = pickle.load(file)
    
    
#%% Toy Simulation of DiaBLE functioning 

import random
from DiaBLE import get_neighbors_and_degrees
from DiaBLE import GraphVisualizerForDiable

Toy_Seed_List = Seed_genes_Symbols[0:9]
Toy_Gene_List = Seed_genes_Symbols[0:120]
Toy_Graph = Graph_of_PPI.subgraph(Toy_Gene_List) 
# maybe for the extraction of the disease interactome I should have used this 
# instruction instead of what I actually used in List_Extractor.
# Fun fact using this type of extraction we won't have no more,
# the genes of the seed list that are ignored by both DIAMOnD and DiaBLE

neighbors, all_degrees = get_neighbors_and_degrees(Toy_Graph)
cluster_nodes = set(Toy_Seed_List)
All_Genes = set(Toy_Gene_List)
neighbors_of_the_seeds = set()

for node in cluster_nodes:
    neighbors_of_the_seeds |= neighbors[node]
# not_in_cluster.discard('MELTF')
seeds_and_neighbors = neighbors_of_the_seeds.copy()
neighbors_of_the_seeds -= cluster_nodes
    
# initialization of DiaBLE Universe
Diable_Universe = set()
for node in seeds_and_neighbors:
    Diable_Universe |= neighbors[node]
new_neighbors = Diable_Universe.copy()
new_neighbors -= seeds_and_neighbors


Total_nodes_to_add = 10;
added_nodes = []

GraphVisualizerForDiable(Toy_Graph, cluster_nodes,
                          neighbors_of_the_seeds, 
                          new_neighbors)

while len(added_nodes) < Total_nodes_to_add:
    next_node = random.choice(list(neighbors_of_the_seeds))
    cluster_nodes.add(next_node)
    added_nodes.append(next_node)
    
    for node in cluster_nodes:
        neighbors_of_the_seeds |= neighbors[node]
    # not_in_cluster.discard('MELTF')
    seeds_and_neighbors = neighbors_of_the_seeds.copy()
    neighbors_of_the_seeds -= cluster_nodes
        
    # initialization of DiaBLE Universe
    for node in seeds_and_neighbors:
        Diable_Universe |= neighbors[node]
    new_neighbors = Diable_Universe.copy()
    new_neighbors -= seeds_and_neighbors
    
    GraphVisualizerForDiable(Toy_Graph, cluster_nodes,
                              neighbors_of_the_seeds, 
                              new_neighbors)
    
#%% Node inference 

list_of_nodes_added_by_DIAMOnD = DIAMOnD(Graph_of_PPI, Seed_genes_Symbols, 
                                         number_of_nodes_inferred, alpha = 1, 
                                         outfile = node_inferred_by_DIAMOnD_file)

list_of_nodes_added_by_DiaBLE = DiaBLE(Graph_of_PPI, Seed_genes_Symbols, 
                                       number_of_nodes_inferred, alpha = 1, 
                                       outfile = node_inferred_by_DiaBLE_file)


nodes_added_by_DIAMOnD = [Generic_tuple[0] 
                         for Generic_tuple in list_of_nodes_added_by_DIAMOnD]
nodes_added_by_DiaBLE = [Generic_tuple[0] 
                         for Generic_tuple in list_of_nodes_added_by_DiaBLE]

Seeds_and_DIAMOnD_Inferred_Seeds = \
   pd.concat([Seed_genes_Symbols, pd.Series(nodes_added_by_DIAMOnD)], 
             ignore_index=True)
   


#%% Network after the DIAMOnD Inference
Disease_Interactome_linked_Actual_and_Inferred_seeds = PPI[
    (PPI[First_gene_symbol_indicator].isin(Seeds_and_DIAMOnD_Inferred_Seeds)) 
    &
    (PPI[Second_gene_symbol_indicator].isin(Seeds_and_DIAMOnD_Inferred_Seeds))
    ] 
   
Seed_genes_with_links = list(
    set(Disease_Interactome_linked_Actual_and_Inferred_seeds[
        First_gene_symbol_indicator
        ].unique()).union(
        Disease_Interactome_linked_Actual_and_Inferred_seeds[
            Second_gene_symbol_indicator
            ].unique())
    )
Seed_genes_without_links = list(
    set(Seed_genes_with_links).difference(Seed_genes_Symbols))

Graph_of_inferred_connected_components_of_disease_interactome = nx.from_pandas_edgelist(
    Disease_Interactome_linked_Actual_and_Inferred_seeds, 
    First_gene_symbol_indicator, Second_gene_symbol_indicator)
print(Graph_of_inferred_connected_components_of_disease_interactome)
nx.draw_networkx(Graph_of_inferred_connected_components_of_disease_interactome, 
                 with_labels = False, node_size = 10, width = 0.5)




#%% 


    
    







    
