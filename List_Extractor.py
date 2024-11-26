# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:03:05 2024

@author: atara
"""

import pandas as pd
import networkx as nx
# import matplotlib.pyplot as plt


human_ID = 9606
interaction_of_interest = "physical"
Protein_A_identifier = "Entrez Gene Interactor A"
Protein_B_identifier = "Entrez Gene Interactor B"


complete_interactome = pd.read_csv("BIOGRID-ALL-4.4.239.tab3.txt", 
                                   delimiter='\t')
human_interactome = complete_interactome[
    (complete_interactome["Organism ID Interactor A"] == human_ID)
    & (complete_interactome["Organism ID Interactor B"] == human_ID)]

human_physical_interactome = human_interactome[
    (human_interactome["Experimental System Type"] == interaction_of_interest)]

## human_physical_interactome.to_csv("human_physical_interactome.csv", index = False)

# Duplicate Removal
PPI_NoDuplicates = human_physical_interactome.drop_duplicates(
    subset = [Protein_A_identifier, Protein_B_identifier])

PPI_NoDuplicates_NoSelfLoop = PPI_NoDuplicates[
    PPI_NoDuplicates[Protein_A_identifier] != PPI_NoDuplicates[Protein_B_identifier]]

PPI = PPI_NoDuplicates_NoSelfLoop
del complete_interactome, human_interactome, human_physical_interactome, PPI_NoDuplicates, \
    PPI_NoDuplicates_NoSelfLoop

PPI.to_csv("PPI_of_Interest.csv", index = False)

Graph_of_PPI = nx.from_pandas_edgelist(PPI, Protein_A_identifier, Protein_B_identifier)
print(Graph_of_PPI)


## To verify
Connected_Components = nx.connected_components(Graph_of_PPI)
 
largest_cc= max(Connected_Components, key=len)

LCC_subgraph = nx.subgraph(Graph_of_PPI, largest_cc).copy()
print(LCC_subgraph)


