import os
import pickle
import networkx as nx
import pandas as pd

# import matplotlib.pyplot as plt


interactome_file_name = "data/BIOGRID-ALL-4.4.240.tab3.txt"  # starter file
seeds_file_name = "data/DISEASES_Summary_GDA_CURATED_C0025202.tsv"
human_ID = 9606
interaction_of_interest = "physical"
protein_A_identifier = "Entrez Gene Interactor A"
protein_B_identifier = "Entrez Gene Interactor B"
PPI_list_variable = "data/PPI_of_interest.pk1"
first_gene_symbol_indicator = "Official Symbol Interactor A"
second_gene_symbol_indicator = "Official Symbol Interactor B"


# Importing and Processing the Interactome

if os.path.exists(PPI_list_variable):
    print("The file already exists! Importing information from the file ...")
    with open(PPI_list_variable, 'rb') as file:
        PPI = pickle.load(file)

else:
    print("The file does not exist yet. Processing... ")
    complete_interactome = pd.read_csv(interactome_file_name, delimiter='\t')

    human_interactome = complete_interactome[
        (complete_interactome["Organism ID Interactor A"] == human_ID)
        & (complete_interactome["Organism ID Interactor B"] == human_ID)]

    human_physical_interactome = human_interactome[
        (human_interactome["Experimental System Type"] == interaction_of_interest)]

    PPI_NoDuplicates = human_physical_interactome.drop_duplicates(
        subset=[first_gene_symbol_indicator, second_gene_symbol_indicator])

    PPI_NoDuplicates_NoSelfLoop = PPI_NoDuplicates[
        PPI_NoDuplicates[first_gene_symbol_indicator]
        !=
        PPI_NoDuplicates[second_gene_symbol_indicator]]
    # Note that with this second version,
    # the number of edges are "slightly" different
    # (9M with the first version vs 8.7 with the second version)
    # I think I have to review Database and GeneSymbols

    PPI = PPI_NoDuplicates_NoSelfLoop
    with open(PPI_list_variable, 'wb') as file:
        pickle.dump(PPI, file)

    del complete_interactome, human_interactome, human_physical_interactome, \
        PPI_NoDuplicates, PPI_NoDuplicates_NoSelfLoop

# %% Create the Interactome LCC
graph_of_PPI = nx.from_pandas_edgelist(PPI, protein_A_identifier, protein_B_identifier)
print(graph_of_PPI)

## To verify
connected_components = nx.connected_components(graph_of_PPI)

largest_cc = max(connected_components, key=len)

LCC_subgraph = nx.subgraph(graph_of_PPI, largest_cc).copy()
print(LCC_subgraph)

# %% Importing and the seed genes
# Note that this will not be the final version, because
Seed_genes_df = pd.read_csv(seeds_file_name, sep='\t')
Seed_genes_Symbols = Seed_genes_df["Gene"]

Disease_Interactome = PPI[
    (PPI[first_gene_symbol_indicator].isin(Seed_genes_Symbols))
    &
    (PPI[second_gene_symbol_indicator].isin(Seed_genes_Symbols))
    ]

Graph_of_disease_interactome = nx.from_pandas_edgelist(
    Disease_Interactome,
    first_gene_symbol_indicator, second_gene_symbol_indicator)

nx.draw(Graph_of_disease_interactome)
# there are some self-loops, we have to debug the initial list
# processing of the PPI; maybe the fact that this graph is determined by the
# official gene symbols could have changed something.
# Yes, definitely with this new version we don't have self loops no more
