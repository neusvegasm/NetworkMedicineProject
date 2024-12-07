import os
import pickle
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# Variables ----------------------------------------------------------------------------------------------

interactome_file_name = "data/BIOGRID-ORGANISM-Homo_sapiens-4.4.240.tab3.txt"  # PPI FILE
seeds_file_name = "data/DISEASES_Summary_GDA_CURATED_C0025202.tsv"
human_ID = 9606
interaction_of_interest = "physical"
protein_A_identifier = "Entrez Gene Interactor A"
protein_B_identifier = "Entrez Gene Interactor B"
PPI_list_variable = "data/PPI_of_interest.pk1"
first_gene_symbol_indicator = "Official Symbol Interactor A"
second_gene_symbol_indicator = "Official Symbol Interactor B"


# Importing and Processing Files ----------------------------------------------------------------


## File already preprocessed
if os.path.exists(PPI_list_variable):
    print("The file already exists! Importing information from the file ...")
    with open(PPI_list_variable, 'rb') as file:
        PPI = pickle.load(file)

## First time reading the file
else:
    print("The file does not exist yet. Processing... ")
    complete_interactome = pd.read_csv(interactome_file_name, delimiter='\t')

    ### Filter out non-human interactions
    human_interactome = complete_interactome[
        (complete_interactome["Organism ID Interactor A"] == human_ID)
        & (complete_interactome["Organism ID Interactor B"] == human_ID)]

    ### Filter out non-physical interactions
    human_physical_interactome = human_interactome[
        (human_interactome["Experimental System Type"] == interaction_of_interest)]

    ### Delete duplicates
    PPI_NoDuplicates = human_physical_interactome.drop_duplicates(
        subset=[first_gene_symbol_indicator, second_gene_symbol_indicator])

    ### Delete self-loops
    PPI_NoDuplicates_NoSelfLoop = PPI_NoDuplicates[
        PPI_NoDuplicates[first_gene_symbol_indicator]
        !=
        PPI_NoDuplicates[second_gene_symbol_indicator]]

    PPI = PPI_NoDuplicates_NoSelfLoop
    with open(PPI_list_variable, 'wb') as file:
        pickle.dump(PPI, file)

    del complete_interactome, human_interactome, human_physical_interactome, \
        PPI_NoDuplicates, PPI_NoDuplicates_NoSelfLoop

# Create the Interactome ---------------------------------------------------------------------------------------

graph_of_PPI = nx.from_pandas_edgelist(PPI, protein_A_identifier, protein_B_identifier)

## Find the Large Connected Component (LCC)

connected_components = nx.connected_components(graph_of_PPI)
largest_cc = max(connected_components, key=len)
LCC_subgraph = nx.subgraph(graph_of_PPI, largest_cc).copy()

## Plotting the LCC subgraph
plt.figure(figsize=(10, 8))
nx.draw(LCC_subgraph, with_labels=True, node_size=50, font_size=8)
plt.title("Largest Connected Component (LCC) of PPI Graph")
plt.show()

# Importing the seed genes ----------------------------------------------------------------------------------
'''
seed_genes_df = pd.read_csv(seeds_file_name, sep='\t')
seed_genes_symbols = seed_genes_df["Gene"]

disease_interactome = PPI[
    (PPI[first_gene_symbol_indicator].isin(seed_genes_symbols))
    &
    (PPI[second_gene_symbol_indicator].isin(seed_genes_symbols))
    ]

graph_of_disease_interactome = nx.from_pandas_edgelist(
    disease_interactome,
    first_gene_symbol_indicator, second_gene_symbol_indicator)

nx.draw(graph_of_disease_interactome)
'''
# there are some self-loops, we have to debug the initial list
# processing of the PPI; maybe the fact that this graph is determined by the
# official gene symbols could have changed something.
# Yes, definitely with this new version we don't have self loops no more
