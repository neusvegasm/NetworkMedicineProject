# -*- coding: utf-8 -*-

import pandas as pd
import networkx as nx
import os
import pickle
import matplotlib.pyplot as plt


Interactome_file_name = "BIOGRID-ALL-4.4.239.tab3.txt" # starter file
Interactome_file_name_corrected = "BIOGRID-ORGANISM-Homo_sapiens-4.4.239.tab3.txt" # read this
Interactome_file_name_corrected_updated = "BIOGRID-ORGANISM-Homo_sapiens-4.4.240.tab3.txt"
Seeds_file_name = "DISEASES_Summary_GDA_CURATED_C0025202.tsv"
human_ID = 9606
interaction_of_interest = "physical"
Protein_A_identifier = "Entrez Gene Interactor A"
Protein_B_identifier = "Entrez Gene Interactor B"
PPI_list_variable = "PPI_of_Interest.pk1"
PPI_list_variable_second_try = "PPI_of_interest_second_try.pk1"
PPI_list_variable_corrected = "PPI_of_interest_corrected" # read this
PPI_list_variable_corrected_updated = "PPI_of_interest_corrected updated"
First_gene_symbol_indicator = "Official Symbol Interactor A"
Second_gene_symbol_indicator = "Official Symbol Interactor B"

#%% Importing and Processing the Interactome

if os.path.exists(PPI_list_variable_corrected_updated):
    print("The file already exists! Importing information from the file ...")
    with open(PPI_list_variable_corrected_updated, 'rb') as file:
        PPI = pickle.load(file)
    
else:
    print("The file does not exist yet. Processing... ")
    complete_interactome = pd.read_csv(Interactome_file_name_corrected_updated, 
                                       delimiter='\t')
    
    human_interactome = complete_interactome[
        (complete_interactome["Organism ID Interactor A"] == human_ID)
        & (complete_interactome["Organism ID Interactor B"] == human_ID)]
    
    human_physical_interactome = human_interactome[
        (human_interactome["Experimental System Type"] == interaction_of_interest)]
    
    ## human_physical_interactome.to_csv("human_physical_interactome.csv", index = False)
    
    # Duplicate Removal
    
    # PPI_NoDuplicates = human_physical_interactome.drop_duplicates(
    #     subset = [Protein_A_identifier, Protein_B_identifier])
    
    # PPI_NoDuplicates_NoSelfLoop = PPI_NoDuplicates[
    #     PPI_NoDuplicates[Protein_A_identifier] != PPI_NoDuplicates[Protein_B_identifier]]
    
    PPI_NoDuplicates = human_physical_interactome.drop_duplicates(
        subset = [First_gene_symbol_indicator, Second_gene_symbol_indicator])
    
    PPI_NoDuplicates_NoSelfLoop = PPI_NoDuplicates[
        PPI_NoDuplicates[First_gene_symbol_indicator] 
        != 
        PPI_NoDuplicates[Second_gene_symbol_indicator]]
    # Note that with this second version, 
    # the number of edges are "slightly" different 
    # (9M with the first version vs 8.7 with the second version)
    # I think I have to review Database and GeneSymbols 
    
    PPI = PPI_NoDuplicates_NoSelfLoop
    with open(PPI_list_variable_corrected_updated, 'wb') as file:
        pickle.dump(PPI, file)
        
    del complete_interactome, human_interactome, human_physical_interactome, \
        PPI_NoDuplicates, PPI_NoDuplicates_NoSelfLoop



#%% Create the Interactome LCC
Graph_of_PPI = nx.from_pandas_edgelist(
    PPI, First_gene_symbol_indicator, Second_gene_symbol_indicator
    )
print("Features of the initial PPI:\n" + str(Graph_of_PPI))


## To verify
Connected_Components = nx.connected_components(Graph_of_PPI)
 
largest_cc= max(Connected_Components, key=len)

LCC_subgraph = nx.subgraph(Graph_of_PPI, largest_cc).copy()
print("Features of the LCC from the initial PPI:\n" + str(LCC_subgraph))



#%% Importing and the seed genes
# Note that this will not be the final version, because
Seed_genes_df = pd.read_csv(Seeds_file_name, sep = '\t')
Seed_genes_Symbols = Seed_genes_df["Gene"]
   
Disease_Interactome_linked_seeds = PPI[
    (PPI[First_gene_symbol_indicator].isin(Seed_genes_Symbols)) 
    &
    (PPI[Second_gene_symbol_indicator].isin(Seed_genes_Symbols))
    ] 
# with this I s of textract from the original PPI just the seeds that have some links
# with other seeds in the seeds list. That means the result will be already 
# the connected components of the disease interactome 

Seed_genes_with_links = list(
    set(Disease_Interactome_linked_seeds[First_gene_symbol_indicator].unique()).union(
        Disease_Interactome_linked_seeds[Second_gene_symbol_indicator].unique())
    )
Seed_genes_without_links = list(
    set(Seed_genes_with_links).symmetric_difference(Seed_genes_Symbols))
# I keep track of the seeds that are not linked with other in the original 
# interactome, because, even if they won't be in the LCC disease interactome
# they will be interesting in the validation procedure and in the inference 
# procedure as well (at least I think that)

Graph_of_connected_components_of_disease_interactome = nx.from_pandas_edgelist(
    Disease_Interactome_linked_seeds, 
    First_gene_symbol_indicator, Second_gene_symbol_indicator)
# in our case there is from the start just one connected component - that is 
# the largest one.
# Anyway to keep the code more general I explicitly isolated the LCC
Connected_components_of_disease_interactome = nx.connected_components(
    Graph_of_connected_components_of_disease_interactome)
largest_cc_of_disease_interactome = max(Connected_components_of_disease_interactome, 
                                        key = len)
Graph_of_LCC_of_disease_interactome = nx.subgraph(
    Graph_of_connected_components_of_disease_interactome, 
    largest_cc_of_disease_interactome).copy()


Complete_graph_of_disease_interactome = Graph_of_connected_components_of_disease_interactome.copy()
for i, missing_seed in enumerate(Seed_genes_without_links):
    Complete_graph_of_disease_interactome.add_node(i, label = missing_seed)


print("Features of the complete Disease Interactome:\n" 
      + str(Complete_graph_of_disease_interactome))

plt.show()
nx.draw_networkx(Complete_graph_of_disease_interactome, 
                 with_labels = False, node_size = 10, width = 0.5)

print("Features of the LCC from the Disease Interactome:\n" 
      + str(Graph_of_LCC_of_disease_interactome))

plt.show()
nx.draw_networkx(Graph_of_LCC_of_disease_interactome, 
                 with_labels = False, node_size = 10, width = 0.5)




#%% Central Measures of the Graph
Degrees_disease_interactome = dict(Graph_of_LCC_of_disease_interactome.degree())
type(Degrees_disease_interactome)
Degrees_disease_interactome_sorted = sorted(Degrees_disease_interactome.items(),
                                                 key=lambda x: x[1], reverse=True)
Seed_genes_Sorted_by_degree = [Generic_sorted_Seed[0] 
                               for Generic_sorted_Seed 
                               in Degrees_disease_interactome_sorted]

Betweenness_disease_interactome = nx.betweenness_centrality(
    Graph_of_LCC_of_disease_interactome)
Eig_Centrality_disease_interactome = nx.eigenvector_centrality(
    Graph_of_LCC_of_disease_interactome)
Collection_of_Centrality_Measures = [(Seed,
                                      Degrees_disease_interactome[Seed],
                                      Betweenness_disease_interactome[Seed],
                                      Eig_Centrality_disease_interactome[Seed],
                                      Betweenness_disease_interactome[Seed]/
                                      Eig_Centrality_disease_interactome[Seed]) 
                                      for Seed in Seed_genes_Sorted_by_degree]

Table_Visualization_of_the_Centrality_Measures = pd.DataFrame(
    Collection_of_Centrality_Measures, 
    columns = ['GeneSymbol','Degree', 'Betweenness', 'EigenVector Centrality', 
               'Ratio Between-Degree'])
print(Table_Visualization_of_the_Centrality_Measures)

# We can optimize all by writing a function or a class,
# that do all this stuff (if we want)

#%% Node Degree and Node Betweeness scatter plot
List_of_Degree_Values = [Centrality_Measure[1] 
                         for Centrality_Measure 
                         in Collection_of_Centrality_Measures]
List_of_Betweenness_Values = [Centrality_Measure[2] 
                         for Centrality_Measure 
                         in Collection_of_Centrality_Measures]

plt.show()
plt.scatter(List_of_Degree_Values, List_of_Betweenness_Values, 
            s = 10)
plt.xlabel('Deegres of the Seed Genes')
plt.ylabel('Betweenness of the Seed Genes')
plt.title("Deegrees vs Betwennesses Centralities")