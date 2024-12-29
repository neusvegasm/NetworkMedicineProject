##%% Importing the libraries

import pickle
import pandas as pd
import random
import numpy as np
from DIAMOnD import DIAMOnD
from DiaBLE import DiaBLE
from DIFFUSION import run_heat_diffusion
# import networkx as nx


#%% Define all the parameters and strings of the script

Graph_of_PPI_file_name = "Graph_of_PPI.pickle"
Graph_of_Disease_Interactome_file_name = "Graph_of_Disease_Interactome.pickle"
Graph_of_LCC_of_disease_interactome_file_name = "Graph_of_LCC_of_disease_interactome.pickle"
Seed_genes_Symbol_file_name = "Seed_genes_Symbol.pickle"
Genes_to_ignore_file = "Genes_ignored_by_DIAMOnD.txt"
PPI_list_variable_corrected_updated = "PPI_of_interest_corrected updated"
node_inferred_by_DIAMOnD_file = 'Node Inferred by DIAMOnD.pickle'
node_inferred_by_DiaBLE_file = 'Node Inferred by DiaBLE'
First_gene_symbol_indicator = "Official Symbol Interactor A"
Second_gene_symbol_indicator = "Official Symbol Interactor B"
File_Name_of_Metrics_of_DiaBLE_inference = "Metrics_of_DiaBLE_inference"
File_Name_of_Metrics_of_DIAMOnD_inference = "Metrics_of_DIAMOnD_inference"
File_Name_of_Metrics_of_Diffusion_inference = "Metrics_of_Diffusion_inference"



k_parameter_of_splitting = 5;
Top_positions_reduced = [25, 50]
Top_positions = [25, 50, 65, 130, 252] 
Diffusion_times_to_try = [0.002, 0.005, 0.01]


Precision_Matrix = [[None for _ in range(k_parameter_of_splitting)]
                    for _ in range(len(Top_positions))]
Recall_Matrix = [[None for _ in range(k_parameter_of_splitting)]
                    for _ in range(len(Top_positions))]
Founded_True_Positive = [[None for _ in range(k_parameter_of_splitting)]
                    for _ in range(len(Top_positions))]
F1_score_Matrix = [[None for _ in range(k_parameter_of_splitting)]
                    for _ in range(len(Top_positions))]
Metrics = [Founded_True_Positive, Precision_Matrix, 
           Recall_Matrix, F1_score_Matrix]



#%% Loading all the data from the directory

with open(Graph_of_PPI_file_name, "rb") as file:
    Graph_of_PPI = pickle.load(file)
with open(Graph_of_Disease_Interactome_file_name, "rb") as file:
    Graph_of_Disease_Interactome = pickle.load(file)
with open(Graph_of_LCC_of_disease_interactome_file_name, "rb") as file:
    Graph_of_LCC_of_disease_interactome = pickle.load(file)
with open(PPI_list_variable_corrected_updated, 'rb') as file:
    PPI = pickle.load(file)
    
Seed_genes_Symbols = pd.read_pickle(Seed_genes_Symbol_file_name)
with open(Genes_to_ignore_file, "r") as file:
    Genes_to_ignore = file.readlines()
Genes_to_ignore = [Gene.strip() for Gene in Genes_to_ignore]
Seed_genes_Symbols = Seed_genes_Symbols[~Seed_genes_Symbols.isin(Genes_to_ignore)]


#%% Defining the function that we will use
# (it seems that all is working, anyway try to do a Debug to verify)

def randomly_split_of_the_list(list_of_gene, n_split):
    random.shuffle(list_of_gene)
    List_of_subarray = np.array_split(list_of_gene, n_split)
    List_of_sublists = [None for _ in range(n_split)]
    for i in range(n_split):
        List_of_sublists[i] = [str(element) 
                               for element
                               in list(List_of_subarray[i])]
    return List_of_sublists

def Calculating_Metrics(nodes_added, Validation_list_of_seeds,
                       it_on_Top_position, it_on_validation):
    True_positive = len(set(nodes_added).intersection(Validation_list_of_seeds))
    Founded_True_Positive[it_on_Top_position][it_on_validation] = True_positive
    Precision_Matrix[it_on_Top_position][it_on_validation] = \
        True_positive / len(nodes_added)
    Recall_Matrix[it_on_Top_position][it_on_validation] = \
        True_positive / len(Validation_list_of_seeds)
    if Precision_Matrix[it_on_Top_position][it_on_validation]\
        + Recall_Matrix[it_on_Top_position][it_on_validation]:
            F1_score_Matrix[it_on_Top_position][it_on_validation] = \
                2 * (Precision_Matrix[it_on_Top_position][it_on_validation]\
                     * Recall_Matrix[it_on_Top_position][it_on_validation]) \
                    / (Precision_Matrix[it_on_Top_position][it_on_validation]\
                       + Recall_Matrix[it_on_Top_position][it_on_validation])  
    else:
        F1_score_Matrix[it_on_Top_position][it_on_validation] = None
                    
    Metrics = [Founded_True_Positive, Precision_Matrix, 
               Recall_Matrix, F1_score_Matrix]
    return Metrics
    

def Implementing_Inference_Algorithms(Choosed_algorithm,
        Graph_of_PPI, Seed_genes_Symbols, Top_positions,
        outfile = node_inferred_by_DIAMOnD_file, Diffusion_time = 0.005):
    
    List_of_sublists_of_seeds = randomly_split_of_the_list(
        list(Seed_genes_Symbols), k_parameter_of_splitting)
    for it_on_validation in range(k_parameter_of_splitting):
        List_of_sublists_of_seeds_copy = List_of_sublists_of_seeds.copy()
        Validation_list_of_seeds = List_of_sublists_of_seeds[it_on_validation]
        List_of_sublists_of_seeds_copy.remove(List_of_sublists_of_seeds[it_on_validation])
        Training_list_of_seeds = [Gene 
                                  for sublist in List_of_sublists_of_seeds_copy
                                  for Gene in sublist]
        match Choosed_algorithm:
            case 'DIAMOnD':
                for it_on_Top_position in range(len(Top_positions)):
                    list_of_nodes_added = DIAMOnD(Graph_of_PPI, Training_list_of_seeds,
                                                  Top_positions[it_on_Top_position],
                                                  alpha = 1, 
                                                  outfile = node_inferred_by_DIAMOnD_file)
                    nodes_added = [Generic_tuple[0] for Generic_tuple 
                                   in list_of_nodes_added]                    
                    Metrics = Calculating_Metrics(nodes_added,
                                                  Validation_list_of_seeds,
                                                  it_on_Top_position, it_on_validation)
                
            case 'DiaBLE':
                for it_on_Top_position in range(len(Top_positions)):
                    list_of_nodes_added = DiaBLE(Graph_of_PPI, Training_list_of_seeds,
                                                  Top_positions[it_on_Top_position],
                                                  alpha = 1, 
                                                  outfile = node_inferred_by_DIAMOnD_file)
                    nodes_added = [Generic_tuple[0] for Generic_tuple 
                                   in list_of_nodes_added]
                    Metrics = Calculating_Metrics(nodes_added,
                                                  Validation_list_of_seeds,
                                                  it_on_Top_position, it_on_validation)
            
            case 'Diffusion':
                for it_on_Top_position in range(len(Top_positions)):
                    nodes_added = run_heat_diffusion(
                        Graph_of_PPI, Training_list_of_seeds, diffusion_time = Diffusion_time,
                        n_positions = Top_positions[it_on_Top_position])
                    Metrics = Calculating_Metrics(nodes_added,
                                                  Validation_list_of_seeds,
                                                  it_on_Top_position, it_on_validation)
            case _:
                print("Invalid Inference Method. \n" +
                      "Please choose one from the following:\n" +
                      "'DIAMOND', 'DiaBLE', 'Diffusion' \n")
                Metrics = []
    return Metrics

def Handle_ZeroDivisionError(array):
    means = []
    stds = []
    for matrix in array:
        matrix_means = []
        matrix_stds = []
        for row in matrix:
            try:
                mean = np.nanmean(row)
                std = np.nanstd(row)
            except ZeroDivisionError:
                mean = np.nan
                std = np.nan
            matrix_means.append(mean)
            matrix_stds.append(std)
        means.append(matrix_means)
        stds.append(matrix_stds)
    return means, stds
            
def Table_Visualization_of_Metrics(Metrics):
    Metrics_array = np.array(Metrics)
    Metrics_array = Metrics_array.astype(np.float64)
    try:
        Mean_of_the_Metrics = \
            np.nanmean(Metrics_array, axis = 2)
        Standard_Deviation_of_the_Metrics = \
            np.nanstd(Metrics_array, axis = 2)
    except ZeroDivisionError:
        Mean_of_the_Metrics, Standard_Deviation_of_the_Metrics = \
            Handle_ZeroDivisionError(Metrics_array)
    Mean_of_the_Metrics = \
        np.round(Mean_of_the_Metrics, 3)        
    Standard_Deviation_of_the_Metrics = \
        np.round(Standard_Deviation_of_the_Metrics, 3)

    # List_of_Tuple_Mean_and_Std =\
    #     [
    #      [(Mean_of_the_Metrics[i, j], 
    #       Standard_Deviation_of_the_Metrics[i, j]) 
    #       for j in range(Standard_Deviation_of_the_Metrics.shape[1])]
    #       for i in range(Mean_of_the_Metrics.shape[0])
    #     ]
        
    List_of_PlusMinus_Mean_and_Std = \
        [ 
          [str(Mean_of_the_Metrics[i, j]) + ' ± ' + 
           str(Standard_Deviation_of_the_Metrics[i, j])
           for j in range(Standard_Deviation_of_the_Metrics.shape[1])
           ]
           for i in range(Mean_of_the_Metrics.shape[0])
        ]
    
    Columns_of_the_Table = ["Top 25", "Top 50", "Top 65", "Top 130", "Top 252"]
    Rows_of_the_Table = ["Precision", "Recall", "F1 Score"]
    Table_Visualization = pd.DataFrame(List_of_PlusMinus_Mean_and_Std,
                                       columns = Columns_of_the_Table,
                                       index =  Rows_of_the_Table)
    return Table_Visualization


#%% Applying the Algorithms
Diffusion_times_to_try = [0.002, 0.005, 0.01]

Diffusion_Metrics_on_various_times = \
    [[[None for _ in range(k_parameter_of_splitting)]
     for _ in range(len(Top_positions))] for _ in range(len(Diffusion_times_to_try))]


# # Creation of the group of Matrices of Metrics
# Choosed_algorithm = 'Diffusion'
# for it_on_diffusion_times in range(len(Diffusion_times_to_try)):
#     Diffusion_Metrics = \
#         Implementing_Inference_Algorithms(Choosed_algorithm, Graph_of_PPI, 
#                                           Seed_genes_Symbols, Top_positions, 
#                                           Diffusion_time = \
#                                               Diffusion_times_to_try[it_on_diffusion_times])
#     Diffusion_Metrics_on_various_times[it_on_diffusion_times] = Diffusion_Metrics

# # Saving of the group of Matrices of Matrix
# with open(File_Name_of_Metrics_of_Diffusion_inference, 'wb') as file:
#     pickle.dump(Diffusion_Metrics_on_various_times, file)
    
# Downloading of the group of Matrices of Metrics
with open(File_Name_of_Metrics_of_Diffusion_inference, 'rb') as file:
    Diffusion_Metrics_on_various_times_download = pickle.load(file)
    
Table_Visualization_of_Diffusion_Metrics_iterated = \
    [None for _ in range(len(Diffusion_times_to_try))]

for it_on_Visualizations in range(len(Diffusion_times_to_try)):
    Table_Visualization_of_Diffusion_Metrics_iterated[it_on_Visualizations] = \
        Table_Visualization_of_Metrics(
            Diffusion_Metrics_on_various_times_download[it_on_Visualizations][1:])
    print(Table_Visualization_of_Diffusion_Metrics_iterated[it_on_Visualizations])
del it_on_Visualizations
    


#%% Implementing DiaBLE

# # Creation of the Matrices of the Metrics
# Choosed_algorithm = 'DiaBLE'
# File_Name_of_Metrics_of_DiaBLE_inference = "Metrics_of_DiaBLE_inference"

# DiaBLE_Metrics = Implementing_Inference_Algorithms(Choosed_algorithm, Graph_of_PPI, 
#                                   Seed_genes_Symbols, Top_positions)

# # Saving the Matrices of the Metrics
# with open(File_Name_of_Metrics_of_DiaBLE_inference, 'wb') as file:
#     pickle.dump(DiaBLE_Metrics, file)
    
# Downloading the Matrices of the Metrics preaviously saved
with open(File_Name_of_Metrics_of_DiaBLE_inference, 'rb') as file:
    DiaBLE_metrics_download = pickle.load(file)
    
Table_Visualization_of_DiaBLE = \
    Table_Visualization_of_Metrics(DiaBLE_metrics_download[1:])
    
#%% Implementing DIAMOnD

# # Creation of the Matrices of the Metrics
# Choosed_algorithm = 'DIAMOnD'
# File_Name_of_Metrics_of_DIAMOnD_inference = "Metrics_of_DiaBLE_inference"

# DIAMOnD_Metrics = Implementing_Inference_Algorithms(Choosed_algorithm, Graph_of_PPI, 
#                                   Seed_genes_Symbols, Top_positions)

# # Saving the Matrices of the Metrics
# with open(File_Name_of_Metrics_of_DiaBLE_inference, 'wb') as file:
#     pickle.dump(DIAMOnD_Metrics, file)
    
# Downloading the Matrices of the Metrics preaviously saved
with open(File_Name_of_Metrics_of_DIAMOnD_inference, 'rb') as file:
    DIAMOnD_metrics_download = pickle.load(file)
    
Table_Visualization_of_DIAMOnD = \
    Table_Visualization_of_Metrics(DIAMOnD_metrics_download[1:])

    
