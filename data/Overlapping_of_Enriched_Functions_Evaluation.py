'''Evaluation of the Overlapping of Functional Modules*'''

#%% Import Libraries

import pandas as pd

#%% Variables

# Name of the files to load 
Seed_GO_Molecular_Function_file = 'Seed_GO_Molecular_Function_2023_table.txt'
Inferred_GO_Molecular_Function_file = 'Inferred_GO_Molecular_Function_2023_table.txt'
#
Seed_GO_Cellular_Component_file = "Seed_GO_Cellular_Component_2023_table.txt"
Inferred_GO_Cellular_Component_file = 'Inferred_GO_Cellular_Component_2023_table.txt'
#
Seed_GO_Biological_Process_file = 'Seeds_GO_Biological_Process_2023_table.txt'
Inferred_GO_Biological_Process_file = 'Inferred_GO_Biological_Process_2023_table.txt'
#
Seed_KEGG_file = 'Seed_KEGG_2021_Human_table.txt'
Inferred_KEGG_file = 'Inferred_KEGG_2021_Human_table.txt'
#
Seed_REACTOME_file = 'Seed_Reactome_Pathways_2024_table.txt'
Inferred_REACTOME_file = 'Inferred_Reactome_Pathways_2024_table.txt'


# Other variables defined
Tuple_of_GO_Molecular_Function_file = (Seed_GO_Molecular_Function_file, 
                                       Inferred_GO_Molecular_Function_file)
Tuple_of_GO_Cellular_Component_file = (Seed_GO_Cellular_Component_file, 
                                       Inferred_GO_Cellular_Component_file)
Tuple_of_GO_Biological_Process_file = (Seed_GO_Biological_Process_file,
                                       Inferred_GO_Biological_Process_file)
Tuple_of_KEGG_file = (Seed_KEGG_file, Inferred_KEGG_file)
Tuple_of_REACTOME_file = (Seed_REACTOME_file, Inferred_REACTOME_file)

List_of_Tuples_of_File_Name = [Tuple_of_GO_Molecular_Function_file, 
                               Tuple_of_GO_Cellular_Component_file,
                               Tuple_of_GO_Biological_Process_file,
                               Tuple_of_KEGG_file, Tuple_of_REACTOME_file]

List_of_Nested_Dictionaries_Names = \
    ['1. Dictionary of GO Molecular Function', '2. Dictionary of GO Cellular Component',
     '3. Dictionary of GO Biological Process', '4. Dictionary of KEGG', '5. Dictionary of REACTOME']
Columns_of_Table_Visualization = ['Enrichement of the Seeds',
                                  'Enrichment of the Inferred Genes',
                                  'Overlapping of the sets', 'Evaluated Metric']
Rows = ['GO Molecular Function', 'GO Cellular Component', 'GO Biological Process',
        'KEGG', 'REACTOME']


#%% Functions to use

def Overlapping_Evaluation(File_Name_Tuple):
    Seed_Enriched_function_file_name = File_Name_Tuple[0]
    Inferred_Enriched_function_file_name = File_Name_Tuple[1]
    Seed_Enriched_function = pd.read_csv(Seed_Enriched_function_file_name, 
                                             delimiter ='\t')
    Seed_Enriched_function_Term = set(Seed_Enriched_function.iloc[:, 0])
    Inferred_Enriched_function = pd.read_csv(Inferred_Enriched_function_file_name, 
                                             delimiter ='\t')
    Inferred_Enriched_function_Term = set(Inferred_Enriched_function.iloc[:, 0])
    Overlapping_Enriched_Function_Term = Seed_Enriched_function_Term.intersection(
             Inferred_Enriched_function_Term)

    Overlapping_Metric = len(Overlapping_Enriched_Function_Term) / \
        min(len(Seed_Enriched_function_Term), len(Inferred_Enriched_function_Term))
    Dictionary_of_Enriched_Functions =\
        {'Enrichment of the Seeds' : Seed_Enriched_function_Term, 
         'Enrichment of the Inferred Genes' : Inferred_Enriched_function_Term, 
         'Overlapping of the sets' : Overlapping_Enriched_Function_Term,
         'Evaluated Metric' : Overlapping_Metric}    
    return Dictionary_of_Enriched_Functions


#%% Overlapping Computation

Dictionary_of_Enriched_Functions = {}
Table_Visualization_of_the_Enrichment_Results = pd.DataFrame(
    columns = Columns_of_Table_Visualization)

for it_on_Enriched_Pathway in range(len(List_of_Nested_Dictionaries_Names)):
    Nested_Dictionary = Overlapping_Evaluation(
            List_of_Tuples_of_File_Name[it_on_Enriched_Pathway])
    Dictionary_of_Enriched_Functions[List_of_Nested_Dictionaries_Names[it_on_Enriched_Pathway]] = \
        Nested_Dictionary
    Data_of_The_Row_of_Table_Visualization = \
        [[len(Nested_Dictionary['Enrichment of the Seeds']),
          len(Nested_Dictionary['Enrichment of the Inferred Genes']),
          len(Nested_Dictionary['Overlapping of the sets']),
          Nested_Dictionary['Evaluated Metric']]]
    Update_of_Table_Visualization =  pd.DataFrame(Data_of_The_Row_of_Table_Visualization, 
                                                  columns = Columns_of_Table_Visualization)
    Table_Visualization_of_the_Enrichment_Results = pd.concat(
        [Table_Visualization_of_the_Enrichment_Results, Update_of_Table_Visualization],
        ignore_index=True)

Table_Visualization_of_the_Enrichment_Results.index = Rows
print(Table_Visualization_of_the_Enrichment_Results)