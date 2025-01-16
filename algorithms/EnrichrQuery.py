
import requests
import pandas as pd
import warnings
import time

#%%
Columns_Names = ['Term', 'Overlap', 'P-value', 'Adjusted P-value',
                  'Old P-value', 'Old Adjusted P-value', 'Odds Ratio',
                  'Combined Score', 'Genes']
warnings.simplefilter(action='ignore', category=FutureWarning)


#%% Function for sending genes to Enrichr
def enrichr_query(gene_list, description):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(gene_list)
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    return response.json()

def fetch_enrichment_results(url):
    for attempt in range(10):  # Prova fino a 5 volte
        try:
            response = requests.get(url)
            response.raise_for_status()  # Solleva un'eccezione per errori HTTP
            return response
        except requests.exceptions.RequestException as e:
            print(f"Errore: {e}. Tentativo {attempt + 1} di 5.")
            time.sleep(2)  # Attendi 2 secondi prima di riprovare
    raise Exception("Errore nel recupero dei risultati di arricchimento dopo 5 tentativi.")

# Function for Produce the results of the enrichment
# def get_enrichment_results(user_list_id, ontology):
#     ENRICHR_URL = \
#         f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={ontology}'
#     response = requests.get(ENRICHR_URL)
#     if not response.ok:
#         raise Exception('Error fetching enrichment results')
#     return response.json()

def get_enrichment_results(user_list_id, ontology):
    ENRICHR_URL = \
        f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={ontology}'
    response = fetch_enrichment_results(ENRICHR_URL)
    if not response.ok:
        raise Exception('Error fetching enrichment results')
    return response.json()

#%% Lista di geni di esempio
# gene_list = ['RELA']
# description = 'Example gene list'

# # Inviare la lista di geni a Enrichr
# response = enrichr_query(gene_list, description)
# user_list_id = response['userListId']

# # Ottenere i risultati per le tre ontologie
# ontologies = ['GO_Biological_Process_2023',
#               'GO_Molecular_Function_2023', 
#               'GO_Cellular_Component_2023']
# results = {}
# for ontology in ontologies:
#     results[ontology] = get_enrichment_results(user_list_id, ontology)

# # Visualizzare i risultati
# Columns_Names = ['Term', 'Overlap', 'P-value', 'Adjusted P-value',
#                  'Old P-value', 'Old Adjusted P-value', 'Odds Ratio',
#                  'Combined Score', 'Genes']

# GeneList_Enriched_Onthology = pd.DataFrame(columns = Columns_Names)
# for ontology, result in results.items():
#     print(f'\nResults for {ontology}:')
#     ith_GeneList_Onthology = pd.DataFrame(result[ontology], 
#                       columns = Columns_Names)
#     print(ith_GeneList_Onthology.head())
#     GeneList_Enriched_Onthology = pd.concat(
#         [GeneList_Enriched_Onthology, ith_GeneList_Onthology],
#         ignore_index = True)

# Set_of_GeneList_Enriched_Onthologies = set(GeneList_Enriched_Onthology.iloc[:, 1])

# for ontology, result in results.items():
#     print(f'\nResults for {ontology}:')
#     df = pd.DataFrame(result[ontology])
#     print(df.head())

#%%

# Seed_GO_Molecular_Function_file = 'Seed_GO_Molecular_Function_2023_table.txt'
# Seed_GO_Cellular_Component_file = "Seed_GO_Cellular_Component_2023_table.txt"
# Seed_GO_Biological_Process_file = 'Seeds_GO_Biological_Process_2023_table.txt'

# Seed_Onthologies_List = [Seed_GO_Molecular_Function_file, Seed_GO_Cellular_Component_file,
#                          Seed_GO_Biological_Process_file]

# Seed_Enriched_Onthologies = pd.DataFrame(columns = Columns_Names)
# for it_on_Onthologies_List in range(len(Seed_Onthologies_List)):
#     ith_Seed_Onthology = pd.read_csv(Seed_Onthologies_List[it_on_Onthologies_List],
#                                        delimiter ='\t')
#     Seed_Enriched_Onthologies =  pd.concat([Seed_Enriched_Onthologies, ith_Seed_Onthology],
#                                            ignore_index = True)
    
# Set_of_Seed_Enriched_Onthologies = set(Seed_Enriched_Onthologies.iloc[:, 0])

# with open('Set_of_Seed_Enriched_Onthologies', 'wb') as file:
#     pickle.dump(Set_of_Seed_Enriched_Onthologies, file)

# with open('Set_of_Seed_Enriched_Onthologies', 'rb') as file:
#   Set_of_Seed_Enriched_Onthologies = pickle.load(file)

def Enrich_Main_List(gene_list):
    description = 'Example gene list'
    response = enrichr_query(gene_list, description)
    user_list_id = response['userListId']
    ontologies = ['GO_Biological_Process_2023',
                  'GO_Molecular_Function_2023', 
                  'GO_Cellular_Component_2023']
    results = {}
    for ontology in ontologies:
        results[ontology] = get_enrichment_results(user_list_id, ontology)
    
    GeneList_Enriched_Onthology = pd.DataFrame(columns = Columns_Names)
    for ontology, result in results.items():
        # print(f'\nResults for {ontology}:')
        ith_GeneList_Onthology = pd.DataFrame(result[ontology], 
                          columns = Columns_Names)
        # print(ith_GeneList_Onthology.head())
        GeneList_Enriched_Onthology = pd.concat(
            [GeneList_Enriched_Onthology, ith_GeneList_Onthology],
            ignore_index = True)

    Set_of_GeneList_Enriched_Onthologies = set(GeneList_Enriched_Onthology.iloc[:, 1])
    
    return Set_of_GeneList_Enriched_Onthologies

############################################################################

def Compute_Enrichment_Weigts(gene_list, Set_of_Seed_Enriched_Functions, 
                              Type_of_Normalization = 'Seed'):
    description = 'Example gene list'

    # Inviare la lista di geni a Enrichr
    response = enrichr_query(gene_list, description)
    user_list_id = response['userListId']

    # Ottenere i risultati per le tre ontologie
    ontologies = ['GO_Biological_Process_2023',
                  'GO_Molecular_Function_2023', 
                  'GO_Cellular_Component_2023']
    results = {}
    for ontology in ontologies:
        results[ontology] = get_enrichment_results(user_list_id, ontology)



    GeneList_Enriched_Onthology = pd.DataFrame(columns = Columns_Names)
    for ontology, result in results.items():
        # print(f'\nResults for {ontology}:')
        ith_GeneList_Onthology = pd.DataFrame(result[ontology], 
                          columns = Columns_Names)
        # print(ith_GeneList_Onthology.head())
        GeneList_Enriched_Onthology = pd.concat(
            [GeneList_Enriched_Onthology, ith_GeneList_Onthology],
            ignore_index = True)

    Set_of_GeneList_Enriched_Onthologies = set(GeneList_Enriched_Onthology.iloc[:, 1])
    
    if Type_of_Normalization == 'Seed':
        Weight = len(Set_of_Seed_Enriched_Functions.intersection(Set_of_GeneList_Enriched_Onthologies))\
            / len(Set_of_Seed_Enriched_Functions)
    elif Type_of_Normalization == 'Self':
        Weight = len(Set_of_Seed_Enriched_Functions.intersection(Set_of_GeneList_Enriched_Onthologies))\
            / (len(Set_of_GeneList_Enriched_Onthologies)+1)
        # the +1 is for avoid ZeroDivision
    else:
        print("Error in Type_of_Normalization! You have to choose between 'Seed' and 'Self'\n")
    
    return Weight
    

