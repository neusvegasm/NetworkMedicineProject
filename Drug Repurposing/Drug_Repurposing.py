# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:11:35 2025

@author: atara
"""

import pandas as pd
# import requests
from pytrials.client import ClinicalTrials

#%% Function to check the clinical trials

target_fields = ['NCT Number', 'Study Title', 'Brief Summary', 'Study Results']

def clinical_trials_query(disease, drug):
    search_terms = f'{disease} AND {drug}'
    ct = ClinicalTrials()
    clinical_trials = ct.get_study_fields(search_expr = search_terms,
                                          fields = target_fields,
                                          max_studies = 1000, 
                                          fmt = "csv")
    column_names = clinical_trials[0]
    rows = clinical_trials[1:]
    Table_Visualization_of_clinical_trials = pd.DataFrame(rows, columns = column_names)
    
    keys = clinical_trials[0]
    keys.insert(0, 'Drug')
    clinical_trials = clinical_trials[1:]
    list_of_dictionaries = []
    for inner_list in clinical_trials:
        inner_list.insert(0, drug)
        dict_item = {keys[i]: inner_list[i] for i in range(len(keys))}
        list_of_dictionaries.append(dict_item)

    return list_of_dictionaries, Table_Visualization_of_clinical_trials


#%% DBIdb Evaluation

File_of_DBIdb = 'gene_interaction_results-04_01_2025.tsv'
Table_Visualization_of_DBIdb = pd.read_csv(File_of_DBIdb, sep = '\t')
Table_Visualization_of_DBIdb_approved = \
    Table_Visualization_of_DBIdb[Table_Visualization_of_DBIdb.iloc[:,2] == 'Approved']
    
Drugs_Identified_with_repetitions = Table_Visualization_of_DBIdb_approved['drug']
# ans = Drugs_Identified.unique()
Drugs_Identified = Drugs_Identified_with_repetitions.value_counts()
New_Order_of_the_Table = Drugs_Identified.index

Table_Visualization_of_DBIdb_approved_reindexed = \
    Table_Visualization_of_DBIdb_approved.set_index('drug').loc[
        New_Order_of_the_Table].reset_index()

Most_Important_Drugs_isolation = Drugs_Identified[Drugs_Identified == Drugs_Identified.max()]
Most_Important_Drugs_list =list(Most_Important_Drugs_isolation.index)


#%%
disease = "melanoma"

clinical_trials_collection = []
counter = 0
for drug in Most_Important_Drugs_list:
    clinical_trial_Dictionary_associated_to_drug, _ = clinical_trials_query(disease, drug)
    clinical_trials_collection.append(clinical_trial_Dictionary_associated_to_drug)

    # clinical_trials_collection.append(clinical_trials_query(disease, drug))
    counter += 1


