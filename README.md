# Gene Inference and Drug Repurposing Framework

This repository provides a comprehensive framework for inferring putative disease genes and evaluating drug repurposing opportunities using protein-protein interaction (PPI) networks and various algorithms.

## Project Structure

- **algorithms/**  
  Contains the implementation of algorithms used for gene inference.

- **data/**  
  Includes essential datasets for preprocessing, gene-disease association (GDA) analysis, overlapping evaluations, and drug repurposing. Key files:
  - **AliasGenes.xlsx**, **GenesToRemove.xlsx**, **GenesToReplace.xlsx**: Files for GDA analysis.
  - **Inferred_GO_Biological_Process_2023_table.txt**, **Inferred_GO_Cellular_Component_2023_table.txt**, **Inferred_GO_Molecular_Function_2023_table.txt**: Functional enrichment files for overlapping evaluations.
  - **Inferred_KEGG_2021_Human_table.txt**, **Inferred_Reactome_Pathways_2024_table.txt**: Pathway enrichment files for overlapping evaluations.
  - **Seed_GO_*_table.txt**, **Seed_KEGG_2021_Human_table.txt**, **Seed_Reactome_Pathways_2024_table.txt**: Seed files for overlapping evaluations.
  - **Overlapping_of_Enriched_Functions_Evaluation.py**: Script for functional overlap evaluation.
  - **gene_interaction_results-04_01_2025.tsv**: Gene interaction results.

- **results/**  
  Stores the results of gene inference and evaluation.

- **1-ListExtractor.ipynb**  
  Preprocesses input data and extracts relevant lists for constructing the PPI graph.

- **2-PPI_GraphFeatures.ipynb**  
  Builds the PPI graph and computes its features.

- **3-GDA.ipynb**  
  Performs Gene-Disease Association (GDA) analysis.

- **4-Inference.ipynb**  
  Evaluates and compares the performance of different inference algorithms.

- **5-GeneInference.ipynb**  
  Predicts genes using the best-performing algorithm.

- **6-OverlappingEvaluation.ipynb**  
  Assesses overlaps in gene predictions across methods.

- **7-DrugReporpusing.ipynb**  
  Conducts drug repurposing analysis using the DGIdb dataset.

- **requirements.txt**  
  Lists the Python dependencies required for this project.

- **.gitattributes**  
  Specifies settings for version control behavior.

- **.gitignore**  
  Excludes unnecessary files from version control.

## Getting Started

### Prerequisites
- Python 3.7 or higher
- Required Python libraries listed in `requirements.txt`

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/neusvegasm/NetworkMedicineProject.git
   cd NetworkMedicineProject
   ```
2. Install dependencies
   ```bash
   pip install -r requirements.txt
   ```
### Data Preparation
1. **Download PPI Data**:
   - Obtain the latest BioGRID release from [BioGRID](https://thebiogrid.org/).
   - Use the `organisms` tab3 file and extract data for **Homo sapiens** only.
   - Save the file in the `data/` directory.

2. **Gather Gene-Disease Associations (GDAs)**:
   - Use the file `DISEASES_Summary_GDA_CURATED_C0025202.tsv` from DisGeNet.

3. **Preprocess Data**:
   - Run `1-ListExtractor.ipynb` to prepare data for constructing the PPI graph.

### Workflow
1. **Build PPI Graph**: Execute `2-PPI_GraphFeatures.ipynb` to construct the graph and compute features.
2. **Analyze GDAs**: Use `3-GDA.ipynb` for Gene-Disease Association analysis.
3. **Algorithm Comparison**: Evaluate different algorithms using `4-Inference.ipynb`.
4. **Gene Prediction**: Infer putative genes with `5-GeneInference.ipynb`.
5. **Evaluate Overlaps**: Run `6-OverlappingEvaluation.ipynb` for functional overlap analysis.
6. **Drug Repurposing**: Perform drug repurposing using `7-DrugReporpusing.ipynb`.

## Results
The results of gene inference and drug repurposing are stored in the `results/` directory.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments
This project utilizes data from BioGRID, DGIdb, and other publicly available resources. Contributions to the algorithms and methodology are greatly appreciated.


