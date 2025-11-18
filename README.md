## Evaluations

- **SVM Model Evaluation:** Cross-validations and benchmark evaluation results and optimization details are available in the `SVM_optimization` folder.
- **Von Heijne Model Evaluation:** Cross-validations, and benchmark evaluation results are available in the `von_heijne_method` folder.

# Signal Peptide Prediction
This repository contains a complete pipeline for **prediction of signal peptides** using sequence based features and Support Vector Machines (SVMs) & Von Heime.

The project consists of **five major stages**:

1. **Data Collection** — UniProt retrieval, filtering, redundancy removal (MMseq2), data splitting
2. **Data Visualization** — Exploratory dataset analysis and biological insight
3. **Von Heime Method** - Firstly Von Heijne method is applied to predict signal peptides.
4. **Feature Extraction, Model Training (SVM) & Optimization** — Numerical representation of sequences, cross-validation, feature selection, SVM training
5. **Benchmark Evaluation & Error Analysis** — Independent testing + biological interpretation

Each stage is documented inside its respective directory, while this README provides a high-level overview.

## Summary of Data Collection Step

The **data_collection/** directory contains all scripts used for building the curated dataset used throughout the project.

Below is a short, script-based summary of the entire workflow.

### 1. Retrieve Raw Positive & Negative Datasets (from UniProt)
- Queries were executed on **17.09.2025**.
- Positive and negative datasets were downloaded using:

  - `data_collection_positive.py`
  - `data_collection_negative.py`

These scripts apply the initial UniProt filters (reviewed, non-fragment, length ≥ 40 aa, SP evidence / no-SP evidence).


### 2. Apply Additional Filtering
Filtering steps include:
- Removing signal peptides shorter than 14 aa
- Filtering out unwanted protein annotations / descriptions

This results in:

- `positive.tsv`, `positive.fasta`
- `negative.tsv`, `negative.fasta`

### 3. Redundancy Reduction via MMseqs2
Clustering was performed to remove redundancy and prevent bias:

- Identity ≥ 30%
- Coverage ≥ 40%

Outputs include cluster tables and FASTA files:
But the important files are these two which contain just representative sequences.
- `pos-cluster-results_rep_seq.fasta`
- `neg-cluster-results_rep_seq.fasta`

### 4. Generate Annotated Non-Redundant Datasets
Two scripts process cluster representatives:

- `filter_reps.py`
  Produces `positive_NR.tsv` and `negative_NR.tsv` with desired UniProt annotations.

- `fasta_sep.py`
  Creates one FASTA file per representative into `data_collection/fasta/`

### 5. Shuffle and Split into Train/Benchmark Sets
Performed using:

- `shuffle_NR.py`

Outputs:

- `train-pos.tsv`, `train-neg.tsv`
- `benchmark-pos.tsv`, `benchmark-neg.tsv`

Train: 80%
Benchmark: 20%


### 6. Build 5-Fold Cross-Validation Sets and Merge them into one whole file
Class-balanced splitting performed by:

- `split_cross_validation.py`
Now we have 5 different cross validation fold.

- `merge_data.py`
We combined all folds for training data and also benchmark data into one single file which contains all information together (for future analysis). But this time we have more columns for `Signal_Peptide`, `Transmembrane_Helix`,`Train`,`Fold`
With using these columns, we can create and test model. In addition, whole `.fasta` files are located into `/fasta` folder. The model go there and read the file.


## Data Visualization Overview

The `data_visualization/` directory contains all scripts used to analyze and visualize the biological and statistical properties of the dataset.
These visualizations help verify data quality, differences between positive and negative data and support biological interpretation before model training.

Note: We did all the same analysis both train set and benchmark set!

### 1. Shuffle Combined Dataset
First of all before analysis, the merged dataset (`combined_dataset.tsv`) is shuffled to avoid fold based ordering of samples.

- **Script:** `file_shuffle.py`
- **Output:** `combined_dataset_final.tsv`
This file is used as input for all further analyses.


#### 2. Protein Length Distribution
Visualizes the length distribution of positive vs. negative sequences.

- **Script:** `proteinlenght.py`

#### 3. Signal Peptide Length Distribution
Compares the cleavage site lengths of signal peptides between training and benchmark sets, and this provides insights into the typical size range of signal peptides.

- **Script:** `signal_peptide_graph.py`

#### 4. Comparative Amino Acid Composition
We compare the amino acid composition of signal peptides (SPs) against the SwissProt background distribution. This highlights which amino acids are enriched or depleted in SPs relative to general proteins.

- **Script:** `comparison_aa.py`

#### 5. Taxonomic Classification
Analyzes the kingdom-level and species-level distribution of proteins. Useful to detect groups with high False Positive and False Negative risk.

- **Script:** `taxonomic_classification.py`

#### 6. Sequence Logos of SP Cleavage Sites
Extracts motifs around cleavage sites using a window of **[-13, +2]** and generates sequence logos (via WebLogo).

- **Script:**
  - `motif_logo.py` -> motif extraction from fasta files. 
After using this script, we got `.tsv` files which are contain motifs for train and benchmark data. Then we use weblogo service to create sequence logo of SP cleavage sites.
It is useful for showing conserved sequence patterns at cleavage sites. 

