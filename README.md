## Evaluations

- **SVM Model Evaluation:** Cross-validations and benchmark evaluation results and optimization details are available in the `SVM_optimization` folder.
- **Von Heijne Model Evaluation:** Cross-validations, and benchmark evaluation results are available in the `von_heine_method` folder.

# Signal Peptide Prediction
This repository contains a complete pipeline for **prediction of signal peptides** using sequence based features and Support Vector Machines (SVMs) & Von Heime.

The project consists of **five major stages**:

1. **Data Collection** — UniProt retrieval, filtering, redundancy removal (MMseq2), data splitting
2. **Data Visualization** — Exploratory dataset analysis and biological insight
3. **Von Heime Method** - Firstly Von Heijne method is applied to predict signal peptides.
4. **Feature Extraction, Model Training (SVM) & Optimization** — Numerical representation of sequences, cross-validation, feature selection, SVM training
5. **Benchmark Evaluation & Error Analysis** — Independent testing + biological interpretation

Each stage is documented inside its respective directory, while this README provides a high-level overview.

## Data Collection Step

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


### 2. Protein Length Distribution
Visualizes the length distribution of positive vs. negative sequences.

- **Script:** `proteinlenght.py`

### 3. Signal Peptide Length Distribution
Compares the cleavage site lengths of signal peptides between training and benchmark sets, and this provides insights into the typical size range of signal peptides.

- **Script:** `signal_peptide_graph.py`

### 4. Comparative Amino Acid Composition
We compare the amino acid composition of signal peptides (SPs) against the SwissProt background distribution. This highlights which amino acids are enriched or depleted in SPs relative to general proteins.

- **Script:** `comparison_aa.py`

### 5. Taxonomic Classification
Analyzes the kingdom-level and species-level distribution of proteins. Useful to detect groups with high False Positive and False Negative risk.

- **Script:** `taxonomic_classification.py`

### 6. Sequence Logos of SP Cleavage Sites
Extracts motifs around cleavage sites using a window of **[-13, +2]** and generates sequence logos (via WebLogo).

- **Script:**
  - `motif_logo.py` -> motif extraction from fasta files. 
After using this script, we got `.tsv` files which are contain motifs for train and benchmark data. Then we use weblogo service to create sequence logo of SP cleavage sites.
It is useful for showing conserved sequence patterns at cleavage sites. 

## von Heijne Model – PSWM-Based Signal Peptide Prediction

This folder contains a full implementation of a Position-Specific Weight Matrix (PSWM) model inspired by the classical von Heijne approach for signal peptide (SP) cleavage site recognition.
The workflow consists of:
1. Motif extraction for PSWM construction
2. PSWM model training (log-odds scoring)
3. Sliding-window scoring of full protein sequences
4. 5-fold cross-validation
5. Benchmark evaluation using the best model
6. Extraction of TP/FN motifs for sequence logo generation

### 1. Motif Extraction for Model Training
Training motifs are extracted directly inside `model_creation.py` using `extract_motifs_from_rows()` function. 
For each positive sequence (Signal_Peptide > 0):
The sequence is obtained from `../data_collection/fasta/Accession.fasta`

A fixed [-13, +2] window around the annotated cleavage site is extracted (total motif length = 15 aa)

These motifs form the input for PSWM construction.

Note: This motif extraction is internal and used exclusively for model training.
It is distinct from the benchmark TP/FN motif extraction scripts, which are used only for visualization.

### 2. PSWM Construction
A Position-Specific Weight Matrix is created using `build_pswm(motifs)` function.

Steps:
- Count matrix generated for all 20 amino acids (with pseudocount = 1)
- Convert counts into probabilities
- Convert probabilities into log2-odds weights using SwissProt background frequencies

Result:

A 20 × 15 matrix where each entry reflects how enriched a residue is at a specific motif position relative to background.

### 3. Scoring Protein Sequences Using Sliding Window
Each full protein sequence is scanned using a sliding 15-aa window with `best_window_score(seq, W, window=15, limit=70)` function.

- Score = sum of PSWM weights for each residue in the window
- The highest-scoring window is taken as the sequence score
- Only positions up to 70 aa are scanned (typical SP range)

This produces one score per protein -> used for classification.


### 4. Cross-Validation
Cross-validation are made with `cross_validation(df)` function inside `model_creation.py`.

For each fold:
- Train folds: used to extract motifs + build PSWM
- Validation fold: threshold optimization
- Test fold: final evaluation with that threshold

Reported metrics per fold:
Accuracy
Recall (Sensitivity)
Precision
F1-score
MCC
Confusion matrix

At the end, mean ± standard error (SE) across folds is printed. (The details at `von_heine_method` directory.

### 5. Benchmark Evaluation - Combined Metrics Visualization
The script generates a three-panel figure:

- ROC curve
- Precision–Recall curve
- Average PSWM heatmap across folds

The best CV configuration (highest MCC) is refit on:
- training folds + validation fold
- !(test fold excluded)

Benchmark sequences are then scored, threshold applied, and metrics computed. Additionally, we create `benchmark_false_negatives.tsv` and `benchmark_true_positives.tsv` to use for sequence logo creation.

### 6. Motif Logo for TP and FN

The following scripts extract motifs from benchmark TP and FN sequences: `false_negative_motifs.py` and `true_positive_motifs.py`. Then we use weblogo service to create sequence logo.

## Feature Extraction, Selection and SVM Optimization

The `SVM_optimization/` directory contains all scripts used to:
1. Encode protein sequences into numerical feature vectors
2. Perform nested cross-validation with Random Forest (RF) feature ranking + SVM
3. Train the final SVM model
4. Evaluate it on an independent benchmark set
5. Perform detailed FP/FN error analysis

### 1. Feature Extraction

**Script:** `feature_extraction.py`
**Inputs:**
- `../data_visualization/combined_dataset_final.tsv`
- FASTA files in `../data_collection/fasta/`

For each protein, the N-terminal region is encoded using amino acid composition and several physicochemical scales (from ExPASy / literature).
Sliding windows of size **5** and **9** are used to capture local context.
**Outputs:**
- `all_features.tsv`
- `np_all_features.npz` (X, y, folds, feature_names, accessions)

A similar pipeline is used to generate benchmark features → `benchmark_features.tsv`, `np_benchmark_features.npz`.

### 2. Nested RF + SVM Cross-Validation

**Script:** `svm_and_feature_selection.py`
This script implements a **nested 5-fold cross-validation**:

- **Outer loop:** evaluates generalization (outer test folds)
- **Inner loop:** for each outer fold:
  - Random Forest provides **Gini-based feature importance**
  - Features are ranked; subsets of increasing size (`k`) are tested
  - For each subset, an SVM with RBF kernel is trained with a grid over:
    - C ∈ {0.1, 1.0, 10.0, 100.0}
    - γ ∈ {`scale`, 0.01, 0.1, 1.0}
  - Best `(k, C, γ)` is selected based on inner validation accuracy

#### Cross-Validation Averages (Mean ± SE)

| Metric    | Mean   | ± SE   |
|-----------|--------|--------|
| Accuracy  | 0.973  | 0.003  |
| MCC       | 0.861  | 0.013  |
| Precision | 0.883  | 0.013  |
| Recall    | 0.870  | 0.013  |
| F1        | 0.876  | 0.012  |

Across inner folds, the most frequently selected hyperparameters were:
- **C = 10.0**
- **γ = 0.01**

These are used to define the final SVM model.

### 3. Final Model Training

**Script:** `final_model_creation.py`
**Input:** `np_all_features.npz` (training set only)

Configuration:

- Model: `SVC(kernel="rbf", C=10.0, gamma=0.01, probability=True)`
- Number of features: **34** (all extracted features, because we used every features at least 2 times during cross validation)

The script:
1. Loads all features and labels
2. Trains the SVM on the full training set
3. Saves the model as:

```bash
final_svm_model.pkl
```
### 4. Benchmark Evaluation
**Scripts:** `final_model_test_benchmark.py` & `test_and_visualize_benchmark.py`

These scripts:
1. Load the trained model (`final_svm_model.pkl`)
2. Load benchmark features from `np_benchmark_features.npz`.
3. Evaluate performance on an independent benchmark set (2006 sequences)
#### **Benchmark Metrics**

| Metric      | Score   |
|-------------|---------|
| Accuracy    | 0.9741  |
| MCC         | 0.8657  |
| Precision   | 0.8884  |
| Recall      | 0.8721  |
| F1-score    | 0.8802  |
#### **Confusion Matrix (Benchmark Set)**

- **TP = 191**
- **FP = 24**
- **FN = 28**
- **TN = 1763**

### False Positive Rate (FPR) Comparison Between Models

| Model        | FPR (Overall) | FPR (Transmembrane) |
|--------------|----------------|----------------------|
| **von Heijne** | **2.41%**       | **13.45%**           |
| **SVM**        | **1.34%**       | **~12.6%**               |


### 5. Error Analysis: TP/FN, Length, Feature Patterns, and Composition

**Script:**
- `test_and_visualize_benchmark.py`

This analysis provides detailed insight into **where the SVM model fails** on the benchmark set.

#### **5.1 Transmembrane-Related Errors**

False positives frequently correspond to **N-terminal transmembrane helices** misinterpreted as SPs.

- For TM-annotated proteins, **FPR ≈ 12.6%**

#### **5.2 Feature-Level Differences (TP vs FN)**
Top 10 feature differences (**TP − FN**) highlight the strongest gaps for:

- **kd_max (Kyte–Doolittle hydrophobicity max)**
- **transmembrane_tendency_ZhaoLondon_max**

**Interpretation:**

- TP sequences have **stronger hydrophobic cores**
- FN sequences tend to have **weaker or less typical hydrophobic regions**

NOTE: Other detailed analysis are on the related directory. You can check them in there.

### **Overall Conclusion**

The SVM model achieves **high accuracy and robustness**, both in cross-validation and benchmark evaluation.
Downstream analyses clarify *why misclassifications occur*, highlighting biological factors such as:

- Hydrophobicity strength
- Transmembrane mimicry
- SP length deviations
- Taxonomic diversity

These insights can guide future model improvements.
