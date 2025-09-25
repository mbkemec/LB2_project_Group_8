# LB2_project_Group_8

## Dataset Preparation

### Data Collection
Protein sequences were retrieved from the *UniProt* database. Based on selected features:

- *Positive dataset query: (taxonomy_id:2759) AND (reviewed:true) AND (existence:1) AND (fragment:false) AND (ft_signal_exp:*) AND (length:[40 TO *])*
- *Negative dataset query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO ]) NOT (ft_signal:) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))AND (reviewed:true) AND (existence:1)*
  
  All queries were executed on [UniProt](https://www.uniprot.org) on *17.09.2025*
  
---

After collecting the preliminary data, the following filtering steps were applied:

- Proteins with signal peptides shorter than 14 residues were removed.  
- Sequences with specific descriptions were filtered out.

These filtering steps were performed using the following Python scripts:
- [`data_collection_positive.py`](data_collection_positive.py)
- [`data_collection_negative.py`](data_collection_negative.py)

The resulting output files are:
- [`positive.tsv`](positive.tsv) – Filtered positive sequences (TSV format)  
- [`positive.fasta`](positive.fasta) – Filtered positive sequences (FASTA format)
  
- [`negative.tsv`](negative.tsv) – Filtered negative sequences (TSV format)  
- [`negative.fasta`](negative.fasta) – Filtered negative sequences (FASTA format)
---
### Results Summary

|  | Positive | Negative |
|----------|----------|----------|
| Number of Data    | 2949 | 20615   |
| Number of Filtered Data    | 2932   | 20615

| Negative Data | True | False |
|----------|----------|----------|
| Transmembrane Helix   | 1384 | 19231   |
---

## Dataset Pre-processing

The datasets are pre-processed for cross-validation and benchmarking.

### Redundancy Reduction

1. **Clustering**  
   Positive and negative sequences were clustered using **MMseqs2** with thresholds of **≥ 30% sequence identity** and **≥ 40% alignment coverage** to ensure independence between training and test sets.  

The resulting output files are:

- [`pos-cluster-results_cluster.tsv`](pos-cluster-results_cluster.tsv) – Clustered positive sequences (TSV format)  
- [`pos-cluster-results_all_seqs.fasta`](pos-cluster-results_all_seqs.fasta) – Clustered positive sequences (FASTA format)  
- [`neg-cluster-results_cluster.tsv`](neg-cluster-results_cluster.tsv) – Clustered negative sequences (TSV format)  
- [`neg-cluster-results_all_seqs.fasta`](neg-cluster-results_all_seqs.fasta) – Clustered negative sequences (FASTA format)

2. **Selecting Representative Sequences**  
   **MMseqs2** also generated a FASTA file of representative sequences, selecting one sequence from each cluster to:  
   - **Prevent data leakage** – protein families often share similar traits, so related sequences could otherwise be over-represented.  
   - **Reduce overfitting** – uneven family sizes can cause large, highly populated families to dominate the dataset.

The resulting output files are:

- [`pos-cluster-results_rep_seq.fasta`](pos-cluster-results_rep_seq.fasta) – Representative clustered positive sequences (FASTA format)  
- [`neg-cluster-results_rep_seq.fasta`](neg-cluster-results_rep_seq.fasta) – Representative clustered negative sequences (FASTA format)
  
3. **Creating Representative TSV and FASTA Files**  
   After MMseqs2 clustering, the representative FASTA files contain only UniProt IDs.  
   Two processing steps are performed to generate complete, usable datasets:

    a. **Generate Annotated TSVs**  
      Using [`filter_reps.py`](filter_reps.py) – **extracts representative sequences and produces TSV tables** with **Organism**, **Kingdom**, **Length**, and **Transmembrane Helix** information.

      **Inputs**
      - [`positive.fasta`](positive.fasta) – Filtered positive sequences  
      - [`negative.fasta`](negative.fasta) – Filtered negative sequences  

      **Compared With**  
      Representative FASTA files produced by MMseqs2 (one sequence per cluster):
      - [`pos-cluster-results_rep_seq.fasta`](pos-cluster-results_rep_seq.fasta) – Positive cluster representatives  
      - [`neg-cluster-results_rep_seq.fasta`](neg-cluster-results_rep_seq.fasta) – Negative cluster representatives  

      **Outputs**
      - [`positive_NR.tsv`](positive_NR.tsv) – Non-redundant positive dataset with full UniProt annotations  
      - [`negative_NR.tsv`](negative_NR.tsv) – Non-redundant negative dataset with full UniProt annotations  

   b. **Split Representative FASTA Files**  
      Using [`fasta_sep.py`](fasta_sep.py) – **creates one FASTA file per representative sequence**.  
      All individual FASTA files are saved inside the [`fasta`](LB2_project_Group_8/data_collection /fasta/) folder,  
      which is located inside the [`data_collection`](data_collection) directory.


---

4. **Splitting Data into Training and Benchmark Sets**  
   Data is randomly shuffled and divided into two sets using  
   [`shuffle_NR.py`](shuffle_NR.py) – **shuffles and splits data into training/benchmark sets**:  
   - **Training set (80 %)** – for model fitting and hyperparameter tuning  
   - **Benchmarking set (20 %)** – for final, unbiased evaluation  

   **Resulting files**
   - [`train-pos.tsv`](train-pos.tsv)  
   - [`train-neg.tsv`](train-neg.tsv)  
   - [`benchmark-pos.tsv`](benchmark-pos.tsv)  
   - [`benchmark-neg.tsv`](benchmark-neg.tsv)  

---

7. **Building Cross-Validation Subsets**  
   The training set is randomly split into **five** subsets while preserving the positive/negative ratio,  
   using [`split_cross_validation.py`](split_cross_validation.py) – **generates 5-fold cross-validation subsets**.

**Resulting Files**

- **Positive Cross-Validation Sets**
  - [`pos-cv1.tsv`](pos-cv1.tsv) 
  - [`pos-cv2.tsv`](pos-cv2.tsv)
  - [`pos-cv3.tsv`](pos-cv3.tsv) 
  - [`pos-cv4.tsv`](pos-cv4.tsv)  
  - [`pos-cv5.tsv`](pos-cv5.tsv)  

- **Negative Cross-Validation Sets**
  - [`neg-cv1.tsv`](neg-cv1.tsv) 
  - [`neg-cv2.tsv`](neg-cv2.tsv)  
  - [`neg-cv3.tsv`](neg-cv3.tsv) 
  - [`neg-cv4.tsv`](neg-cv4.tsv) 
  - [`neg-cv5.tsv`](neg-cv5.tsv) 
---

### Results Summary  
The final dataset sizes after clustering and splitting are summarized below:

| Dataset                      | Total Samples | Negative | Positive |
|------------------------------|--------------:|--------:|--------:|
| Merged (Train)               |  8021 |  7141 |  874 |
| Merged (Benchmark)           |  2006 |  1787 |  219 |
| Negative (Before Clustering) | 20615 | 20615 |    0 |
| Positive (Before Clustering) |  2932 |     0 | 2932 |
| Negative (After Clustering)  |  8934 |  8934 |    0 |
| Positive (After Clustering)  |  1093 |     0 | 1093 |


---
