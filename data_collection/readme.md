# LB2_project_Group_8

## Dataset Preparation

### Data Collection
Protein sequences were retrieved from the *UniProt* database. Based on selected features:

- *Positive dataset query: (taxonomy_id:2759) AND (reviewed:true) AND (existence:1) AND (fragment:false) AND (ft_signal_exp:*) AND (length:[40 TO *])*
- *Negative dataset query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO ]) NOT (ft_signal:) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))AND (reviewed:true) AND (existence:1)*
  
  All queries were executed on [UniProt](https://www.uniprot.org) on *17.09.2025*
  
---
#### Preprocessing Dataset

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

  
4. **Creating a New TSV File**  
   A new TSV file containing only the representative sequences is generated. After MMseq2 clustering, the representative sequences FASTA files were only containing UniProt IDs. TSV files with full sequence information for these representatives are generated using the Python script [`filter_reps.py`](filter_reps.py)

Inputs:

-[`positive.fasta`](positive.fasta) – Filtered positive sequences 
-[`negative.fasta`](negative.fasta) – Filtered negative sequences
compared with
-[`pos-cluster-results_rep_seq.fasta`](pos-cluster-results_rep_seq.fasta) – Representative clustered positive sequences
-[`neg-cluster-results_rep_seq.fasta`](neg-cluster-results_rep_seq.fasta) – Representative clustered negative sequences

Outputs:

New files containing only representative sequences and their information obtained from UniProt is created.

-[`positive_NR.tsv`](positive_NR.tsv)
-[`negative_NR.tsv`](negative_NR.tsv)


6. **Splitting Data into Training and Benchmark Sets**  
   Data is randomly divided into:  
   - **Training set (80 %)** – for model fitting and hyperparameter tuning.  
   - **Benchmarking set (20 %)** – for final, unbiased evaluation.

7. **Building Cross-Validation Subsets**  
   The training set is randomly split into **five** subsets while preserving the overall positive/negative ratio for each fold.
---
