# LB2_project_Group_8

## Dataset Preparation

### Data Collection
Protein sequences were retrieved from the *UniProt* database. Based on selected features:

- *Positive dataset query: (taxonomy_id:2759) AND (reviewed:true) AND (existence:1) AND (fragment:false) AND (ft_signal_exp:*) AND (length:[40 TO *])*
- *Negative dataset query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO ]) NOT (ft_signal:) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))AND (reviewed:true) AND (existence:1)*
  
  All queries were executed on [UniProt](https://www.uniprot.org) on *17.09.2025*

---
#### Preprocessing Dataset
After getting the preliminary data, the following filtering steps were applied:

- Proteins signal peptide shorter than 14 residues are removed.  
- Peptides with the description are filtered.
  
- [positive.py](positive.py) – collects and processes positive sequences.
- [negative.py](negative.py) – collects and processes negative sequences.
---

## Results Summary


|  | Positive | Negative |
|----------|----------|----------|
| Number of Data    | 2949 | 20615   |
| Number of Filtered Data    | 2932   | 20615

| Negative Data | True | False |
|----------|----------|----------|
| Transmembrane Helix   | 1384 | 19231   |


## Dataset Pre-processing

The datasets are pre-processed for cross-validation and benchmarking.

### Redundancy Reduction

1. **Clustering**  
   Positive and negative sequences are clustered with **MMseqs2** using thresholds of **≥ 30 % sequence identity** and **≥ 40 % alignment coverage** to ensure independence between training and test sets.

2. **Selecting Representative Sequences**  
   Representative sequences are selected from each cluster to:  
   - **Prevent data leakage** – protein families often share similar traits, so related sequences could otherwise be over-represented.  
   - **Reduce overfitting** – uneven family sizes can cause large, highly populated families to dominate the dataset.

3. **Creating a New TSV File**  
   A new TSV file containing only the representative sequences is generated.

4. **Splitting Data into Training and Benchmark Sets**  
   Data is randomly divided into:  
   - **Training set (80 %)** – for model fitting and hyperparameter tuning.  
   - **Benchmarking set (20 %)** – for final, unbiased evaluation.

5. **Building Cross-Validation Subsets**  
   The training set is randomly split into **five** subsets while preserving the overall positive/negative ratio for each fold.
