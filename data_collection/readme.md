# LB2_project_Group_8

## Dataset Preparation

### Data Collection
Protein sequences were retrieved from the *UniProt* database. Based on selected features:

- *Positive dataset query: (taxonomy_id:2759) AND (reviewed:true) AND (existence:1) AND (fragment:false) AND (ft_signal_exp:*) AND (length:[40 TO *])*
- *Negative dataset query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO *]) NOT (ft_signal:*) AND ("(cc_scl_term_exp:SL-0091) " OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039)) AND (reviewed:true) AND (existence:1)*
  
  All queries were executed on [UniProt](https://www.uniprot.org) on *17.09.2025*

---
#### Preprocessing Dataset
After getting the preliminary data, the following filtering steps were applied:

- Proteins signal peptide shorter than 14 residues are removed.  
- Peptides with the description are filtered.
  

---

## Results Summary


|  | Positive | Negative |
|----------|----------|----------|
| Number of Data    | 2949 | 20010   |
| Number of Filtered Data    | 2932   | 20010

|  | True | False |
|----------|----------|----------|
| Transmembrane Helix   | 2476 | 17534   |
