# Feature Extraction, Selection and SVM

## Feature Extraction
Protein sequences were numerically encoded based on physicochemical properties obtained from the [ProtScale (ExPASy) database](https://web.expasy.org/protscale/).

For each amino acid, the following features were extracted:

- **Kyte–Doolittle Hydrophobicity**  
  Reflects the hydrophobic or hydrophilic nature of amino acids.

- **Polarity (Zimmerman scale)**  
  Measures polarity, contributing to solubility and hydrogen bonding.

- **Transmembrane Tendency (Zhao & London)**  
  Estimates the likelihood of residues being part of transmembrane segments.

- **Flexibility (Bhaskaran & Ponnuswamy)**  
  Represents local conformational flexibility of amino acid residues.

- **Secondary Structure Propensities:**  
  - **Alpha-helix (Chou & Fasman)**  
  - **Coil (Deleage & Roux)**
  - **Beta-sheet (Chou & Fasman)** 

Each protein sequence was converted into a **feature dataframe**, and sliding windows 5 and 9 were applied to capture **local context** across residues.

| Category | Feature Names |
|-----------|----------------|
| **Amino Acid Composition (20)** | `comp_A`, `comp_C`, `comp_D`, `comp_E`, `comp_F`, `comp_G`, `comp_H`, `comp_I`, `comp_K`, `comp_L`, `comp_M`, `comp_N`, `comp_P`, `comp_Q`, `comp_R`, `comp_S`, `comp_T`, `comp_V`, `comp_W`, `comp_Y` |
| **Hydrophobicity (Kyte–Doolittle)** | `kd_mean`, `kd_max` |
| **Polarity (Zimmerman)** | `polarity_Zimmerman_mean`, `polarity_Zimmerman_max` |
| **Transmembrane Tendency (Zhao & London)** | `transmembrane_tendency_ZhaoLondon_mean`, `transmembrane_tendency_ZhaoLondon_max` |
| **Flexibility (Bhaskaran & Ponnuswamy)** | `flexibility_BhaskaranPonnuswamy_mean`, `flexibility_BhaskaranPonnuswamy_max` |
| **Secondary Structure Propensity** | `helix_ChouFasman_mean`, `helix_ChouFasman_max`, `coil_DeleageRoux_mean`, `coil_DeleageRoux_max`, `beta_ChouFasman_mean`, `beta_ChouFasman_max` |


## Feature Selection & SVM

### Overview
This project implements a **nested 5-fold cross-validation** pipeline combining **Random Forest (RF)** feature ranking and **SVM** classification.  
RF provides **Gini-based feature importances**, and SVM hyperparameters (C, γ) are optimized through a focused grid search.

### Fold-wise Performance Metrics
| Outer Fold | Accuracy | MCC | Precision | Recall | F1 |
|:-----------:|:---------:|:----:|:----------:|:--------:|:--------:|
| 0 | 0.964 | 0.813 | 0.841 | 0.826 | 0.833 |
| 1 | 0.978 | 0.887 | 0.915 | 0.884 | 0.899 |
| 2 | 0.975 | 0.870 | 0.905 | 0.864 | 0.884 |
| 3 | 0.977 | 0.883 | 0.887 | 0.904 | 0.895 |
| 4 | 0.971 | 0.852 | 0.866 | 0.871 | 0.868 |

