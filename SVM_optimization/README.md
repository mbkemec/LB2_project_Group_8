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

Each protein sequence was converted into a **feature dataframe**, and sliding windows 5 and 9 were applied to capture **local context** across residues.

| Category | Feature Names |
|-----------|----------------|
| **Amino Acid Composition (20)** | `comp_A`, `comp_C`, `comp_D`, `comp_E`, `comp_F`, `comp_G`, `comp_H`, `comp_I`, `comp_K`, `comp_L`, `comp_M`, `comp_N`, `comp_P`, `comp_Q`, `comp_R`, `comp_S`, `comp_T`, `comp_V`, `comp_W`, `comp_Y` |
| **Hydrophobicity (Kyte–Doolittle)** | `kd_mean`, `kd_max` |
| **Polarity (Zimmerman)** | `polarity_Zimmerman_mean`, `polarity_Zimmerman_max` |
| **Transmembrane Tendency (Zhao & London)** | `transmembrane_tendency_ZhaoLondon_mean`, `transmembrane_tendency_ZhaoLondon_max` |
| **Flexibility (Bhaskaran & Ponnuswamy)** | `flexibility_BhaskaranPonnuswamy_mean`, `flexibility_BhaskaranPonnuswamy_max` |
| **Secondary Structure Propensity** | `helix_ChouFasman_mean`, `helix_ChouFasman_max`, `coil_DeleageRoux_mean`, `coil_DeleageRoux_max` |


