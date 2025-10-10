#!/usr/bin/env python3

import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import seaborn as sns

scales = {
    "kd": {  # Kyte–Doolittle hydrophobicity
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
        'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    },
    "polarity_Zimmerman": {  # Polarity
        'A': 0.000, 'R': 52.000, 'N': 3.380, 'D': 49.700, 'C': 1.480, 'Q': 3.530,
        'E': 49.900, 'G': 0.000, 'H': 51.600, 'I': 0.130, 'L': 0.130, 'K': 49.500,
        'M': 1.430, 'F': 0.350, 'P': 1.580, 'S': 1.670, 'T': 1.660, 'W': 2.100,
        'Y': 1.610, 'V': 0.130
    },
    "transmembrane_tendency_ZhaoLondon": {  # Transmembrane tendency
        'A': 0.380, 'R': -2.570, 'N': -1.620, 'D': -3.270, 'C': -0.300, 'Q': -1.840,
        'E': -2.900, 'G': -0.190, 'H': -1.440, 'I': 1.970, 'L': 1.820, 'K': -3.460,
        'M': 1.400, 'F': 1.980, 'P': -1.440, 'S': -0.530, 'T': -0.320, 'W': 1.530,
        'Y': 0.490, 'V': 1.460
    },
    "flexibility_BhaskaranPonnuswamy": {  # Flexibility
        'A': 0.360, 'R': 0.530, 'N': 0.460, 'D': 0.510, 'C': 0.350, 'Q': 0.490,
        'E': 0.500, 'G': 0.540, 'H': 0.320, 'I': 0.460, 'L': 0.370, 'K': 0.470,
        'M': 0.300, 'F': 0.310, 'P': 0.510, 'S': 0.510, 'T': 0.440, 'W': 0.310,
        'Y': 0.420, 'V': 0.390
    },
    "helix_ChouFasman": {  # Chou & Fasman helix propensity
        'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70, 'Q': 1.11, 'E': 1.51,
        'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13,
        'P': 0.57, 'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06
    },
    "coil_DeleageRoux": {  # Conformational parameter for coil
        'A': 0.824, 'R': 0.893, 'N': 1.167, 'D': 1.197, 'C': 0.953, 'Q': 0.947,
        'E': 0.761, 'G': 1.251, 'H': 1.068, 'I': 0.886, 'L': 0.810, 'K': 0.897,
        'M': 0.810, 'F': 0.797, 'P': 1.540, 'S': 1.130, 'T': 1.148, 'W': 0.941,
        'Y': 1.109, 'V': 0.772
          
    },
    
    "beta_ChouFasman": {  # Chou & Fasman beta-sheet propensity
    'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 0.54, 'C': 1.19, 'Q': 1.10, 'E': 0.37,
    'G': 0.75, 'H': 0.87, 'I': 1.60, 'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38,
    'P': 0.55, 'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70
},

}

window_map = {"kd":5,"polarity_Zimmerman":5,"transmembrane_tendency_ZhaoLondon":9,"flexibility_BhaskaranPonnuswamy":5,"helix_ChouFasman":9,"coil_DeleageRoux":5, "beta_ChouFasman":9}


def aa_composition(seq, k=22):
    subseq = seq[:k]
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    comp_vector = [subseq.count(aa) / k for aa in amino_acids]
    return comp_vector


def scale_features(seq, scale_dict, n=40, window=5):
    subseq = seq[:n]
    d = window // 2
    padded = "X" * d + subseq + "X" * d
    analysed = ProteinAnalysis(padded)
    scores = analysed.protein_scale(scale_dict, window=window)
    scores = scores[d:-d]
    return np.mean(scores), np.max(scores)

def read_fasta_sequence(fasta_path):
    seq_lines = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq_lines.append(line)
    return "".join(seq_lines)

def extract_features_from_fasta(accession, fasta_dir):
    fasta_file = fasta_dir + accession + ".fasta"
    seq = read_fasta_sequence(fasta_file)

    comp_feats = aa_composition(seq, k=22)

    scale_feats = []
    for scale_name, scale_dict in scales.items():
        w = window_map.get(scale_name,5)
        mean_val, max_val = scale_features(seq, scale_dict, n=40, window=w)
        scale_feats.extend([mean_val, max_val])

    return comp_feats + scale_feats


def build_feature_dataframe(tsv_path, fasta_dir):
    df = pd.read_csv(tsv_path, sep="\t")
    train_df = df[(df["Train"] == True) & (df["Fold"] != -1)]

    all_features = []
    y_labels = []
    folds = []
    accessions = []

    for _, row in train_df.iterrows():
        acc = row["Accession"]
        feats = extract_features_from_fasta(acc, fasta_dir)
        all_features.append(feats)
        y_labels.append(1 if row["Signal_Peptide"] > 0 else 0)
        folds.append(row["Fold"])
        accessions.append(acc)

    #Column names
    comp_cols = [f"comp_{aa}" for aa in "ACDEFGHIKLMNPQRSTVWY"]
    scale_cols = []
    for scale_name in scales.keys():
        scale_cols.extend([f"{scale_name}_mean", f"{scale_name}_max"])

    all_cols = comp_cols + scale_cols

    features_df = pd.DataFrame(all_features, columns=all_cols)
    features_df["Accession"] = accessions
    features_df["Fold"] = folds
    features_df["y"] = y_labels
    return features_df


if __name__ == "__main__":
    tsv_path = "../data_visualization/combined_dataset_final.tsv"
    fasta_dir = "../data_collection/fasta/"
    df_feats = build_feature_dataframe(tsv_path, fasta_dir)
    print(df_feats.head())
    print("\nFinal shape:", df_feats.shape)
    feature_cols = [c for c in df_feats.columns if c not in ["Accession", "Fold", "y"]]
    print("Total features:", len(feature_cols))
    df_feats.to_csv("all_features.tsv", sep="\t", index=False)
    feature_names = np.array(feature_cols)
    print("Unique Accession number:", df_feats["Accession"].nunique())
    X = df_feats[feature_cols].to_numpy()
    y = df_feats["y"].to_numpy()
    folds = df_feats["Fold"].to_numpy()
    accessions = df_feats["Accession"].to_numpy()

    np.savez("np_all_features.npz",X=X,y=y,folds=folds,feature_names=np.array(feature_cols),accessions=accessions)

    print("X shape:",X.shape)
    print("y shape:",y.shape)
    print("fold shape",folds.shape)
    print("accessions shape",accessions.shape)
    print("features shape",feature_names)
    data = np.load("np_all_features.npz",allow_pickle=True)
    print("Keys:",data.files)

    X = data["X"]
    y = data["y"]
    folds=data["folds"]
    accessions=data["accessions"]
    feature_names=data["feature_names"]

    print("First 5 X lines:",X[:5])
    print("First 5 y lines:",y[:5])
    print("First 5 fold lines:",folds[:5])
    print("First 5 accessions lines:",accessions[:5])
    print("First 5 features name lines:",feature_names[:5])




