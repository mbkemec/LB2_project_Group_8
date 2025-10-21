#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
import joblib



data = np.load("np_all_features.npz", allow_pickle=True)
X = data["X"]
y = data["y"]
folds = data["folds"]
feature_names = data["feature_names"]
accessions = data["accessions"]

meta = pd.read_csv("../data_visualization/combined_dataset_final.tsv", sep="\t")


meta_filtered = meta[meta["Accession"].isin(accessions)].copy()
# Reorder npz file with combined dataset
meta_filtered = meta_filtered.set_index("Accession").loc[accessions].reset_index()


train_mask = meta_filtered["Train"] == True
X_train, y_train = X[train_mask.values], y[train_mask.values]
print(f"Train set: {X_train.shape}")

feature_freq = {
    "comp_A": 20,
    "comp_D": 20,
    "beta_ChouFasman_mean": 20,
    "comp_E": 20,
    "comp_G": 20,
    "comp_L": 20,
    "comp_V": 20,
    "helix_ChouFasman_mean": 20,
    "coil_DeleageRoux_max": 20,
    "coil_DeleageRoux_mean": 20,
    "flexibility_BhaskaranPonnuswamy_mean": 20,
    "flexibility_BhaskaranPonnuswamy_max": 20,
    "transmembrane_tendency_ZhaoLondon_max": 20,
    "transmembrane_tendency_ZhaoLondon_mean": 20,
    "polarity_Zimmerman_mean": 20,
    "polarity_Zimmerman_max": 20,
    "kd_mean": 20,
    "kd_max": 20,
    "helix_ChouFasman_max": 20,
    "beta_ChouFasman_max": 20,
    "comp_R": 19,
    "comp_C": 18,
    "comp_P": 18,
    "comp_N": 17,
    "comp_K": 16,
    "comp_I": 15,
    "comp_F": 13,
    "comp_S": 13,
    "comp_T": 7,
    "comp_Y": 6,
    "comp_Q": 5,
    "comp_W": 5,
    "comp_M": 2,
    "comp_H": 2
}


# Use all available features
selected_features = list(feature_freq.keys())
print(f"Using all {len(selected_features)} features for final model.")

selected_idx = [np.where(feature_names == f)[0][0] for f in selected_features]
X_train_sel = X_train[:, selected_idx]


best_C = 10.0
best_gamma = 0.01

final_model = Pipeline([
    ("scaler", StandardScaler()),
    ("svm", SVC(kernel="rbf", C=best_C, gamma=best_gamma, probability=True, random_state=42))
])

print("\nTraining final SVM model with training set")
final_model.fit(X_train_sel, y_train)



joblib.dump({
    "model": final_model,
    "selected_features": selected_features,
    "best_C": best_C,
    "best_gamma": best_gamma
}, "final_svm_model.pkl")

print("\nFinal model saved as 'final_svm_model.pkl'")

