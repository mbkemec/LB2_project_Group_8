#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, matthews_corrcoef, precision_score, recall_score, f1_score
import joblib



model_data = joblib.load("final_svm_model.pkl")
final_model = model_data["model"]
selected_features = model_data["selected_features"]

data = np.load("np_benchmark_features.npz", allow_pickle=True)
X = data["X"]
y = data["y"]
feature_names = data["feature_names"]
accessions = data["accessions"]

meta = pd.read_csv("../data_visualization/combined_dataset_final.tsv", sep="\t")


meta_filtered = meta[meta["Accession"].isin(accessions)].copy()
meta_filtered = meta_filtered.set_index("Accession").loc[accessions].reset_index()



bench_mask = meta_filtered["Train"] == False
X_bench, y_bench = X[bench_mask.values], y[bench_mask.values]

selected_idx = [np.where(feature_names == f)[0][0] for f in selected_features]
X_bench_sel = X_bench[:, selected_idx]

print(f"Benchmark set shape: {X_bench_sel.shape}")


print("\nUsing model on benchmark dataset...")

y_pred = final_model.predict(X_bench_sel)
acc = accuracy_score(y_bench, y_pred)
mcc = matthews_corrcoef(y_bench, y_pred)
prec = precision_score(y_bench, y_pred)
rec = recall_score(y_bench, y_pred)
f1 = f1_score(y_bench, y_pred)

print("\n/// BENCHMARK RESULTS ///")
print(f"Accuracy : {acc:.4f}")
print(f"MCC       : {mcc:.4f}")
print(f"Precision : {prec:.4f}")
print(f"Recall    : {rec:.4f}")
print(f"F1-score  : {f1:.4f}")


results = {
    "Accuracy": acc,
    "MCC": mcc,
    "Precision": prec,
    "Recall": rec,
    "F1-score": f1
}

