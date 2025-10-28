#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix, ConfusionMatrixDisplay
)
import joblib
import matplotlib.pyplot as plt
import seaborn as sns


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


cm = confusion_matrix(y_bench, y_pred)
tn, fp, fn, tp = cm.ravel()
print("\nConfusion Matrix:\n", cm)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=["Non-SP", "SP"])
disp.plot(cmap="Blues")
plt.title("Confusion Matrix (Benchmark set)")
plt.xlabel("Predicted label")
plt.ylabel("True label")
plt.tight_layout()
plt.savefig("confusion_matrix.png",dpi=300)
plt.show()
plt.close()

fp_idx = np.where((y_bench == 0) & (y_pred == 1))[0]
fn_idx = np.where((y_bench == 1) & (y_pred == 0))[0]
tp_idx = np.where((y_bench == 1) & (y_pred == 1))[0]

meta_fp = meta_filtered.iloc[fp_idx]
meta_fn = meta_filtered.iloc[fn_idx]
meta_tp = meta_filtered.iloc[tp_idx]

print(f"\nFP count: {len(meta_fp)}, FN count: {len(meta_fn)}, TP count: {len(meta_tp)}")


tm_mask = meta_filtered["Transmembrane_Helix"] == True
non_tm_mask = meta_filtered["Transmembrane_Helix"] == False
neg_mask = y_bench == 0
fp_mask = (y_bench == 0) & (y_pred == 1)

fpr_tm = fp_mask[tm_mask].sum() / neg_mask[tm_mask].sum()
fpr_non_tm = fp_mask[non_tm_mask].sum() / neg_mask[non_tm_mask].sum()

print(f"\nFPR (Transmembrane)     = {fpr_tm*100:.2f}%")
#print(f"FPR (Non-Transmembrane) = {fpr_non_tm*100:.2f}%")

plt.style.use("seaborn-v0_8-darkgrid")

if "Signal_Peptide" in meta.columns:
    tp_len = pd.to_numeric(meta_tp["Signal_Peptide"], errors="coerce").dropna()
    fn_len = pd.to_numeric(meta_fn["Signal_Peptide"], errors="coerce").dropna()

    tp_len = tp_len[tp_len > 0]
    fn_len = fn_len[fn_len > 0]

    print(f"\nTP mean length: {tp_len.mean():.2f} aa")
    print(f"FN mean length: {fn_len.mean():.2f} aa")

    len_df = pd.concat([
        pd.DataFrame({"Length": tp_len, "Group": "True Positive"}),
        pd.DataFrame({"Length": fn_len, "Group": "False Negative"})
    ], ignore_index=True)


    plt.figure(figsize=(8, 6))
    sns.histplot(data=len_df,x="Length",hue="Group",stat="probability",common_norm=False,bins=20)

    plt.xlim(0, 60)
    plt.xlabel("Signal peptide length (aa)")
    plt.ylabel("Probability")
    plt.title("Signal peptide length distribution (TP vs FN)")
    plt.legend(title="Group",labels=["False Negative","True Positive"])
    plt.tight_layout()
    plt.savefig("signal_length_comparison.png",dpi=300)
    plt.show()
    plt.close()

feat_df = pd.DataFrame(X_bench_sel, columns=selected_features)
group_means = {
    "TP": feat_df.iloc[tp_idx].mean(),
    "FP": feat_df.iloc[fp_idx].mean(),
    "FN": feat_df.iloc[fn_idx].mean(),
}
group_means_df = pd.DataFrame(group_means)
print("\n Feature-level mean values (TP vs FP vs FN)")
print(group_means_df.head(10))

# visualize key features
print("\nIt means that negative but model predict positive!!")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
plt.style.use("seaborn-v0_8-darkgrid")

diff_fp_fn = (group_means_df["FP"] - group_means_df["FN"]).sort_values(ascending=False)
top_diff = diff_fp_fn.head(10)
top_diff.plot(kind="barh", ax=axes[0], color="#56B4E9")
axes[0].set_title("Top 10 Feature Differences (FP − FN)")
axes[0].set_xlabel("Mean Difference")
axes[0].invert_yaxis() # That is for the ranking of the features

if "kd_max" in feat_df.columns:
    kd_tp = feat_df.iloc[tp_idx]["kd_max"]
    kd_fp = feat_df.iloc[fp_idx]["kd_max"]
    kd_box = pd.DataFrame({
        "Value": np.concatenate([kd_tp, kd_fp]),
        "Group": ["True Positive"] * len(kd_tp) + ["False Positive"] * len(kd_fp)
    })
    sns.boxplot(data=kd_box, x="Group", y="Value", palette="Set2", ax=axes[1])
    sns.stripplot(data=kd_box, x="Group", y="Value", color="black", alpha=0.4, jitter=True, ax=axes[1])
    axes[1].set_title("Hydrophobicity (Kyte–Doolittle max)")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("kd_max")

if "transmembrane_tendency_ZhaoLondon_max" in feat_df.columns:
    tm_tp = feat_df.iloc[tp_idx]["transmembrane_tendency_ZhaoLondon_max"]
    tm_fp = feat_df.iloc[fp_idx]["transmembrane_tendency_ZhaoLondon_max"]
    tm_box = pd.DataFrame({
        "Value": np.concatenate([tm_tp, tm_fp]),
        "Group": ["True Positive"] * len(tm_tp) + ["False Positive"] * len(tm_fp)
    })
    sns.boxplot(data=tm_box, x="Group", y="Value", palette="Set2", ax=axes[2])
    sns.stripplot(data=tm_box, x="Group", y="Value", color="black", alpha=0.4, jitter=True, ax=axes[2])
    axes[2].set_title("Transmembrane Tendency (Zhao–London max)")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("TM tendency (max)")

plt.tight_layout()
plt.savefig("feature_comparison.png",dpi=300)
plt.show()
plt.close(fig)

print("\nFP/FN Distribution by Kingdom")

print("FP by Kingdom (%):")
print(meta_fp["Kingdom"].value_counts(normalize=True) * 100)

print("\nFN by Kingdom (%):")
print(meta_fn["Kingdom"].value_counts(normalize=True) * 100)
print("-"*50)
print("\nFP/FN Rates by Kingdom")
kingdoms = meta_filtered["Kingdom"].unique()
for k in kingdoms:
    mask = meta_filtered["Kingdom"] == k

    neg_k = (y_bench == 0) & mask
    pos_k = (y_bench == 1) & mask

    fp_k = ((y_pred == 1) & neg_k).sum()
    fn_k = ((y_pred == 0) & pos_k).sum()

    fpr_k = fp_k / neg_k.sum() if neg_k.sum() > 0 else np.nan
    fnr_k = fn_k / pos_k.sum() if pos_k.sum() > 0 else np.nan

    print(f"{k:<10} | FPR={fpr_k*100:5.2f}% | FNR/Recall={fnr_k*100:5.2f}%")

results = {
    "Accuracy": acc,
    "MCC": mcc,
    "Precision": prec,
    "Recall": rec,
    "F1-score": f1,
    "FPR_TM": fpr_tm,
    "FPR_nonTM": fpr_non_tm,
    "FP_count": len(meta_fp),
    "FN_count": len(meta_fn),
    "TP_count": len(meta_tp),
}
#print(results)
