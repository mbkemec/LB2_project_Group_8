#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, confusion_matrix, roc_curve, auc, average_precision_score
import math
import matplotlib.pyplot as plt
import seaborn as sns


def best_window_score(seq, W, window=15, limit=70):
    best_score = float("-inf")
    limit = min(limit, len(seq))  # if sequence shorter than limit, use seq lenght
    for i in range(0, limit - window + 1):
        window_seq = seq[i:i+window]
        score = 0
        for j, aa in enumerate(window_seq):
            if aa in W:
                score += W[aa][j]
        if score > best_score:
            best_score = score
    return best_score


def build_pswm(motifs, pseudocount=1):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    motif_len = len(motifs[0])
    counts = {aa: [pseudocount] * motif_len for aa in aas}
    for motif in motifs:
        for j, aa in enumerate(motif):
            if aa in counts:
                counts[aa][j] += 1
    total = len(motifs) + 20 * pseudocount
    M = {}
    for aa in aas:
        values = []
        for j in range(motif_len):
            val = counts[aa][j] / total
            values.append(val)
        M[aa] = values

    background = {'A': 0.0825, 'R': 0.0553, 'N': 0.0406, 'D': 0.0545, 'C': 0.0137,
                  'Q': 0.0393, 'E': 0.0675, 'G': 0.0707, 'H': 0.0227, 'I': 0.0596,
                  'L': 0.0966, 'K': 0.0584, 'M': 0.0242, 'F': 0.0386, 'P': 0.0470,
                  'S': 0.0656, 'T': 0.0534, 'W': 0.0108, 'Y': 0.0292, 'V': 0.0687}

    W = {}
    for aa in aas:
      values = []
      for j in range(motif_len):
        val = np.log2(M[aa][j] / background[aa])
        values.append(val)
      W[aa] = values


    return W

def cross_validation(df, fasta_dir="../data_collection/fasta/"):
    folds = sorted(df[df["Train"] == True]["Fold"].unique())
    print("Total fold:", len(folds))

    results = []

    for i, test_fold in enumerate(folds):
        val_fold = folds[(i + 1) % len(folds)]
        train_folds = []
        for f in folds:
            if f not in [test_fold, val_fold]:
                train_folds.append(f)


        print("/"*50)
        print(f"Run {i+1}: Test fold = {test_fold}, Validation fold = {val_fold}, Train folds = {train_folds}")

        # TRAIN motifs
        train_df = df[(df["Train"] == True) & (df["Fold"].isin(train_folds))]
        train_motifs = []
        for _, row in train_df.iterrows():
            if row["Signal_Peptide"] > 0:
                fasta_file = fasta_dir + row["Accession"] + ".fasta"
                with open(fasta_file, "r") as f:
                    seq = "".join([line.strip() for line in f if not line.startswith(">")])
                start = max(0, int(row["Signal_Peptide"]) - 13 - 1)
                end = min(len(seq), int(row["Signal_Peptide"]) + 2 - 1)
                motif = seq[start:end]
                if motif:
                    train_motifs.append(motif)
        print("Training motif count:", len(train_motifs))
        if len(train_motifs) == 0:
            continue
        W = build_pswm(train_motifs)

        # VALIDATION
        val_df = df[(df["Train"] == True) & (df["Fold"] == val_fold)]
        y_val, val_scores = [], []
        for _, row in val_df.iterrows():
            fasta_file = fasta_dir + row["Accession"] + ".fasta"
            with open(fasta_file, "r") as f:
                seq = "".join([line.strip() for line in f if not line.startswith(">")])
            if row["Signal_Peptide"] > 0:
                val_scores.append(best_window_score(seq, W, window=15, limit=70))
                y_val.append(1)
            else:
                val_scores.append(best_window_score(seq, W, window=15, limit=70))
                y_val.append(0)
        precision, recall, thresholds = precision_recall_curve(y_val, val_scores)
        fscore = (2 * precision * recall) / (precision + recall)
        best_index = np.argmax(fscore)
        best_threshold = thresholds[best_index]
        print("Threshold (from validation):", best_threshold)

        # Test
        test_df = df[(df["Train"] == True) & (df["Fold"] == test_fold)]
        y_test, test_scores = [], []
        for _, row in test_df.iterrows():
            fasta_file = fasta_dir + row["Accession"] + ".fasta"
            with open(fasta_file, "r") as f:
                seq = "".join([line.strip() for line in f if not line.startswith(">")])
            if row["Signal_Peptide"] > 0:
                test_scores.append(best_window_score(seq, W, window=15, limit=70))
                y_test.append(1)
            else:
                test_scores.append(best_window_score(seq, W, window=15, limit=70))
                y_test.append(0)

        y_pred = []
        for s in test_scores:
            if s >= best_threshold:
                y_pred.append(1)
            else:
                y_pred.append(0)



        cm = confusion_matrix(y_test, y_pred, labels=[1, 0])
        TP, FN, FP, TN = cm[0][0], cm[0][1], cm[1][0], cm[1][1]

        #print(cm)
        print("Confusion matrix (test fold):")
        print("TP:", TP, "TN:", TN, "FP:", FP, "FN:", FN)
        total = TP + TN + FP + FN
        acc = (TP + TN) / total
        rec = TP / (TP + FN)
        spec = TN / (TN + FP)
        prec = TP / (TP + FP)
        f1 = (2 * prec * rec) / (prec + rec + 1e-9)
        mcc = ((TP*TN) - (FP*FN)) / math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        print("Accuracy:", round(acc, 4))
        print("Recall:", round(rec, 4))
        print("Specificity:", round(spec, 4))
        print("Precision:", round(prec, 4))
        print("F1 Score:", round(f1, 4))
        print("MCC:", round(mcc, 4))
        results.append({"fold":test_fold,"W":W,"precision":prec,"run": i,"TP": TP,"TN":TN,"FP":FP,"FN":FN,"threshold":best_threshold,"y_test":y_test,"test_scores":test_scores,"MCC":mcc,"accuracy":acc,"f1_score":f1,"recall":rec})

    best = results[0]
    for r in results:
        if r["MCC"] > best["MCC"]:
            best = r

    print("/"*50)
    print("Best fold for test:", best["fold"],"Run number:",best["run"]+1, "with MCC =", round(best["MCC"], 4))


    metrics = ["accuracy", "recall", "precision", "MCC", "threshold","f1_score"]
    print("\n\nCross-validation average performance  (mean [ ± SE ]):")
    n = len(results) # that is for how many fold we use
    for m in metrics:
        values = []
        for r in results:
            values.append(r[m])

        mean = np.mean(values)
        std = np.std(values)
        se = std / np.sqrt(n)        # standard error
        print(f"{m}: {mean:.4f} [± {se:.4f}]")



    all_y = []
    all_scores = []
    for r in results:
        all_y.extend(r["y_test"])
        all_scores.extend(r["test_scores"])

# ROC curve
    fpr, tpr, _ = roc_curve(all_y, all_scores)
    roc_auc = auc(fpr, tpr)

# PR curve
    precision, recall, _ = precision_recall_curve(all_y, all_scores)
    avg_prec = average_precision_score(all_y, all_scores)
# Heatmap
    all_W_dfs = [pd.DataFrame(r["W"]).T for r in results]
    df_pswm_avg = sum(all_W_dfs) / len(all_W_dfs)
    sns.set_theme(style="darkgrid")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].plot(fpr, tpr, color="darkorange", lw=2, label=f"AUC = {roc_auc:.2f}")
    axes[0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0].set_xlabel("False Positive Rate")
    axes[0].set_ylabel("True Positive Rate")
    axes[0].set_title("ROC Curve")
    axes[0].legend(loc="lower right")
    axes[0].grid(True, linestyle="--", alpha=0.6)

    axes[1].plot(recall, precision, color="blue", lw=2, label=f"AP = {avg_prec:.2f}")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title("Precision-Recall Curve")
    axes[1].legend(loc="lower left")
    axes[1].grid(True, linestyle="--", alpha=0.6)

    sns.heatmap(df_pswm_avg, cmap="YlGnBu", center=0,annot=True,fmt=".2f",annot_kws={"size":6},
            cbar_kws={'label': 'Weight (log2 base)'}, ax=axes[2])
    axes[2].set_xlabel("Motif Position")
    axes[2].set_ylabel("Amino Acid")
    axes[2].set_title("PSWM Heatmap")

    plt.tight_layout()
    plt.savefig("combined_model_metrics.png",dpi=300)
    plt.show()


if __name__ == "__main__":
    df = pd.read_csv("../data_visualization/combined_dataset_final.tsv", sep="\t")
    cross_validation(df, fasta_dir="../data_collection/fasta/")


