#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, confusion_matrix, roc_curve, auc, average_precision_score
import math
import matplotlib.pyplot as plt
import seaborn as sns


def best_window_score(seq, W, window=15, limit=70):
    best_score = float("-inf")
    limit = min(limit, len(seq))  # if sequence shorter than limit, use seq length
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


def extract_motifs_from_rows(rows, fasta_dir):
    motifs = []
    for _, row in rows.iterrows():
        if row["Signal_Peptide"] > 0:
            fasta_file = fasta_dir + row["Accession"] + ".fasta"
            with open(fasta_file, "r") as f:
                seq = "".join([line.strip() for line in f if not line.startswith(">")])
            start = max(0, int(row["Signal_Peptide"]) - 13 - 1)
            end = min(len(seq), int(row["Signal_Peptide"]) + 2 - 1)
            motif = seq[start:end]
            if motif:
                motifs.append(motif)
    return motifs


def score_rows(rows, W, fasta_dir, window=15, limit=70):
    y_true = []
    scores = []
    for _, row in rows.iterrows():
        fasta_file = fasta_dir + row["Accession"] + ".fasta"
        with open(fasta_file, "r") as f:
            seq = "".join([line.strip() for line in f if not line.startswith(">")])
        s = best_window_score(seq, W, window=window, limit=limit)
        scores.append(s)
        y_true.append(1 if row["Signal_Peptide"] > 0 else 0)
    return y_true, scores


def refit_best_and_eval_benchmark(df, best, fasta_dir="../data_collection/fasta/"):

    test_fold = best["test_fold"]
    val_fold = best["val_fold"]
    train_folds = best["train_folds"]

    retrain_folds = list(train_folds) + [val_fold]

    retrain_df = df[(df["Fold"].isin(retrain_folds)) & (df["Train"] == True)]
    retrain_motifs = extract_motifs_from_rows(retrain_df, fasta_dir)
    print("\n" + "/"*50)
    print("Refitting best model on folds:", retrain_folds,
          "(excluding test fold", test_fold, ")")
    print("Total motifs used for refit:", len(retrain_motifs))

    if len(retrain_motifs) == 0:
        print("No motifs found for refit. Skipping benchmark eval.")
        return

    W_refit = build_pswm(retrain_motifs)

    # benchmark set
    benchmark_df = df[(df["Train"] == False) | (df["Fold"] == -1)]
    y_bench, bench_scores = score_rows(benchmark_df, W_refit, fasta_dir)

    # apply threshold from best CV run
    thr = best["threshold"]
    y_pred_bench = [1 if s >= thr else 0 for s in bench_scores]

    cm = confusion_matrix(y_bench, y_pred_bench, labels=[1, 0])
    TP, FN, FP, TN = cm[0][0], cm[0][1], cm[1][0], cm[1][1]

    total = TP + TN + FP + FN
    acc = (TP + TN) / total if total > 0 else float("nan")
    rec = TP / (TP + FN) if (TP + FN) > 0 else float("nan")
    spec = TN / (TN + FP) if (TN + FP) > 0 else float("nan")
    prec = TP / (TP + FP) if (TP + FP) > 0 else float("nan")
    f1 = (2 * prec * rec) / (prec + rec + 1e-9)
    mcc = ((TP*TN) - (FP*FN)) / math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)+1e-9)
    fpr = FP / (FP + TN) if (FP + TN) > 0 else float("nan")  # NEW

    print("\nBenchmark evaluation:")
    print("Threshold used:", thr)
    print("TP:", TP, "TN:", TN, "FP:", FP, "FN:", FN)
    print("Accuracy:", round(acc, 4))
    print("Recall (Sensitivity):", round(rec, 4))
    print("Specificity:", round(spec, 4))
    print("Precision:", round(prec, 4))
    print("F1 Score:", round(f1, 4))
    print("MCC:", round(mcc, 4))
    print("FPR:", round(fpr, 4))  # NEW
    # Transmembrane FPR
    if "Transmembrane_Helix" in benchmark_df.columns:
        tm_mask = (benchmark_df["Transmembrane_Helix"] == True).values
        non_tm_mask = (benchmark_df["Transmembrane_Helix"] == False).values

        y_bench_arr = np.array(y_bench)
        y_pred_arr = np.array(y_pred_bench)

        neg_mask = (y_bench_arr == 0)       
        fp_mask  = (y_bench_arr == 0) & (y_pred_arr == 1) 

        fpr_tm = fp_mask[tm_mask].sum() / neg_mask[tm_mask].sum() if neg_mask[tm_mask].sum() > 0 else float("nan")
        fpr_non_tm = fp_mask[non_tm_mask].sum() / neg_mask[non_tm_mask].sum() if neg_mask[non_tm_mask].sum() > 0 else float("nan")

        print(f"\nTransmembrane-specific FPRs:")
        print(f"  FPR (Transmembrane)     = {fpr_tm*100:.2f}%")
        print(f"  FPR (Non-Transmembrane) = {fpr_non_tm*100:.2f}%\n")

    # For False Negatives
    fn_mask = [(yt == 1 and yp == 0) for yt, yp in zip(y_bench, y_pred_bench)]
    fn_rows = benchmark_df[fn_mask].copy()
    fn_rows["score"] = np.array(bench_scores)[fn_mask]
    fn_rows["y_true"] = np.array(y_bench)[fn_mask]
    fn_rows["y_pred"] = np.array(y_pred_bench)[fn_mask]

    out_file = "benchmark_false_negatives.tsv"
    fn_rows.to_csv(out_file, sep="\t", index=False)
    print(f"\nSaved false negatives to {out_file} ({len(fn_rows)} sequences)\n")
    print("/"*50 + "\n")

    # For True Positives
    tp_mask = [(yt == 1 and yp == 1) for yt, yp in zip(y_bench, y_pred_bench)]
    tp_rows = benchmark_df[tp_mask].copy()
    tp_rows["score"] = np.array(bench_scores)[tp_mask]
    tp_rows["y_true"] = np.array(y_bench)[tp_mask]
    tp_rows["y_pred"] = np.array(y_pred_bench)[tp_mask]

    out_file_tp = "benchmark_true_positives.tsv"
    tp_rows.to_csv(out_file_tp, sep="\t", index=False)
    print(f"\nSaved true positives to {out_file_tp} ({len(tp_rows)} sequences)\n")
    print("/"*50 + "\n")

	
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

        train_df = df[(df["Train"] == True) & (df["Fold"].isin(train_folds))]
        train_motifs = extract_motifs_from_rows(train_df, fasta_dir)

        print("Training motif count:", len(train_motifs))
        if len(train_motifs) == 0:
            continue
        W = build_pswm(train_motifs)

        val_df = df[(df["Train"] == True) & (df["Fold"] == val_fold)]
        y_val, val_scores = score_rows(val_df, W, fasta_dir)
        precision, recall, thresholds = precision_recall_curve(y_val, val_scores)
        fscore = (2 * precision * recall) / (precision + recall)
        best_index = np.argmax(fscore)
        best_threshold = thresholds[best_index]
        print("Threshold (from validation):", best_threshold)

        test_df = df[(df["Train"] == True) & (df["Fold"] == test_fold)]
        y_test, test_scores = score_rows(test_df, W, fasta_dir)

        y_pred = [1 if s >= best_threshold else 0 for s in test_scores]

        cm = confusion_matrix(y_test, y_pred, labels=[1, 0])
        TP, FN, FP, TN = cm[0][0], cm[0][1], cm[1][0], cm[1][1]

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

        results.append({
            "fold": test_fold,
            "test_fold": test_fold,
            "val_fold": val_fold,
            "train_folds": train_folds,
            "W": W,
            "precision": prec,
            "run": i,
            "TP": TP,
            "TN": TN,
            "FP": FP,
            "FN": FN,
            "threshold": best_threshold,
            "y_test": y_test,
            "test_scores": test_scores,
            "MCC": mcc,
            "accuracy": acc,
            "f1_score": f1,
            "recall": rec
        })

    # pick best by MCC
    best = results[0]
    for r in results:
        if r["MCC"] > best["MCC"]:
            best = r

    print("/"*50)
    print("Best fold for test:", best["test_fold"],
          "Run number:", best["run"]+1,
          "with MCC =", round(best["MCC"], 4))

    metrics = ["accuracy", "recall", "precision", "MCC", "threshold", "f1_score"]
    print("\n\nCross-validation average performance  (mean [ ± SE ]):")
    n = len(results)
    for m in metrics:
        values = [r[m] for r in results]
        mean = np.mean(values)
        std = np.std(values)
        se = std / np.sqrt(n)
        print(f"{m}: {mean:.4f} [± {se:.4f}]")


    all_y = []
    all_scores = []
    for r in results:
        all_y.extend(r["y_test"])
        all_scores.extend(r["test_scores"])


    fpr_curve, tpr_curve, _ = roc_curve(all_y, all_scores)
    roc_auc = auc(fpr_curve, tpr_curve)

    precision_curve, recall_curve, _ = precision_recall_curve(all_y, all_scores)
    avg_prec = average_precision_score(all_y, all_scores)

    # Heatmap of average W
    all_W_dfs = [pd.DataFrame(r["W"]).T for r in results]
    df_pswm_avg = sum(all_W_dfs) / len(all_W_dfs)
    sns.set_theme(style="darkgrid")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].plot(fpr_curve, tpr_curve, color="darkorange", lw=2, label=f"AUC = {roc_auc:.2f}")
    axes[0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0].set_xlabel("False Positive Rate")
    axes[0].set_ylabel("True Positive Rate")
    axes[0].set_title("ROC Curve")
    axes[0].legend(loc="lower right")
    axes[0].grid(True, linestyle="--", alpha=0.6)

    axes[1].plot(recall_curve, precision_curve, color="blue", lw=2, label=f"AP = {avg_prec:.2f}")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title("Precision-Recall Curve")
    axes[1].legend(loc="lower left")
    axes[1].grid(True, linestyle="--", alpha=0.6)

    sns.heatmap(
        df_pswm_avg,
        cmap="YlGnBu",
        center=0,
        annot=True,
        fmt=".2f",
        annot_kws={"size":6},
        cbar_kws={'label': 'Weight (log2 base)'},
        ax=axes[2]
    )
    axes[2].set_xlabel("Motif Position")
    axes[2].set_ylabel("Amino Acid")
    axes[2].set_title("PSWM Heatmap")

    plt.tight_layout()
    plt.savefig("combined_model_metrics.png", dpi=300)
    plt.show()

    refit_best_and_eval_benchmark(df, best, fasta_dir=fasta_dir)


if __name__ == "__main__":
    df = pd.read_csv("../data_visualization/combined_dataset_final.tsv", sep="\t")
    cross_validation(df, fasta_dir="../data_collection/fasta/")


