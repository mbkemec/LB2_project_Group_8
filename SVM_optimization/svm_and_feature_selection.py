#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, matthews_corrcoef, precision_score, recall_score, f1_score
from joblib import Parallel, delayed
from collections import Counter


data = np.load("np_all_features.npz", allow_pickle=True)
X = data["X"]
y = data["y"]
folds = data["folds"]
feature_names = data["feature_names"]

print(f"X shape: {X.shape}")
print(f"y shape: {y.shape}")
print(f"folds shape: {folds.shape}")
print(f"Total features: {len(feature_names)}")

def svm_pipeline(C, gamma):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("svm", SVC(kernel="rbf", C=C, gamma=gamma))
    ])

def rf_feature_ranking(X_train, y_train, feature_names):
    rf = RandomForestClassifier(n_estimators=400, random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    importances = rf.feature_importances_
    gini_df = pd.DataFrame({
        "feature": feature_names,
        "importance": importances
    }).sort_values("importance", ascending=False).reset_index(drop=True)
    return gini_df


C_grid = [0.1, 1.0, 10.0, 100.0]
gamma_grid = ["scale", 0.01, 0.1, 1.0]

nested_results = []


def process_inner_fold(outer_fold, inner_fold, X, y, folds, feature_names):
    # Split data for nested CV
    X_test = X[folds == outer_fold]
    y_test = y[folds == outer_fold]
    X_not_test = X[folds != outer_fold]
    y_not_test = y[folds != outer_fold]
    folds_not_test = folds[folds != outer_fold]

    X_val = X_not_test[folds_not_test == inner_fold]
    y_val = y_not_test[folds_not_test == inner_fold]
    X_train = X_not_test[folds_not_test != inner_fold]
    y_train = y_not_test[folds_not_test != inner_fold]

    # RF feature rank
    gini_df = rf_feature_ranking(X_train, y_train, feature_names)
    ranked_features = gini_df["feature"].values

    best_acc = -np.inf
    best_k = None
    best_params = None
    best_features = None

    for k in range(5, X_train.shape[1] + 1):
        subset = ranked_features[:k]
        idx = [np.where(feature_names == f)[0][0] for f in subset]
        Xtr_sel = X_train[:, idx]
        Xva_sel = X_val[:, idx]

        for C in C_grid:
            for gamma in gamma_grid:
                pipe = svm_pipeline(C, gamma)
                pipe.fit(Xtr_sel, y_train)
                val_acc = pipe.score(Xva_sel, y_val)
                if val_acc > best_acc:
                    best_acc = val_acc
                    best_k = k
                    best_params = {"C": C, "gamma": gamma}
                    best_features = subset

    # Retrain with train+val on selected top-k features
    top_idx = [np.where(feature_names == f)[0][0] for f in best_features]
    X_train_val_sel = np.vstack([X_train[:, top_idx], X_val[:, top_idx]])
    y_train_val = np.concatenate([y_train, y_val])
    X_test_sel = X_test[:, top_idx]

    final_pipe = svm_pipeline(best_params["C"], best_params["gamma"])
    final_pipe.fit(X_train_val_sel, y_train_val)
    y_pred_test = final_pipe.predict(X_test_sel)

    acc = accuracy_score(y_test, y_pred_test)
    mcc = matthews_corrcoef(y_test, y_pred_test)
    prec = precision_score(y_test, y_pred_test)
    rec = recall_score(y_test, y_pred_test)
    f1 = f1_score(y_test, y_pred_test)

    print(f"[Outer={outer_fold} | Inner={inner_fold}] Best k={best_k} | Best params={best_params}")

    return {
        "outer_test_fold": outer_fold,
        "inner_val_fold": inner_fold,
        "C": best_params["C"],
        "gamma": best_params["gamma"],
        "ACC": acc,
        "MCC": mcc,
        "Precision": prec,
        "Recall": rec,
        "F1": f1,
        "selected_features": list(best_features)
    }


for outer_fold in range(5):
    print(f"\n--- Outer test fold: {outer_fold} ---")
    inner_folds = [f for f in range(5) if f != outer_fold]
    results = Parallel(n_jobs=-1)(
        delayed(process_inner_fold)(outer_fold, inner_fold, X, y, folds, feature_names)
        for inner_fold in inner_folds
    )
    nested_results.extend(results)


results_df = pd.DataFrame(nested_results)
print("\n--- Fold Results ---")
print(results_df)

avg_results = results_df.groupby("outer_test_fold")[["ACC", "MCC", "Precision", "Recall", "F1"]].mean()
print("\n--- Average Results (Outer Fold) ---")
print(avg_results)


all_selected_features = []
for entry in nested_results:
    all_selected_features.extend(entry["selected_features"])

feature_freq = Counter(all_selected_features)
freq_df = pd.DataFrame.from_dict(feature_freq, orient="index", columns=["Selection_Count"])
freq_df = freq_df.reindex(feature_names)  # preserve original order
freq_df = freq_df.fillna(0).astype(int)
freq_df = freq_df.sort_values("Selection_Count", ascending=False)

print("\n--- Feature Selection Frequency (all features ranked) ---")
print(freq_df)

param_df = results_df[["C", "gamma"]]
print("\n--- Hyperparameter Frequency ---")
print("C frequency:\n", param_df["C"].value_counts())
print("\ngamma frequency:\n", param_df["gamma"].value_counts())

best_C = param_df["C"].mode()[0]
best_gamma = param_df["gamma"].mode()[0]
print(f"\nMost frequently selected params: C={best_C}, gamma={best_gamma}")
