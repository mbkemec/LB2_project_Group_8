#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv("combined_dataset_final.tsv", sep="\t")

df["Signal_Peptide"] = pd.to_numeric(df["Signal_Peptide"])

train_sp = df[(df["Train"] == True) & (df["Signal_Peptide"] > 0)]["Signal_Peptide"]
bench_sp = df[(df["Train"] == False) & (df["Signal_Peptide"] > 0)]["Signal_Peptide"]

print("Training SP count:", len(train_sp))
print(train_sp.describe())
print("Benchmark SP count:", len(bench_sp))
print(bench_sp.describe())


plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
sns.histplot(train_sp, bins=30, stat="probability", alpha=0.3, color="blue", label="Train")
sns.kdeplot(train_sp, color="blue", linewidth=2)
plt.title("Signal Peptide Lengths (Training)")
plt.xlabel("Signal Peptide Length (aa)")
plt.ylabel("Probability")
plt.text(53,0.105,f"Training count: {len(train_sp)}",color="blue")

plt.subplot(1,2,2)
sns.histplot(bench_sp, bins=30, stat="probability", alpha=0.3, color="green", label="Benchmark")
sns.kdeplot(bench_sp, color="green", linewidth=2)
plt.title("Signal Peptide Lengths (Benchmark)")
plt.xlabel("Signal Peptide Length (aa)")
plt.ylabel("Probability")
plt.text(45,0.1,f"Benchmark count: {len(bench_sp)}",color="green")

plt.tight_layout()
plt.savefig("signal_peptide_graph.png", dpi=300)
plt.show()
