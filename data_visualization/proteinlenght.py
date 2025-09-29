import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- Load file ---
df = pd.read_csv("combined_dataset_final.tsv", sep="\t")

# --- Labels ---
df["Label"] = df["Transmembrane_Helix"].map({True: "Positive", False: "Negative"})
df["Set"]   = df["Train"].map({True: "Train", False: "Benchmark"})

# --- Seaborn theme ---
sns.set_theme(style="whitegrid")

def plot_normalized_hist_subplot(df_subset, title, pos):
    plt.subplot(1,2,pos)
    for label, color in zip(["Positive", "Negative"], ["blue", "green"]):
        data = df_subset.loc[df_subset["Label"] == label, "Length"]
        weights = np.ones_like(data) / len(data)   # normalize class-wise
        sns.histplot(
            x=data,
            bins=150,
            element="bars",    
            fill=True,           
            edgecolor=color,      
            linewidth=0.7,       
            alpha=0.4,           
            kde=False,         
            color=color,        
            stat="count",
            weights=weights,
            label=label
        )
    plt.title(title, fontsize=14, weight="bold")
    plt.xlabel("Protein Length (aa)", fontsize=12)
    plt.ylabel("Probability (within class)", fontsize=12)
    plt.xlim(0, 3000)
    plt.legend(title="Label")

# --- Create subplot figure ---
plt.figure(figsize=(14,5))

# Left: Train
plot_normalized_hist_subplot(df[df["Train"] == True], "(Train) Distribution of Protein Lengths", 1)

# Right: Benchmark
plot_normalized_hist_subplot(df[df["Train"] == False], "(Benchmark) Distribution of Protein Lengths", 2)

plt.tight_layout()
plt.savefig('protein_length_comparison.png', dpi=300)
plt.show()

