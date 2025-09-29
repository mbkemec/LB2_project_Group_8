import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

INPUT_TSV = "combined_dataset_final.tsv"
FASTA_FOLDER = "../data_collection/fasta"
WINDOW_COLUMN = "Signal_Peptide"
TRAIN_COLUMN = "Train" 

# Load dataset and filter valid SPs
df = pd.read_csv(INPUT_TSV, sep="\t")
df = df[df[WINDOW_COLUMN].notnull() & (df[WINDOW_COLUMN] > 0)]

aa_list = list("ACDEFGHIKLMNPQRSTVWY")

def aa_freq(df_subset):
    counts = {aa:0 for aa in aa_list}
    total = 0
    
    for _, row in df_subset.iterrows():
        fasta_file = f"{FASTA_FOLDER}/{row['Accession']}.fasta"
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            sequence = "".join([l.strip() for l in lines if not l.startswith(">")])
        sp_seq = sequence[:int(row[WINDOW_COLUMN])]
        
        for aa in sp_seq:
            if aa in counts:
                counts[aa] += 1
                total += 1
                
    freqs = [counts[aa]/total if total>0 else 0 for aa in aa_list]
    return freqs

# Frequencies
train_freq = aa_freq(df[df[TRAIN_COLUMN]==True])
bench_freq = aa_freq(df[df[TRAIN_COLUMN]==False])
swissprot_freq = [0.078,0.013,0.053,0.063,0.039,0.073,0.022,0.062,0.058,0.092,
                  0.023,0.045,0.051,0.039,0.052,0.070,0.058,0.065,0.013,0.032]

# Plot grouped bars
x = np.arange(len(aa_list))
width = 0.25
fig, ax = plt.subplots(figsize=(12,6))
ax.bar(x - width, train_freq, width, label="Training SPs")
ax.bar(x, bench_freq, width, label="Benchmark SPs")
ax.bar(x + width, swissprot_freq, width, label="SwissProt")
ax.set_xticks(x)
ax.set_xticklabels(aa_list)
ax.set_ylabel("Fraction")
ax.set_xlabel("Amino Acid")
ax.set_title("Comparative Amino-Acid Composition of SPs")
ax.legend()
plt.savefig('comparative_aa_for_sp.png', dpi =300)
plt.show()

