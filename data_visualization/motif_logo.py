import pandas as pd


INPUT_TSV = "combined_dataset_final.tsv"   
FASTA_FOLDER = "../data_collection/fasta"
OUTPUT_TRAIN_TSV = "motifs_training.tsv"
OUTPUT_BENCH_TSV = "motifs_benchmark.tsv"
WINDOW_LEFT = 13
WINDOW_RIGHT = 2


# Load dataset
df = pd.read_csv(INPUT_TSV, sep="\t")

# Keep only rows with a valid Signal_Peptide
df = df[df['Signal_Peptide'] > 0]

# Separate training and benchmark
df_train = df[df['Train'] == True]
df_bench = df[df['Train'] == False]

def extract_motifs(df_subset):
    motifs = []
    for index, row in df_subset.iterrows():
        cleavage_site = int(row['Signal_Peptide'])
        fasta_file = f"{FASTA_FOLDER}/{row['Accession']}.fasta"

        # Read FASTA
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            sequence = "".join([l.strip() for l in lines if not l.startswith(">")])

        cleavage_index = cleavage_site - 1
        start = max(0, cleavage_index - WINDOW_LEFT)
        end = cleavage_index + WINDOW_RIGHT
        motif = sequence[start:end].ljust(WINDOW_LEFT + WINDOW_RIGHT, "-")
        motifs.append(motif)
    return motifs

# Extract motifs
motifs_train = extract_motifs(df_train)
motifs_bench = extract_motifs(df_bench)

# Save to separate TSVs 
pd.DataFrame(motifs_train, columns=["Motif"]).to_csv(OUTPUT_TRAIN_TSV, sep="\t", index=False)
pd.DataFrame(motifs_bench, columns=["Motif"]).to_csv(OUTPUT_BENCH_TSV, sep="\t", index=False)

print(f" Saved {len(motifs_train)} training motifs to {OUTPUT_TRAIN_TSV}")
print(f" Saved {len(motifs_bench)} benchmark motifs to {OUTPUT_BENCH_TSV}")
