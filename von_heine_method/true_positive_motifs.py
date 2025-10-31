import pandas as pd


INPUT_TSV = "benchmark_true_positives.tsv"
FASTA_FOLDER = "../data_collection/fasta"

OUTPUT_BENCH_TSV = "motifs_true_positives.tsv"
WINDOW_LEFT = 12
WINDOW_RIGHT = 3


# Load dataset
df_fn = pd.read_csv(INPUT_TSV, sep="\t")


def extract_motifs(df_subset):
    motifs = []
    for index, row in df_subset.iterrows():
        cleavage_site = int(row['Signal_Peptide'])
        fasta_file = f"{FASTA_FOLDER}/{row['Accession']}.fasta"

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
motifs_bench = extract_motifs(df_fn)

# Save
pd.DataFrame(motifs_bench, columns=["Motif"]).to_csv(OUTPUT_BENCH_TSV, sep="\t", index=False)

print(f" Saved {len(motifs_bench)} false negatives from benchmark motifs to {OUTPUT_BENCH_TSV}")
