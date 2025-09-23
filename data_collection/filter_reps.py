#!/usr/bin/env python3

def load_fasta_accessions(fasta_file):

    accessions = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                acc = line[1:].split()[0] 
                if acc not in accessions:  
                    accessions.append(acc)
    return accessions

def filter_tsv(tsv_file, output_file, keep_accessions):

    with open(tsv_file, "r") as fin, open(output_file, "w") as fout:
        header = fin.readline()
        fout.write(header)  
        for line in fin:
            acc = line.strip().split("\t")[0]
            if acc in keep_accessions:
                fout.write(line)

def main():
    # Pozitives
    pos_fasta = "pos-cluster-results_rep_seq.fasta"
    pos_tsv = "positive.tsv"
    pos_out = "positive_NR.tsv"

    pos_accessions = load_fasta_accessions(pos_fasta)
    filter_tsv(pos_tsv, pos_out, pos_accessions)
    print(f"{len(pos_accessions)} representative positives written to {pos_out}")

    # Negatives
    neg_fasta = "neg-cluster-results_rep_seq.fasta"
    neg_tsv = "negative.tsv"
    neg_out = "negative_NR.tsv"

    neg_accessions = load_fasta_accessions(neg_fasta)
    filter_tsv(neg_tsv, neg_out, neg_accessions)
    print(f"{len(neg_accessions)} representative negatives written to {neg_out}")

if __name__ == "__main__":
    main()
