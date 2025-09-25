def process_one_file(file_path, train_status, fold_number, file_type):
    rows_from_this_file = []
    with open(file_path, "r") as file:
        header_line = file.readline()
        for line in file:
            columns = line.strip().split("\t")
            accession = columns[0]
            organism = columns[1]
            kingdom = columns[2]
            length = columns[3]

            if file_type == "positive":
                signal_peptide = columns[4] if len(columns) > 4 else "0"
                transmembrane_helix = "None"
            else:
                signal_peptide = "0"
                transmembrane_helix = columns[4] if len(columns) > 4 else "False"

            final_row = [
                accession,
                organism,
                kingdom,
                length,
                signal_peptide,
                transmembrane_helix,
                train_status,
                str(fold_number)
            ]
            rows_from_this_file.append(final_row)
    return rows_from_this_file

def main():
    all_rows_combined = []

    for i in range(5):
        fold_number = i + 1
        pos_file_name = f"pos-cv{fold_number}.tsv"
        neg_file_name = f"neg-cv{fold_number}.tsv"

        pos_rows = process_one_file(pos_file_name, "True", i, "positive")
        all_rows_combined.extend(pos_rows)

        neg_rows = process_one_file(neg_file_name, "True", i, "negative")
        all_rows_combined.extend(neg_rows)

    pos_bench_rows = process_one_file("benchmark-pos.tsv", "False", -1, "positive")
    all_rows_combined.extend(pos_bench_rows)

    neg_bench_rows = process_one_file("benchmark-neg.tsv", "False", -1, "negative")
    all_rows_combined.extend(neg_bench_rows)

    output_file_name = "combined_dataset.tsv"
    final_header = [
        "Accession", "Organism", "Kingdom", "Length",
        "Signal_Peptide", "Transmembrane_Helix", "Train", "Fold"
    ]

    with open(output_file_name, "w") as f:
        f.write("\t".join(final_header) + "\n")
        for row in all_rows_combined:
            f.write("\t".join(map(str, row)) + "\n")

    print(f"Finished! All data has been written to {output_file_name}")
    print(f"A total of {len(all_rows_combined)} rows were processed.")

if __name__ == "__main__":
    main()

