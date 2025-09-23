#!/usr/bin/env python3

import random

random.seed(42)

def split_into_folds(input_file, prefix, num_folds=5):

    with open(input_file, "r", encoding="utf-8") as f:
        header = f.readline().strip()
        lines = []
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)

    total_count = len(lines)
    print(f"{input_file} : {total_count} Total line")


    random.shuffle(lines)


    fold_size = total_count // num_folds


    for i in range(num_folds):
        start = i * fold_size

        if i == num_folds - 1:
            end = total_count
        else:
            end = (i + 1) * fold_size

        fold_lines = lines[start:end]


        output_file = f"{prefix}-cv{i+1}.tsv"

        with open(output_file, "w", encoding="utf-8") as out_f:
            out_f.write(header + "\n")
            for line in fold_lines:
                out_f.write(line + "\n")

        print(f"{output_file} : {len(fold_lines)} Line")



if __name__ == "__main__":
    split_into_folds("train-pos.tsv", "pos")

    split_into_folds("train-neg.tsv", "neg")
