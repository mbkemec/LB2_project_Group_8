#!/usr/bin/env python3
import random

random.seed(42)
def split_dataset(input_file, train_output, benchmark_output):

    with open(input_file, "r", encoding="utf-8") as f:
        header = f.readline().strip()
        lines = []

        for line in f:
            line = line.strip()
            if line:
                lines.append(line)


    total_count = len(lines)
    print(f"{input_file} : {total_count} total line length")


    random.shuffle(lines)


    split_index = int(total_count * 0.8)

    train_lines = lines[:split_index]
    benchmark_lines = lines[split_index:]


    with open(train_output, "w", encoding="utf-8") as train_f:
        train_f.write(header + "\n")
        for line in train_lines:
            train_f.write(line + "\n")


    with open(benchmark_output, "w", encoding="utf-8") as bench_f:
        bench_f.write(header + "\n")
        for line in benchmark_lines:
            bench_f.write(line + "\n")

    print(f"Training set: {len(train_lines)} line -> {train_output}")
    print(f"Benchmark set: {len(benchmark_lines)} line -> {benchmark_output}")



if __name__ == "__main__":
    split_dataset(input_file="positive_NR.tsv",train_output="train-pos.tsv", benchmark_output="benchmark-pos.tsv")

    split_dataset(input_file="negative_NR.tsv",train_output="train-neg.tsv", benchmark_output="benchmark-neg.tsv")
