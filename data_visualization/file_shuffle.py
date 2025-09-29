#!/usr/bin/env python3

import random

random.seed(42)

def main():

    with open("../data_collection/combined_dataset.tsv","r") as f:
        header = f.readline()
        lines = f.readlines()

    random.shuffle(lines)

    with open("combined_dataset_final.tsv","w") as f:
        f.write(header)
        f.writelines(lines)

if __name__ == "__main__":
    main()
