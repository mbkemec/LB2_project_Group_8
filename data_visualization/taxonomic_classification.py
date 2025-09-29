
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('combined_dataset_final.tsv', sep="\t")

train_df = df[df['Train'] == True]
bench_df = df[df['Train'] == False]

def plot_side_by_side(data_train, data_bench):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

   # Training Distribution
    kingdom_counts = data_train['Kingdom'].value_counts()
    axes[0].pie(
        kingdom_counts,
        labels=kingdom_counts.index,
        autopct='%1.1f%%',
        colors=plt.cm.Blues(np.linspace(0.4, 0.9, len(kingdom_counts)))
    )
    axes[0].set_title('Training Set – Kingdom Distribution')

    species_counts = data_train['Organism'].value_counts()
    top_species = species_counts.head(5)
    others = species_counts.iloc[5:].sum()
    species_data = pd.concat([top_species, pd.Series({'Others': others})])

    axes[1].pie(
        species_data,
        labels=species_data.index,
        autopct='%1.1f%%',
        colors=plt.cm.Greens(np.linspace(0.4, 0.9, len(species_data)))
    )
    axes[1].set_title('Training Set – Species Distribution')

    plt.tight_layout()
    plt.savefig("taxonomic_classification_training.png", dpi=300)
    plt.show()

    # Benchmarking Distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    kingdom_counts = data_bench['Kingdom'].value_counts()
    axes[0].pie(
        kingdom_counts,
        labels=kingdom_counts.index,
        autopct='%1.1f%%',
        colors=plt.cm.Blues(np.linspace(0.4, 0.9, len(kingdom_counts)))
    )
    axes[0].set_title('Benchmarking Set – Kingdom Distribution')

    species_counts = data_bench['Organism'].value_counts()
    top_species = species_counts.head(5)
    others = species_counts.iloc[5:].sum()
    species_data = pd.concat([top_species, pd.Series({'Others': others})])

    axes[1].pie(
        species_data,
        labels=species_data.index,
        autopct='%1.1f%%',
        colors=plt.cm.Greens(np.linspace(0.4, 0.9, len(species_data)))
    )
    axes[1].set_title('Benchmarking Set – Species Distribution')

    plt.tight_layout()
    plt.savefig("taxonomic_classification_benchmark.png", dpi=300)
    plt.show()

plot_side_by_side(train_df, bench_df)
