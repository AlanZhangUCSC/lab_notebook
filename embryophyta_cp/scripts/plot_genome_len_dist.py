import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "misc_data/genome_len.tsv"
    df = pd.read_csv(path, sep="\t", header=None, names=["accession", "length"])

    lengths = df["length"]
    mean = lengths.mean()
    median = lengths.median()

    fig, ax = plt.subplots(figsize=(10, 5))

    sns.histplot(lengths, kde=True, ax=ax, color="steelblue", edgecolor="white", linewidth=0.4)

    ax.axvline(mean, color="tomato", linestyle="--", linewidth=1.5, label=f"Mean: {mean:,.0f}")
    ax.axvline(median, color="gold", linestyle="--", linewidth=1.5, label=f"Median: {median:,.0f}")

    ax.set_xlabel("Genome Length (bp)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Genome Lengths")
    ax.legend()

    sns.despine()
    plt.tight_layout()

    out = "genome_length_distribution.png"
    plt.savefig(out, dpi=150)
    print(f"Saved: {out}")

if __name__ == "__main__":
    main()