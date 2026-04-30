import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import os


def parse_identity_matrix(filepath):
    sample_names = []
    rows = []
    in_matrix = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not in_matrix:
                if line.startswith('## Identity sequences matrix'):
                    in_matrix = True
                continue
            if line == '':
                break
            fields = line.split()
            sample_names.append(fields[0])
            rows.append([float(x) for x in fields[1:]])

    df = pd.DataFrame(rows, index=sample_names, columns=sample_names)
    return df


def reorder_by_genus(df, metadata_path):
    meta = pd.read_csv(metadata_path, sep="\t")

    accession_to_genus = dict(zip(meta["accession"], meta["genus"]))

    def get_genus(label):
        for acc, genus in accession_to_genus.items():
            if label.startswith(acc):
                return genus
        return label

    genus_labels = [(label, get_genus(label)) for label in df.index]
    genus_labels_sorted = sorted(genus_labels, key=lambda x: x[1])
    sorted_labels = [x[0] for x in genus_labels_sorted]
    display_labels = [f"{x[0]} ({x[1]})" for x in genus_labels_sorted]

    df_sorted = df.loc[sorted_labels, sorted_labels]
    df_sorted.index = display_labels
    df_sorted.columns = display_labels

    return df_sorted


def plot_heatmap(df, output_path, title=None):
    n = len(df)
    fig_size = max(10, n * 0.18)
    font_size = max(4, min(9, 120 / n))

    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.88))

    sns.heatmap(
        df,
        ax=ax,
        cmap="Blues",
        vmin=0.0,
        vmax=1.0,
        xticklabels=True,
        yticklabels=True,
        linewidths=0,
        cbar_kws={"label": "Sequence Identity", "shrink": 0.6},
    )

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=font_size, rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=font_size, rotation=0)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Sequence Identity", fontsize=9)

    if title:
        ax.set_title(title, fontsize=11, pad=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot a trimAl sequence identity matrix as a heatmap."
    )
    parser.add_argument("matrix", help="Path to trimAl identity matrix (.txt)")
    parser.add_argument(
        "--metadata", "-m", help="Path to metadata TSV with accession and genus columns"
    )
    parser.add_argument(
        "--output", "-o", help="Output image path (default: <matrix_basename>.png)"
    )
    parser.add_argument("--title", "-t", help="Plot title (optional)")
    args = parser.parse_args()

    if not os.path.isfile(args.matrix):
        sys.exit(f"Error: matrix file not found: {args.matrix}")

    output_path = args.output or os.path.splitext(args.matrix)[0] + ".png"

    print(f"Parsing matrix: {args.matrix}")
    df = parse_identity_matrix(args.matrix)
    print(f"  {len(df)} sequences")

    if args.metadata:
        if not os.path.isfile(args.metadata):
            sys.exit(f"Error: metadata file not found: {args.metadata}")
        print(f"Reordering by genus using: {args.metadata}")
        df = reorder_by_genus(df, args.metadata)

    plot_heatmap(df, output_path, title=args.title)


if __name__ == "__main__":
    main()