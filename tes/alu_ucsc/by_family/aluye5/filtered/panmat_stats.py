#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

METRICS = [
  ("nodes", "Nodes"),
  ("samples", "Samples"),
  ("substitutions", "Substitutions"),
  ("insertions", "Insertions"),
  ("deletions", "Deletions"),
  ("max_depth", "Max depth"),
  ("mean_depth", "Mean depth"),
]

def load(path):
  df = pd.read_csv(path, sep="\t")
  missing = [c for c, _ in METRICS if c not in df.columns]
  if missing:
    sys.exit(f"missing columns in {path}: {missing}")
  df["label"] = df["file"].map(lambda p: Path(p).stem)
  return df

def plot_grid(df, out_path):
  n = len(METRICS)
  ncols = 2
  nrows = (n + ncols - 1) // ncols
  fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 2.8 * nrows))
  axes = axes.flatten()

  x = np.arange(len(df))
  colors = plt.cm.tab10(np.linspace(0, 1, max(len(df), 2)))

  for ax, (col, title) in zip(axes, METRICS):
    vals = df[col].to_numpy()
    ax.bar(x, vals, color=colors[: len(df)])
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=30, ha="right", fontsize=8)
    ax.margins(x=0.02)
    if vals.max() >= 1000:
      ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for xi, v in zip(x, vals):
      ax.text(xi, v, f"{v:g}", ha="center", va="bottom", fontsize=7)
    ymax = vals.max()
    if ymax > 0:
      ax.set_ylim(0, ymax * 1.15)

  for ax in axes[n:]:
    ax.set_visible(False)

  fig.suptitle("PanMAN tree comparison", fontsize=13)
  fig.tight_layout(rect=[0, 0, 1, 0.97])
  fig.savefig(out_path, dpi=150, bbox_inches="tight")
  print(f"wrote {out_path}", file=sys.stderr)

def plot_normalized(df, out_path):
  fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(METRICS)), 4.5))
  width = 0.8 / max(len(df), 1)
  x = np.arange(len(METRICS))

  norm = {}
  for col, _ in METRICS:
    m = df[col].max()
    norm[col] = df[col] / m if m > 0 else df[col] * 0

  for i, (_, row) in enumerate(df.iterrows()):
    offsets = x - 0.4 + width * (i + 0.5)
    vals = [norm[col].iloc[i] for col, _ in METRICS]
    ax.bar(offsets, vals, width=width, label=row["label"])

  ax.set_xticks(x)
  ax.set_xticklabels([t for _, t in METRICS], rotation=20, ha="right")
  ax.set_ylabel("Value (normalized to per-metric max)")
  ax.set_title("PanMAN tree comparison (normalized)")
  ax.legend(fontsize=8, loc="upper right")
  fig.tight_layout()
  fig.savefig(out_path, dpi=150, bbox_inches="tight")
  print(f"wrote {out_path}", file=sys.stderr)

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("tsv", help="TSV from the awk command")
  ap.add_argument("-o", "--out-prefix", default="panman_compare")
  args = ap.parse_args()

  df = load(args.tsv)
  plot_grid(df, f"{args.out_prefix}_grid.png")
  plot_normalized(df, f"{args.out_prefix}_normalized.png")

if __name__ == "__main__":
  main()