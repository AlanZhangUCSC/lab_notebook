#!/usr/bin/env python3
"""Filter Alu annotations by bitscore / aligned-HMM-length density, per family."""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_annotations(path: Path) -> pd.DataFrame:
  df = pd.read_csv(path, sep="\t")
  df["aligned_hmm_len"] = df["hmm_end"] - df["hmm_start"] + 1
  df["score_density"] = df["bitscore"] / df["aligned_hmm_len"]
  return df


def compute_thresholds(df: pd.DataFrame, quantile: float, min_family_size: int,
                       global_fallback_quantile: float) -> dict:
  thresholds = {}
  global_cutoff = df["score_density"].quantile(global_fallback_quantile)
  for family, sub in df.groupby("model_name"):
    if len(sub) >= min_family_size:
      thresholds[family] = sub["score_density"].quantile(quantile)
    else:
      thresholds[family] = global_cutoff
  return thresholds


def apply_filter(df: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
  cutoffs = df["model_name"].map(thresholds)
  return df[df["score_density"] >= cutoffs].copy()


def plot_family_distributions(df: pd.DataFrame, thresholds: dict, out_path: Path,
                              max_cols: int = 4) -> None:
  families = sorted(df["model_name"].unique())
  n = len(families)
  ncols = min(max_cols, n)
  nrows = (n + ncols - 1) // ncols

  fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows),
                           squeeze=False)

  for ax, family in zip(axes.flat, families):
    vals = df.loc[df["model_name"] == family, "score_density"].values
    cutoff = thresholds[family]
    ax.hist(vals, bins=50, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(cutoff, color="crimson", linestyle="--", linewidth=1.2,
               label=f"cutoff = {cutoff:.3f}")
    ax.set_title(f"{family}  (n = {len(vals)})", fontsize=10)
    ax.set_xlabel("bitscore / aligned HMM length")
    ax.set_ylabel("count")
    ax.legend(fontsize=8)

  for ax in axes.flat[n:]:
    ax.set_visible(False)

  fig.tight_layout()
  fig.savefig(out_path, dpi=600)
  plt.close(fig)


def plot_scatter(df: pd.DataFrame, kept_mask: pd.Series, thresholds: dict,
                 out_path: Path, max_cols: int = 4) -> None:
  families = sorted(df["model_name"].unique())
  n = len(families)
  ncols = min(max_cols, n)
  nrows = (n + ncols - 1) // ncols

  fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows),
                           squeeze=False, sharex=True, sharey=True)

  xmax = df["aligned_hmm_len"].max()
  xs = np.linspace(1, xmax, 100)

  for ax, family in zip(axes.flat, families):
    fam_mask = df["model_name"] == family
    fam = df[fam_mask]
    fam_kept = kept_mask & fam_mask
    fam_drop = (~kept_mask) & fam_mask

    ax.scatter(df.loc[fam_drop, "aligned_hmm_len"],
               df.loc[fam_drop, "bitscore"],
               s=5, alpha=0.5, color="lightgray", label="discarded")
    ax.scatter(df.loc[fam_kept, "aligned_hmm_len"],
               df.loc[fam_kept, "bitscore"],
               s=5, alpha=0.6, color="steelblue", label="retained")

    cutoff = thresholds[family]
    ax.plot(xs, cutoff * xs, color="crimson", linestyle="--", linewidth=1.0,
            label=f"density = {cutoff:.3f}")

    ax.set_title(f"{family}  (n = {len(fam)})", fontsize=10)
    ax.set_xlabel("aligned HMM length")
    ax.set_ylabel("bitscore")
    ax.tick_params(axis="both", labelbottom=True, labelleft=True)
    ax.legend(fontsize=7, markerscale=1.8, loc="lower right")

  for ax in axes.flat[n:]:
    ax.set_visible(False)

  fig.tight_layout()
  fig.savefig(out_path, dpi=600)
  plt.close(fig)


def main():
  p = argparse.ArgumentParser(description=__doc__)
  p.add_argument("input", type=Path, help="input .tsv (famdb-style)")
  p.add_argument("-o", "--output", type=Path, required=True,
                 help="filtered annotations output (.tsv)")
  p.add_argument("--hist", type=Path, default=Path("score_density_by_family.png"),
                 help="per-family histogram figure")
  p.add_argument("--scatter", type=Path, default=Path("bitscore_vs_alnlen.png"),
                 help="scatter plot of bitscore vs aligned HMM length")
  p.add_argument("--quantile", type=float, default=0.05,
                 help="discard hits below this per-family score-density "
                      "quantile (default 0.05)")
  p.add_argument("--min-family-size", type=int, default=30,
                 help="families with fewer hits use the global fallback cutoff")
  p.add_argument("--global-quantile", type=float, default=0.05,
                 help="fallback quantile computed across all families")
  p.add_argument("--summary", type=Path, default=None,
                 help="optional per-family summary table (.tsv)")
  args = p.parse_args()

  df = load_annotations(args.input)
  thresholds = compute_thresholds(df, args.quantile, args.min_family_size,
                                  args.global_quantile)

  filtered = apply_filter(df, thresholds)
  kept_mask = df.index.isin(filtered.index)

  out_cols = [c for c in df.columns if c not in ("aligned_hmm_len", "score_density")]
  filtered[out_cols].to_csv(args.output, sep="\t", index=False)

  plot_family_distributions(df, thresholds, args.hist)
  plot_scatter(df, pd.Series(kept_mask, index=df.index), thresholds, args.scatter)

  if args.summary is not None:
    rows = []
    for family, sub in df.groupby("model_name"):
      kept = sub.index.isin(filtered.index).sum()
      rows.append({
        "family": family,
        "n_total": len(sub),
        "n_retained": kept,
        "n_discarded": len(sub) - kept,
        "cutoff_score_density": thresholds[family],
        "median_score_density": sub["score_density"].median(),
      })
    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)

  print(f"input:    {len(df):>8d} hits across {df['model_name'].nunique()} families",
        file=sys.stderr)
  print(f"retained: {len(filtered):>8d} ({100 * len(filtered) / len(df):.1f}%)",
        file=sys.stderr)


if __name__ == "__main__":
  main()