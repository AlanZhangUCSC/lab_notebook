import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args():
  p = argparse.ArgumentParser(description="Plot a length distribution from a counts/length TSV.")
  p.add_argument("input", help="Input file: two whitespace-separated columns (count, length).")
  p.add_argument("output", help="Output image path (extension determines format: .png, .pdf, .svg).")
  p.add_argument("--bin-width", type=int, default=10, help="Bin width for histogram (default: 10).")
  p.add_argument("--xlabel", default="Length")
  p.add_argument("--ylabel", default="Count")
  p.add_argument("--title", default="Length distribution")
  p.add_argument("--log-y", action="store_true", help="Use log scale on the y-axis.")
  p.add_argument("--vlines", type=int, nargs="*", default=[50, 100, 150, 200, 250, 300, 350],
                 help="Length values at which to draw vertical reference lines. Pass with no values to disable.")
  p.add_argument("--dpi", type=int, default=300)
  p.add_argument("--survival-suffix", default="_survival",
                 help="Suffix appended to output filename for the survival plot (default: _survival).")
  p.add_argument("--survival-title", default="Proportion of reads with length \u2265 L")
  return p.parse_args()


def load_counts(path):
  lengths, counts = [], []
  with open(path) as fh:
    for line in fh:
      line = line.strip()
      if not line:
        continue
      c, l = line.split()
      lengths.append(int(l))
      counts.append(int(c))
  return np.array(lengths), np.array(counts)


def aggregate_by_length(lengths, counts):
  order = np.argsort(lengths)
  ls = lengths[order]
  cs = counts[order]
  uniq, idx = np.unique(ls, return_inverse=True)
  agg = np.zeros_like(uniq, dtype=np.int64)
  np.add.at(agg, idx, cs)
  return uniq, agg


def survival_curve(lengths, counts):
  uniq, agg = aggregate_by_length(lengths, counts)
  total = agg.sum()
  ge = np.cumsum(agg[::-1])[::-1]
  return uniq, ge / total


def plot_histogram(args, lengths, counts):
  lo = (lengths.min() // args.bin_width) * args.bin_width
  hi = (lengths.max() // args.bin_width + 1) * args.bin_width
  edges = np.arange(lo, hi + args.bin_width, args.bin_width)
  binned, _ = np.histogram(lengths, bins=edges, weights=counts)
  centers = (edges[:-1] + edges[1:]) / 2
  df = pd.DataFrame({"length": centers, "count": binned})

  fig, ax = plt.subplots(figsize=(10, 5))
  sns.barplot(data=df, x="length", y="count", color=sns.color_palette()[0], ax=ax)
  for v in args.vlines:
    if lo <= v < hi:
      pos = (v - lo) / args.bin_width - 0.5
      ax.axvline(pos, color="firebrick", linestyle="--", linewidth=1, alpha=0.7)
      ax.text(pos, ax.get_ylim()[1], f" {v}", color="firebrick",
              va="top", ha="left", fontsize=9)
  ax.set_xlabel(args.xlabel)
  ax.set_ylabel(args.ylabel)
  ax.set_title(args.title)
  if args.log_y:
    ax.set_yscale("log")
  step = max(1, len(centers) // 15)
  ax.set_xticks(range(0, len(centers), step))
  ax.set_xticklabels([f"{int(c)}" for c in centers[::step]], rotation=45, ha="right")
  ax.yaxis.grid(True, linestyle="-", linewidth=0.5, alpha=0.7)
  ax.xaxis.grid(False)
  ax.set_axisbelow(True)
  fig.tight_layout()
  fig.savefig(args.output, dpi=args.dpi)
  plt.close(fig)


def plot_survival(args, lengths, counts):
  xs, ys = survival_curve(lengths, counts)

  fig, ax = plt.subplots(figsize=(10, 5))
  sns.lineplot(x=xs, y=ys, ax=ax, linewidth=2)
  for v in args.vlines:
    if xs[0] <= v <= xs[-1]:
      ax.axvline(v, color="firebrick", linestyle="--", linewidth=1, alpha=0.7)
      ax.text(v, 1.0, f" {v}", color="firebrick", va="top", ha="left", fontsize=9)
  ax.set_xlabel(args.xlabel)
  ax.set_ylabel("Proportion with length \u2265 L")
  ax.set_title(args.survival_title)
  ax.set_ylim(0, 1.02)
  ax.set_yticks(np.arange(0, 1.01, 0.1))
  if args.log_y:
    ax.set_yscale("log")
    ax.set_ylim(max(ys[ys > 0].min() * 0.5, 1e-6), 1.02)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
  ax.yaxis.grid(True, linestyle="-", linewidth=0.5, alpha=0.7)
  ax.xaxis.grid(False)
  ax.set_axisbelow(True)

  root, ext = os.path.splitext(args.output)
  out_path = f"{root}{args.survival_suffix}{ext}"
  fig.tight_layout()
  fig.savefig(out_path, dpi=args.dpi)
  plt.close(fig)


def main():
  args = parse_args()
  sns.set_theme(style="whitegrid", context="talk")
  lengths, counts = load_counts(args.input)
  plot_histogram(args, lengths, counts)
  plot_survival(args, lengths, counts)


if __name__ == "__main__":
  main()