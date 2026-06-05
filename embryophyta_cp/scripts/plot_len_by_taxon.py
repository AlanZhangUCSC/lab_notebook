import sys
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def parse_args():
  p = argparse.ArgumentParser(description="Plot genome length distribution per taxon.")
  p.add_argument("meta_file")
  p.add_argument("-o", "--output", default="genome_len_boxplot.png")
  p.add_argument("--box-width", type=float, default=0.6)
  p.add_argument("--ymin", type=float, default=None)
  p.add_argument("--ymax", type=float, default=None)
  p.add_argument("--range-limit", type=int, default=10000)
  p.add_argument("--no-outliers", action="store_true")
  p.add_argument("--figsize", type=float, nargs=2, default=None,
                 help="Width and height in inches. Defaults scale with taxon count.")
  p.add_argument("--log-y", action="store_true")
  p.add_argument("--sort-by", choices=["median", "count", "name"], default="median")
  p.add_argument("--xtick-size", type=float, default=10.0)
  p.add_argument("--xtick-length", type=float, default=4.0)
  p.add_argument("--xtick-width", type=float, default=1.0)
  return p.parse_args()

def load_data(path, range_limit):
  taxon_lens = defaultdict(list)
  with open(path) as f:
    for line in f:
      fields = line.rstrip("\n").split("\t")
      if len(fields) < 8:
        continue
      try:
        genome_len = int(fields[2])
      except ValueError:
        continue
      taxon_lens[fields[7]].append(genome_len)
  return {k: v for k, v in taxon_lens.items() if len(v) > 1 and max(v) - min(v) > range_limit}

def to_long_df(taxon_lens):
  rows = [(t, l) for t, lens in taxon_lens.items() for l in lens]
  return pd.DataFrame(rows, columns=["taxon", "genome_len"])

def order_taxa(df, mode):
  if mode == "median":
    return df.groupby("taxon")["genome_len"].median().sort_values().index.tolist()
  if mode == "count":
    return df["taxon"].value_counts().index.tolist()
  return sorted(df["taxon"].unique())

def main():
  args = parse_args()
  taxon_lens = load_data(args.meta_file, args.range_limit)
  print(f"Loaded {len(taxon_lens)} taxa with >1 genome and genome length range > {args.range_limit} bp", file=sys.stderr)
  if not taxon_lens:
    sys.exit("No taxa with >1 genome found.")
  df = to_long_df(taxon_lens)
  order = order_taxa(df, args.sort_by)

  if args.figsize is None:
    width = max(8, 0.35 * len(order))
    figsize = (width, 6)
  else:
    figsize = tuple(args.figsize)

  # Set theme WITHOUT grid lines (use "white" rather than "whitegrid")
  sns.set_theme(style="white", context="talk")
  fig, ax = plt.subplots(figsize=figsize)

  sns.boxplot(
    data=df, x="taxon", y="genome_len", order=order,
    width=args.box_width,
    showfliers=not args.no_outliers,
    fliersize=3,
    linewidth=1.0,
    palette="viridis",
    ax=ax,
  )

  if args.ymin is not None or args.ymax is not None:
    ax.set_ylim(bottom=args.ymin, top=args.ymax)
  if args.log_y:
    ax.set_yscale("log")

  ax.grid(False)

  ax.set_xlabel("Taxon")
  ax.set_ylabel("Genome length (bp)")
  ax.set_title(f"Genome length distribution per taxon (n={len(order)} taxa)")
  ax.tick_params(
    axis="x",
    which="major",
    bottom=True,
    top=False,
    direction="out",
    length=args.xtick_length,
    width=args.xtick_width,
    labelsize=args.xtick_size,
  )
  plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

  fig.tight_layout()
  fig.savefig(args.output, dpi=600)
  print(f"Wrote {args.output}", file=sys.stderr)

if __name__ == "__main__":
  main()