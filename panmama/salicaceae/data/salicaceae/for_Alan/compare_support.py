#!/usr/bin/env python3
import sys
import argparse
import dendropy
import matplotlib.pyplot as plt
import numpy as np


# python3 compare_support.py Salicaceae_cp_final.root_pruned.nwk ../salicaceae_regions/mummer_regions/concatenated/iqtree_results/partition.txt.root_pruned.nwk   --label1 "Trimming A" --label2 "Trimming B" --kde --output support_comparison.png --metric  both

def parse_support_label(label, metric):
  if label is None:
    return None
  label = label.strip()
  if "/" in label:
    parts = label.split("/")
    if metric == "ufboot" and len(parts) >= 2:
      try:
        return float(parts[1])
      except ValueError:
        return None
    elif metric == "alrt" and len(parts) >= 1:
      try:
        return float(parts[0])
      except ValueError:
        return None
  else:
    try:
      return float(label)
    except ValueError:
      return None
  return None

def extract_support_values(tree_path, metric):
  tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
  values = []
  for node in tree.preorder_internal_node_iter():
    if node.parent_node is None or node.is_leaf():
      continue
    v = parse_support_label(node.label, metric)
    if v is not None:
      values.append(v)
  return np.array(values)

def main():
  parser = argparse.ArgumentParser(description="Compare internal node support distributions of two trees.")
  parser.add_argument("tree1", help="First Newick tree file")
  parser.add_argument("tree2", help="Second Newick tree file")
  parser.add_argument("--label1", default="Tree 1", help="Label for first tree")
  parser.add_argument("--label2", default="Tree 2", help="Label for second tree")
  parser.add_argument("--metric", choices=["ufboot", "alrt", "both"], default="ufboot",
                      help="Which support metric to plot (IQ-TREE labels are SH-aLRT/UFBoot)")
  parser.add_argument("--output", "-o", default="support_comparison.png", help="Output image path")
  parser.add_argument("--bins", type=int, default=30, help="Number of histogram bins")
  parser.add_argument("--kde", action="store_true", help="Overlay kernel density estimates")
  args = parser.parse_args()

  metrics_to_plot = ["alrt", "ufboot"] if args.metric == "both" else [args.metric]
  thresholds = {"ufboot": 95, "alrt": 80}
  metric_names = {"ufboot": "UFBoot", "alrt": "SH-aLRT"}

  n_metrics = len(metrics_to_plot)
  fig, axes = plt.subplots(n_metrics, 2, figsize=(12, 5 * n_metrics),
                            gridspec_kw={"width_ratios": [3, 1]}, squeeze=False)

  for row, metric in enumerate(metrics_to_plot):
    v1 = extract_support_values(args.tree1, metric)
    v2 = extract_support_values(args.tree2, metric)

    if len(v1) == 0 or len(v2) == 0:
      sys.exit(f"Error: no parseable {metric_names[metric]} values in one or both trees.")

    ax_hist = axes[row, 0]
    ax_box = axes[row, 1]

    all_vals = np.concatenate([v1, v2])
    bin_edges = np.linspace(all_vals.min(), all_vals.max(), args.bins + 1)

    ax_hist.hist(v1, bins=bin_edges, alpha=0.5, density=True,
                 label=f"{args.label1} (n={len(v1)})", color="#1f77b4",
                 edgecolor="black", linewidth=0.5)
    ax_hist.hist(v2, bins=bin_edges, alpha=0.5, density=True,
                 label=f"{args.label2} (n={len(v2)})", color="#ff7f0e",
                 edgecolor="black", linewidth=0.5)

    if args.kde:
      from scipy.stats import gaussian_kde
      x_grid = np.linspace(all_vals.min(), all_vals.max(), 200)
      for v, color in [(v1, "#1f77b4"), (v2, "#ff7f0e")]:
        if len(v) > 1 and v.std() > 0:
          kde = gaussian_kde(v)
          ax_hist.plot(x_grid, kde(x_grid), color=color, linewidth=2)

    thr = thresholds[metric]
    ax_hist.axvline(thr, color="gray", linestyle="--", linewidth=1, alpha=0.7,
                    label=f"{metric_names[metric]} = {thr}")
    ax_hist.set_xlabel(f"{metric_names[metric]} support")
    ax_hist.set_ylabel("Density")
    ax_hist.set_title(f"Distribution of {metric_names[metric]} support")
    ax_hist.legend(loc="upper left")
    ax_hist.grid(axis="y", alpha=0.3)

    bp = ax_box.boxplot([v1, v2], labels=[args.label1, args.label2],
                         patch_artist=True, widths=0.6)
    for patch, color in zip(bp["boxes"], ["#1f77b4", "#ff7f0e"]):
      patch.set_facecolor(color)
      patch.set_alpha(0.5)
    ax_box.set_ylabel(f"{metric_names[metric]} support")
    ax_box.set_title("Summary")
    ax_box.grid(axis="y", alpha=0.3)

    print(f"\n=== {metric_names[metric]} ===", file=sys.stderr)
    print(f"{args.label1}: n={len(v1)}, mean={v1.mean():.2f}, median={np.median(v1):.2f}, "
          f"<{thr}: {(v1 < thr).sum()} ({100*(v1 < thr).mean():.1f}%)", file=sys.stderr)
    print(f"{args.label2}: n={len(v2)}, mean={v2.mean():.2f}, median={np.median(v2):.2f}, "
          f"<{thr}: {(v2 < thr).sum()} ({100*(v2 < thr).mean():.1f}%)", file=sys.stderr)

    from scipy.stats import mannwhitneyu, ks_2samp
    u_stat, u_p = mannwhitneyu(v1, v2, alternative="two-sided")
    ks_stat, ks_p = ks_2samp(v1, v2)
    print(f"Mann-Whitney U: stat={u_stat:.1f}, p={u_p:.4g}", file=sys.stderr)
    print(f"Kolmogorov-Smirnov: stat={ks_stat:.3f}, p={ks_p:.4g}", file=sys.stderr)

  plt.tight_layout()
  plt.savefig(args.output, dpi=200, bbox_inches="tight")
  print(f"\nPlot saved to {args.output}", file=sys.stderr)

if __name__ == "__main__":
  main()