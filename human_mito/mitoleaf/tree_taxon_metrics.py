#!/usr/bin/env python3
import argparse
import csv
import math
import sys
from collections import Counter

import dendropy


def parse_args():
  p = argparse.ArgumentParser(
    description="Tree-vs-taxonomy congruence: Retention Index, per-taxon F-measure, entropy/MI congruence.")
  p.add_argument("-t", "--tree", required=True, help="Newick file")
  p.add_argument("-m", "--metadata", required=True, help="Delimited metadata file with a header row")
  p.add_argument("-c", "--column", required=True, help="Taxonomy column to evaluate (header name or 0-based index)")
  p.add_argument("--id-column", default=None, help="Sample-ID column (header name or index); default: first column")
  p.add_argument("--delimiter", default="\t", help="Metadata delimiter (default: tab)")
  p.add_argument("--worst", type=int, default=10, help="How many lowest-F taxa to print")
  p.add_argument("--output", default=None, help="Optional TSV: per-taxon F and clade-fragment counts")
  return p.parse_args()


def resolve_column(header, spec):
  if spec in header:
    return header.index(spec)
  if spec.isdigit() and int(spec) < len(header):
    return int(spec)
  sys.exit(f"Column '{spec}' not found in header: {header}")


def load_metadata(path, delimiter, id_spec, tax_spec):
  if delimiter == "\\t":
    delimiter = "\t"
  with open(path, newline="") as fh:
    reader = csv.reader(fh, delimiter=delimiter)
    header = next(reader)
    id_idx = 0 if id_spec is None else resolve_column(header, id_spec)
    tax_idx = resolve_column(header, tax_spec)
    id2label = {}
    for row in reader:
      if len(row) <= max(id_idx, tax_idx):
        continue
      sid = row[id_idx].strip()
      lab = row[tax_idx].strip()
      if sid and lab:
        id2label[sid] = lab
  return id2label


def load_and_prune_tree(path, id2label):
  tree = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
  remove = [lf.taxon for lf in tree.leaf_node_iter()
            if lf.taxon is None or lf.taxon.label.strip() not in id2label]
  dropped = len(remove)
  if remove:
    tree.prune_taxa(remove, suppress_unifurcations=True)
  leaf_label = {}
  for lf in tree.leaf_node_iter():
    leaf_label[lf] = id2label[lf.taxon.label.strip()]
  return tree, leaf_label, dropped


def retention_index(tree, leaf_state, n, state_counts):
  steps = 0
  sset = {}
  for node in tree.postorder_node_iter():
    kids = node.child_nodes()
    if not kids:
      sset[node] = {leaf_state[node]}
      continue
    # Generalized (soft-polytomy) Fitch: a state's tally is the number of
    # child sets containing it; parent keeps the majority states, and the
    # added length is (#children - majority tally).
    tally = {}
    for c in kids:
      for st in sset[c]:
        tally[st] = tally.get(st, 0) + 1
      del sset[c]
    maxc = max(tally.values())
    steps += len(kids) - maxc
    sset[node] = {st for st, v in tally.items() if v == maxc}
  m = len(state_counts) - 1
  g = n - max(state_counts.values())
  denom = g - m
  ri = float("nan") if denom == 0 else (g - steps) / denom
  return ri, steps, g, m


def per_taxon_fmeasure(tree, leaf_label):
  total = Counter(leaf_label.values())
  best = {lab: 0.0 for lab in total}
  counts = {}
  size = {}
  for node in tree.postorder_node_iter():
    kids = node.child_nodes()
    if not kids:
      d = {leaf_label[node]: 1}
      counts[node] = d
      size[node] = 1
      touched = d
    else:
      # Smaller-to-larger merge: reuse the biggest child's dict so total
      # merge work is O(N log N). Only labels whose count changed here can
      # beat their value at a child (denominator only grows upward), so we
      # rescore just those, keeping best-updates O(N log N) too.
      big = max(kids, key=lambda c: len(counts[c]))
      d = counts[big]
      s = size[big]
      touched = {}
      for c in kids:
        if c is big:
          continue
        for lab, ct in counts[c].items():
          d[lab] = d.get(lab, 0) + ct
          touched[lab] = d[lab]
        s += size[c]
        del counts[c]
      counts[node] = d
      size[node] = s
    ssize = size[node]
    for lab in touched:
      f = 2.0 * d[lab] / (ssize + total[lab])
      if f > best[lab]:
        best[lab] = f
  return best, total


def entropy(sizes, n):
  h = 0.0
  for c in sizes:
    p = c / n
    h -= p * math.log(p)
  return h


def congruence(tree, leaf_label, n):
  for node in tree.postorder_node_iter():
    kids = node.child_nodes()
    if not kids:
      node._pl = leaf_label[node]
    else:
      labs = {c._pl for c in kids}
      node._pl = labs.pop() if (len(labs) == 1 and None not in labs) else None
  # A maximal pure clade is a pure node whose parent is not itself inside a
  # pure clade; every tip lands in exactly one such clade (singletons allowed).
  block = {}
  next_id = 0
  leaf_block = {}
  for node in tree.preorder_node_iter():
    par = node.parent_node
    if par is not None and block.get(par) is not None:
      block[node] = block[par]
    elif node._pl is not None:
      block[node] = next_id
      next_id += 1
    else:
      block[node] = None
    if not node.child_nodes():
      leaf_block[node] = block[node]

  tax_sizes = Counter(leaf_label.values())
  blk_sizes = Counter(leaf_block.values())
  frag = Counter()
  seen = {}
  for lf, b in leaf_block.items():
    lab = leaf_label[lf]
    if seen.get((lab, b)) is None:
      seen[(lab, b)] = True
      frag[lab] += 1

  h_tax = entropy(tax_sizes.values(), n)
  h_phy = entropy(blk_sizes.values(), n)
  h_max = math.log(n)
  mi = h_tax  # blocks refine taxa, so I(tax;blocks) = H(tax)
  delta_h = 1.0 if (h_max - h_tax) < 1e-12 else (h_max - h_phy) / (h_max - h_tax)
  nmi_max = 1.0 if h_phy < 1e-12 else mi / h_phy
  nmi_arith = 1.0 if (h_tax + h_phy) < 1e-12 else 2.0 * mi / (h_tax + h_phy)
  return {
    "delta_h": delta_h, "completeness": nmi_max, "nmi": nmi_arith,
    "h_tax": h_tax, "h_phy": h_phy, "h_max": h_max,
    "n_blocks": next_id, "frag": frag}


def main():
  args = parse_args()
  id2label = load_metadata(args.metadata, args.delimiter, args.id_column, args.column)
  if not id2label:
    sys.exit("No usable (id, label) pairs parsed from metadata.")
  tree, leaf_label, dropped = load_and_prune_tree(args.tree, id2label)

  n = len(leaf_label)
  labels = list(leaf_label.values())
  distinct = set(labels)
  if n < 2 or len(distinct) < 2:
    sys.exit(f"Need >=2 tips and >=2 labels after matching (got {n} tips, {len(distinct)} labels).")

  label2id = {lab: i for i, lab in enumerate(distinct)}
  leaf_state = {node: label2id[lab] for node, lab in leaf_label.items()}
  state_counts = Counter(leaf_state.values())

  ri, steps, g, m = retention_index(tree, leaf_state, n, state_counts)
  best, total = per_taxon_fmeasure(tree, leaf_label)
  con = congruence(tree, leaf_label, n)

  weighted = sum(best[lab] * total[lab] for lab in total) / n
  macro = sum(best.values()) / len(best)
  worst_lab = min(best, key=lambda l: best[l])

  print(f"Tips evaluated : {n}  (dropped, unmatched: {dropped})")
  print(f"Taxa at column : {len(distinct)}  ('{args.column}')")
  print()
  print("== Retention Index ==")
  print(f"  RI                 : {ri:.4f}   (steps s={steps}, min m={m}, max g={g})")
  print()
  print("== Per-taxon F-measure ==")
  print(f"  Overall (size-wtd) : {weighted:.4f}")
  print(f"  Macro (unweighted) : {macro:.4f}")
  print(f"  Minimum            : {best[worst_lab]:.4f}  ({worst_lab})")
  k = min(args.worst, len(best))
  print(f"  Lowest {k} taxa:")
  for lab in sorted(best, key=lambda l: (best[l], -total[l]))[:k]:
    print(f"    {best[lab]:.4f}  n={total[lab]:<5d} frags={con['frag'][lab]:<4d} {lab}")
  print()
  print("== Entropy / mutual-information congruence (cut-free) ==")
  print(f"  Delta-H (0..1)     : {con['delta_h']:.4f}")
  print(f"  Completeness/NMImax: {con['completeness']:.4f}")
  print(f"  NMI (arithmetic)   : {con['nmi']:.4f}")
  print(f"  Monophyletic blocks: {con['n_blocks']} (vs {len(distinct)} taxa)")
  print(f"  H_tax={con['h_tax']:.4f}  H_clade={con['h_phy']:.4f}  H_max={con['h_max']:.4f}")

  if args.output:
    with open(args.output, "w", newline="") as fh:
      w = csv.writer(fh, delimiter="\t")
      w.writerow(["taxon", "n_tips", "f_measure", "clade_fragments"])
      for lab in sorted(best):
        w.writerow([lab, total[lab], f"{best[lab]:.6f}", con["frag"][lab]])
    print(f"\nPer-taxon table written to {args.output}")


if __name__ == "__main__":
  main()