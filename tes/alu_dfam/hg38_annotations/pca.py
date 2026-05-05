#!/usr/bin/env python3
"""
K-mer PCA for Alu (or other short) sequences.

Reads a FASTA file, computes canonical k-mer frequency vectors for each
sequence, runs PCA, and writes PC scores (TSV) and a 2D scatter plot.
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize


def revcomp_int(kmer_int, k):
  rc = 0
  for _ in range(k):
    rc = (rc << 2) | (3 - (kmer_int & 3))
    kmer_int >>= 2
  return rc


def build_canonical_lut(k):
  num_kmers = 1 << (2 * k)
  canonical_set = set()
  for x in range(num_kmers):
    canonical_set.add(min(x, revcomp_int(x, k)))
  canonical_sorted = sorted(canonical_set)
  to_col = {c: i for i, c in enumerate(canonical_sorted)}
  lut = np.empty(num_kmers, dtype=np.int32)
  for x in range(num_kmers):
    c = min(x, revcomp_int(x, k))
    lut[x] = to_col[c]
  return lut, len(canonical_sorted)


_BASE_LUT = np.full(256, 255, dtype=np.uint8)
for _b, _v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
  _BASE_LUT[_b] = _v


def kmer_vector(seq, k, canonical_lut, dim, strip_polya):
  if strip_polya:
    seq = seq.rstrip("Aa").lstrip("Tt")
  if len(seq) < k:
    return np.zeros(dim, dtype=np.float32)
  arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
  digits = _BASE_LUT[arr]
  valid = digits != 255

  if valid.all():
    runs = [digits]
  else:
    runs = []
    in_run = False
    start = 0
    for i in range(len(valid)):
      if valid[i] and not in_run:
        start = i
        in_run = True
      elif not valid[i] and in_run:
        if i - start >= k:
          runs.append(digits[start:i])
        in_run = False
    if in_run and len(valid) - start >= k:
      runs.append(digits[start:])

  counts = np.zeros(dim, dtype=np.float32)
  shifts = (np.arange(k - 1, -1, -1) * 2).astype(np.uint32)
  for run in runs:
    if len(run) < k:
      continue
    windows = np.lib.stride_tricks.sliding_window_view(run, k).astype(np.uint32)
    kmers = (windows << shifts).sum(axis=1)
    cols = canonical_lut[kmers]
    counts += np.bincount(cols, minlength=dim).astype(np.float32)
  return counts


def parse_family(seq_id, pattern):
  m = re.search(pattern, seq_id)
  return m.group(1) if m else "Unknown"


def build_matrix(fasta_path, k, strip_polya, family_pattern):
  canonical_lut, dim = build_canonical_lut(k)
  ids, families, vectors = [], [], []
  skipped = 0
  for rec in SeqIO.parse(fasta_path, "fasta"):
    seq = str(rec.seq)
    if len(seq) < k:
      skipped += 1
      continue
    vec = kmer_vector(seq, k, canonical_lut, dim, strip_polya)
    if vec.sum() == 0:
      skipped += 1
      continue
    ids.append(rec.id)
    families.append(parse_family(rec.description, family_pattern))
    vectors.append(vec)
  if not vectors:
    sys.exit("No usable sequences found.")
  X = np.vstack(vectors)
  X = normalize(X, norm="l1", axis=1)
  return ids, families, X, dim, skipped


def plot_pca(scores, families, var_ratio, out_path, pc_x, pc_y):
  ix, iy = pc_x - 1, pc_y - 1
  unique_fams = sorted(set(families))
  has_labels = unique_fams != ["Unknown"]
  fig, ax = plt.subplots(figsize=(9, 7))
  if has_labels:
    cmap = plt.get_cmap("tab20", max(len(unique_fams), 3))
    for i, fam in enumerate(unique_fams):
      mask = np.array(families) == fam
      ax.scatter(scores[mask, ix], scores[mask, iy],
                 s=6, alpha=0.6, label=fam, color=cmap(i))
    ax.legend(markerscale=2, fontsize=8, loc="best", frameon=False)
  else:
    ax.scatter(scores[:, ix], scores[:, iy], s=6, alpha=0.5)
  ax.set_xlabel(f"PC{pc_x} ({var_ratio[ix] * 100:.1f}%)")
  ax.set_ylabel(f"PC{pc_y} ({var_ratio[iy] * 100:.1f}%)")
  ax.set_title("k-mer PCA")
  fig.tight_layout()
  fig.savefig(out_path, dpi=150)
  plt.close(fig)


def write_scores(ids, families, scores, var_ratio, out_path):
  n_pcs = scores.shape[1]
  header = ["id", "family"] + [f"PC{i+1}" for i in range(n_pcs)]
  with open(out_path, "w") as fh:
    fh.write("\t".join(header) + "\n")
    for sid, fam, row in zip(ids, families, scores):
      fh.write(sid + "\t" + fam + "\t")
      fh.write("\t".join(f"{v:.6f}" for v in row))
      fh.write("\n")
  with open(str(out_path) + ".varexp", "w") as fh:
    for i, v in enumerate(var_ratio):
      fh.write(f"PC{i+1}\t{v:.6f}\n")


def main():
  ap = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("fasta", type=Path, help="Input FASTA file")
  ap.add_argument("-k", "--kmer-size", type=int, default=4,
                  help="K-mer size (default: 4)")
  ap.add_argument("-n", "--n-components", type=int, default=10,
                  help="Number of PCs to compute (default: 10)")
  ap.add_argument("-o", "--out-prefix", type=Path, default=Path("alu_pca"),
                  help="Output prefix (default: alu_pca)")
  ap.add_argument("--strip-polya", action="store_true",
                  help="Strip trailing A's (and leading T's) before counting")
  ap.add_argument("--family-regex", default=r"(AluY[A-Za-z0-9]*|AluS[A-Za-z0-9]*|AluJ[A-Za-z0-9]*)",
                  help="Regex (group 1 = family) applied to FASTA headers")
  ap.add_argument("--pc-x", type=int, default=1, help="PC on x-axis (default: 1)")
  ap.add_argument("--pc-y", type=int, default=2, help="PC on y-axis (default: 2)")
  args = ap.parse_args()

  if args.kmer_size < 2 or args.kmer_size > 8:
    sys.exit("k must be between 2 and 8")
  if max(args.pc_x, args.pc_y) > args.n_components:
    sys.exit("--pc-x / --pc-y exceed --n-components")

  print(f"Reading {args.fasta} ...", file=sys.stderr)
  ids, families, X, dim, skipped = build_matrix(
      args.fasta, args.kmer_size, args.strip_polya, args.family_regex)
  print(f"  sequences: {len(ids)}  features: {dim}  skipped: {skipped}",
        file=sys.stderr)

  n_comp = min(args.n_components, min(X.shape) - 1)
  pca = PCA(n_components=n_comp, svd_solver="auto")
  scores = pca.fit_transform(X)
  print(f"PCA: top {n_comp} PCs explain "
        f"{pca.explained_variance_ratio_.sum() * 100:.1f}% of variance",
        file=sys.stderr)

  scores_path = args.out_prefix.with_suffix(".tsv")
  plot_path = args.out_prefix.with_suffix(".png")
  write_scores(ids, families, scores, pca.explained_variance_ratio_, scores_path)
  plot_pca(scores, families, pca.explained_variance_ratio_,
           plot_path, args.pc_x, args.pc_y)
  print(f"Wrote {scores_path} and {plot_path}", file=sys.stderr)


if __name__ == "__main__":
  main()