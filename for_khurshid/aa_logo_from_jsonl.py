"""
Sequence logo of amino acid frequencies across leaf nodes of a Taxonium JSONL,
for a user-specified gene and residue range.

Algorithm:
  1. Read JSONL header to extract the AA mutation dictionary, filtered to the
     requested gene and residue window.
  2. Read all node lines into a dict keyed by node_id, recording parent_id,
     name, num_tips, and the window-relevant mutations on the branch leading
     to each node.
  3. Verify that node_1 (root) has the user-supplied reference AA sequence
     in the window.
  4. For each leaf (num_tips == 0 and name does not start with 'node_'),
     walk parent links to the root, reverse the path, then apply mutations
     in order from the root to reconstruct the leaf's window AA sequence.
  5. Tally per-position AA counts with explicit for loops.
  6. Convert to per-position probabilities and emit the logo.
"""

import argparse
import gzip
import json
import sys

import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

AAS = list("ACDEFGHIKLMNPQRSTVWY*")


def parse_args():
  p = argparse.ArgumentParser(
    description="Generate an amino acid sequence logo from a Taxonium "
                "JSONL tree over a specified gene region.")
  p.add_argument("-i", "--input-jsonl", required=True,
                 help="Input Taxonium tree (.jsonl or .jsonl.gz)")
  p.add_argument("-G", "--gene", required=True,
                 help="Gene name as stored in the JSONL (e.g. G, F, N)")
  p.add_argument("-s", "--start", type=int, required=True,
                 help="First residue position, 1-based inclusive")
  p.add_argument("-e", "--end", type=int, required=True,
                 help="Last residue position, 1-based inclusive")
  p.add_argument("-o", "--output-prefix", default="logo",
                 help="Prefix for output files (default: logo)")
  p.add_argument("--include-stop", action="store_true",
                 help="Include stop codons (*) in the logo (default: drop)")
  p.add_argument("--ref-residues", default=None,
                 help="Reference AA string for the window. "
                      "Length must equal end - start + 1.")
  p.add_argument("--ref-protein", default=None,
                 help="Reference AA sequence for the full gene.")
  p.add_argument("--ref-fasta", default=None,
                 help="FASTA file of reference protein sequences. "
                      "Record matching --gene by header is used.")
  p.add_argument("--color-scheme", default="chemistry",
                 help="Logomaker color scheme (default: chemistry)")
  args = p.parse_args()

  if args.end < args.start:
    p.error("--end must be >= --start")
  ref_sources = sum(x is not None for x in
                    (args.ref_residues, args.ref_protein, args.ref_fasta))
  if ref_sources != 1:
    p.error("Specify exactly one of --ref-residues, --ref-protein, --ref-fasta")
  if args.ref_residues is not None and \
     len(args.ref_residues) != args.end - args.start + 1:
    p.error(f"--ref-residues length ({len(args.ref_residues)}) must equal "
            f"window length ({args.end - args.start + 1})")
  return args


def open_jsonl(path):
  if path.endswith(".gz"):
    return gzip.open(path, "rt")
  return open(path, "r")


def load_ref_from_fasta(path, gene):
  with open(path) as fh:
    seqs = {}
    name = None
    buf = []
    for line in fh:
      line = line.rstrip()
      if line.startswith(">"):
        if name is not None:
          seqs[name] = "".join(buf)
        name = line[1:].strip()
        buf = []
      else:
        buf.append(line)
    if name is not None:
      seqs[name] = "".join(buf)

  for header, seq in seqs.items():
    first_token = header.split()[0] if header else ""
    if first_token == gene or gene in header.split():
      return seq.replace("*", "").upper()
  raise ValueError(f"No FASTA record matching gene '{gene}' found in {path}. "
                   f"Available headers: {list(seqs.keys())}")


def slice_full_protein(full_seq, start, end):
  full_seq = full_seq.replace("*", "").upper()
  if end > len(full_seq):
    raise ValueError(f"Window end ({end}) exceeds reference protein length "
                     f"({len(full_seq)})")
  return full_seq[start - 1:end]


def resolve_reference(args):
  if args.ref_residues is not None:
    return args.ref_residues.upper(), "--ref-residues"
  if args.ref_protein is not None:
    return slice_full_protein(args.ref_protein, args.start, args.end), \
           "--ref-protein"
  full = load_ref_from_fasta(args.ref_fasta, args.gene)
  return slice_full_protein(full, args.start, args.end), \
         f"--ref-fasta ({args.gene})"


def print_reference(args, ref_str, source):
  print("", file=sys.stderr)
  print(f"  Reference {args.gene}:{args.start}-{args.end} (from {source}):",
        file=sys.stderr)
  block = 10
  pos_line = "    "
  seq_line = "    "
  for i in range(0, len(ref_str), block):
    chunk = ref_str[i:i + block]
    pos_label = str(args.start + i)
    pos_line += pos_label.ljust(len(chunk) + 1)
    seq_line += chunk + " "
  print(pos_line.rstrip(), file=sys.stderr)
  print(seq_line.rstrip(), file=sys.stderr)
  print("", file=sys.stderr)


def main():
  args = parse_args()
  window_len = args.end - args.start + 1
  ref_str, ref_source = resolve_reference(args)
  if len(ref_str) != window_len:
    sys.exit(f"Reference length ({len(ref_str)}) != window length "
             f"({window_len}). Check input.")
  ref_window = list(ref_str)
  print_reference(args, ref_str, ref_source)

  print(f"Loading {args.input_jsonl}...", file=sys.stderr)

  relevant_mut = {}
  nodes = {}

  with open_jsonl(args.input_jsonl) as fh:
    header_line = fh.readline()
    header = json.loads(header_line)
    for m in header.get("mutations", []):
      if m.get("type") != "aa":
        continue
      if m.get("gene") != args.gene:
        continue
      pos = m["residue_pos"]
      if args.start <= pos <= args.end:
        idx = pos - args.start
        relevant_mut[m["mutation_id"]] = (idx, m["previous_residue"],
                                          m["new_residue"])
    print(f"  {len(relevant_mut)} window-relevant mutation IDs in header",
          file=sys.stderr)

    for line in fh:
      obj = json.loads(line)
      node_id = obj["node_id"]
      parent_id = obj.get("parent_id")
      name = obj.get("name", "")
      num_tips = obj.get("num_tips", 0)
      mut_ids = obj.get("mutations", [])

      window_muts = []
      for mid in mut_ids:
        hit = relevant_mut.get(mid)
        if hit is not None:
          window_muts.append(hit)

      nodes[node_id] = {
        "parent_id": parent_id,
        "name": name,
        "num_tips": num_tips,
        "muts": window_muts,
      }

  print(f"  Loaded {len(nodes)} nodes", file=sys.stderr)

  root_id = None
  for nid, info in nodes.items():
    if info["name"] == "node_1":
      root_id = nid
      break
  if root_id is None:
    sys.exit("ERROR: no node with name 'node_1' found in JSONL.")
  root = nodes[root_id]
  print(f"  Root 'node_1' has node_id={root_id}, parent_id={root['parent_id']}",
        file=sys.stderr)
  if root["parent_id"] != root_id:
    print(f"  WARNING: root 'node_1' is not self-parented "
          f"(parent_id={root['parent_id']}, expected {root_id}).",
          file=sys.stderr)

  print("Verifying node_1 (root) reconstructs to reference...", file=sys.stderr)
  root_window = list(ref_window)
  for idx, prev_aa, new_aa in root["muts"]:
    root_window[idx] = new_aa
  if root_window != ref_window:
    diffs = [(args.start + i, ref_window[i], root_window[i])
             for i in range(window_len) if root_window[i] != ref_window[i]]
    print(f"  WARNING: root has {len(diffs)} mutations in the window:",
          file=sys.stderr)
    for pos, r, n in diffs:
      print(f"    {args.gene}:{r}{pos}{n}", file=sys.stderr)
    print(f"  Root window after applying root mutations: "
          f"{''.join(root_window)}", file=sys.stderr)
    print("  This usually means the supplied reference does not match "
          "the tree's root sequence. Continuing with supplied reference.",
          file=sys.stderr)
  else:
    print("  OK: root window matches supplied reference.", file=sys.stderr)

  leaf_ids = []
  for nid, info in nodes.items():
    if info["num_tips"] == 1 and not info["name"].startswith("node_"):
      leaf_ids.append(nid)
  print(f"  Identified {len(leaf_ids)} leaves (num_tips==1 and "
        f"name doesn't start with 'node_')", file=sys.stderr)

  print("Reconstructing leaf sequences...", file=sys.stderr)
  leaf_seqs = []

  inconsistencies = 0
  inconsistency_examples = []

  for k, leaf_id in enumerate(leaf_ids):
    if k > 0 and k % 5000 == 0:
      print(f"    {k}/{len(leaf_ids)}", file=sys.stderr)

    path = []
    cur = leaf_id
    while cur != root_id:
      path.append(cur)
      parent = nodes[cur]["parent_id"]
      if parent == cur or parent is None:
        break
      cur = parent
    path.append(root_id)
    path.reverse()

    seq = list(ref_window)
    for nid in path:
      for idx, prev_aa, new_aa in nodes[nid]["muts"]:
        if seq[idx] != prev_aa:
          inconsistencies += 1
          if len(inconsistency_examples) < 5:
            inconsistency_examples.append(
              (nodes[leaf_id]["name"], nid, args.start + idx,
               seq[idx], prev_aa, new_aa))
        seq[idx] = new_aa
    leaf_seqs.append((nodes[leaf_id]["name"], seq))

  print(f"  Reconstructed {len(leaf_seqs)} leaf sequences", file=sys.stderr)
  if inconsistencies > 0:
    print(f"  WARNING: {inconsistencies} mutations had previous_residue "
          f"that did not match the current sequence. First few:",
          file=sys.stderr)
    for leaf_name, nid, pos, found, expected, new in inconsistency_examples:
      print(f"    leaf={leaf_name} node_id={nid} {args.gene}:{pos} "
            f"current={found} expected_prev={expected} new={new}",
            file=sys.stderr)
  else:
    print("  OK: all mutations had consistent previous_residue.",
          file=sys.stderr)

  fa_path = f"{args.output_prefix}.fa"
  with open(fa_path, "w") as fh:
    for name, seq in leaf_seqs:
      fh.write(f">{name}\n")
      fh.write("".join(seq) + "\n")
  print(f"  Wrote {fa_path}", file=sys.stderr)

  print("Counting AAs per position...", file=sys.stderr)
  counts = [[0] * len(AAS) for _ in range(window_len)]
  aa_index = {aa: i for i, aa in enumerate(AAS)}

  for _, seq in leaf_seqs:
    for pos in range(window_len):
      aa = seq[pos]
      j = aa_index.get(aa)
      if j is not None:
        counts[pos][j] += 1

  print("Computing proportions...", file=sys.stderr)
  freqs = [[0.0] * len(AAS) for _ in range(window_len)]
  for pos in range(window_len):
    total = 0
    for j in range(len(AAS)):
      total += counts[pos][j]
    if total > 0:
      for j in range(len(AAS)):
        freqs[pos][j] = counts[pos][j] / total

  freq_df = pd.DataFrame(freqs, columns=AAS,
                         index=range(args.start, args.end + 1))
  freq_df.index.name = "position"
  freq_path = f"{args.output_prefix}_freq.tsv"
  freq_df.to_csv(freq_path, sep="\t", float_format="%.6f")
  print(f"  Wrote {freq_path}", file=sys.stderr)

  print("Generating logo...", file=sys.stderr)
  logo_df = freq_df if args.include_stop else freq_df.drop(
    columns=[c for c in ['*'] if c in freq_df.columns])

  fig, ax = plt.subplots(figsize=(max(8, window_len * 0.3), 3))
  logomaker.Logo(logo_df, ax=ax, color_scheme=args.color_scheme,
                 stack_order="big_on_top")
  tick_start = ((args.start + 4) // 5) * 5
  ax.set_xticks(list(range(tick_start, args.end + 1, 5)))
  ax.set_ylabel("Probability")
  ax.set_xlabel(f"{args.gene} residue")
  ax.set_ylim(0, 1)
  plt.tight_layout()
  png_path = f"{args.output_prefix}.png"
  pdf_path = f"{args.output_prefix}.pdf"
  plt.savefig(png_path, dpi=300)
  plt.savefig(pdf_path)
  print(f"Done. Wrote {png_path}, {pdf_path}, {freq_path}, {fa_path}",
        file=sys.stderr)


if __name__ == "__main__":
  main()