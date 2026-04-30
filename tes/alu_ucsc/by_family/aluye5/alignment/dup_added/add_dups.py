#!/usr/bin/env python3
"""
Expand alignment files by re-adding sequences that were removed by `seqkit rmdup`.

Inputs:
  - rep_fa:    FASTA of representative (deduplicated) sequences that were aligned.
  - dups_fa:   FASTA containing only the sequences removed during dedup.
  - aln_files: one or more aligned FASTA files derived from rep_fa.

Each duplicate is matched to its representative by exact ungapped, case-insensitive
sequence identity, then written into <input>.expanded.aln carrying the
representative's gapped alignment row.
"""

from __future__ import annotations
import argparse
import os
import sys
from typing import Dict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def canonical(seq: str) -> str:
  return seq.replace("-", "").replace(".", "").replace(" ", "").upper()


def build_seq_to_name(rep_fa: str) -> Dict[str, str]:
  m: Dict[str, str] = {}
  collisions = 0
  for record in SeqIO.parse(rep_fa, "fasta"):
    key = canonical(str(record.seq))
    if key in m and m[key] != record.id:
      collisions += 1
    m[key] = record.id
  if collisions:
    print(f"[warn] {collisions} representative entries shared identical sequences; "
          "later entries overwrote earlier ones in the lookup map.", file=sys.stderr)
  return m


def expand_one(aln_path: str,
               seq_to_name: Dict[str, str],
               dups_fa: str,
               out_path: str) -> tuple[int, int]:
  aligned = SeqIO.to_dict(SeqIO.parse(aln_path, "fasta"))
  records = list(aligned.values())

  added = 0
  unmatched = 0
  for dup in SeqIO.parse(dups_fa, "fasta"):
    key = canonical(str(dup.seq))
    rep_id = seq_to_name.get(key)
    if rep_id is None:
      unmatched += 1
      print(f"[warn] {os.path.basename(aln_path)}: duplicate '{dup.id}' "
            "did not match any representative sequence.", file=sys.stderr)
      continue
    rep_record = aligned.get(rep_id)
    if rep_record is None:
      unmatched += 1
      print(f"[warn] {os.path.basename(aln_path)}: representative '{rep_id}' "
            f"for duplicate '{dup.id}' is absent from this alignment.", file=sys.stderr)
      continue
    records.append(SeqRecord(Seq(str(rep_record.seq)),
                             id=dup.id,
                             description=""))
    added += 1

  with open(out_path, "wt") as fh:
    SeqIO.write(records, fh, "fasta")
  if added == 0:
    print(f"[warn] {os.path.basename(aln_path)}: no duplicates matched any "
          "representative; output is identical to input.", file=sys.stderr)
  return added, unmatched


def main() -> int:
  ap = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("--rep-fa", required=True,
                  help="Deduplicated FASTA used as input to alignment "
                       "(e.g. aluye5_to_consensus.mapped.fa).")
  ap.add_argument("--dups-fa", required=True,
                  help="FASTA containing the sequences that rmdup removed "
                       "(e.g. deduped.mapped.renamed.fa).")
  ap.add_argument("--aln", nargs="+", required=True,
                  help="One or more aligned FASTA files derived from --rep-fa.")
  ap.add_argument("--suffix", default=".expanded.aln",
                  help="Suffix appended to each input alignment for output "
                       "(default: .expanded.aln).")
  args = ap.parse_args()

  for p in [args.rep_fa, args.dups_fa, *args.aln]:
    if not os.path.isfile(p):
      print(f"[error] file not found: {p}", file=sys.stderr)
      return 2

  print(f"[info] indexing representatives from {args.rep_fa} ...", file=sys.stderr)
  seq_to_name = build_seq_to_name(args.rep_fa)
  print(f"[info] indexed {len(seq_to_name)} unique representative sequences.",
        file=sys.stderr)

  total_added = 0
  total_unmatched = 0
  for aln in args.aln:
    base, _ = os.path.splitext(aln)
    out_path = base + args.suffix
    print(f"[info] expanding {aln} -> {out_path}", file=sys.stderr)
    added, unmatched = expand_one(aln, seq_to_name, args.dups_fa, out_path)
    print(f"[info]   added {added} duplicates ({unmatched} unmatched).",
          file=sys.stderr)
    total_added += added
    total_unmatched += unmatched

  print(f"[done] added {total_added} duplicate entries across "
        f"{len(args.aln)} alignment(s); {total_unmatched} unmatched.",
        file=sys.stderr)
  if total_added == 0:
    print("[warn] no duplicates matched any representative in any alignment. "
          "Check that --dups-fa contains the sequences removed from --rep-fa "
          "(not the post-dedup file) and that case/gap conventions match.",
          file=sys.stderr)
  return 0 if total_unmatched == 0 else 1


if __name__ == "__main__":
  sys.exit(main())