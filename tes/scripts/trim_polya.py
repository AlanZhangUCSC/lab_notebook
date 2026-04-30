#!/usr/bin/env python3
"""Trim the 2nd sequence in an MSA at the poly-A start of the reference."""

import sys
import argparse


def parse_msa(path):
  headers = []
  seqs = []
  current = []
  with open(path) as f:
    for line in f:
      line = line.rstrip("\n")
      if not line:
        continue
      if line.startswith(">"):
        if current:
          seqs.append("".join(current))
          current = []
        headers.append(line[1:])
      else:
        current.append(line)
    if current:
      seqs.append("".join(current))
  return headers, seqs


def find_polyA_start(seq, min_run=10):
  s = seq.lower()
  n = len(s)
  gaps = set("-.")
  i = n - 1
  while i >= 0 and s[i] in gaps:
    i -= 1
  a_count = 0
  last_a = i
  while i >= 0 and (s[i] == "a" or s[i] in gaps):
    if s[i] == "a":
      a_count += 1
      last_a = i
    i -= 1
  if a_count < min_run:
    return -1
  return last_a


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("msa", help="input MSA fasta with 2 sequences")
  ap.add_argument("-o", "--output", default=None,
                  help="output file (default: stdout)")
  ap.add_argument("-m", "--min-run", type=int, default=10,
                  help="minimum poly-A run length (default: 10)")
  ap.add_argument("-w", "--wrap", type=int, default=60,
                  help="line wrap width for output (default: 60)")
  args = ap.parse_args()

  headers, seqs = parse_msa(args.msa)
  if len(seqs) < 2:
    sys.exit("error: expected at least 2 sequences")

  ref, query = seqs[0], seqs[1]
  if len(ref) != len(query):
    sys.exit(f"error: sequences not aligned (lengths {len(ref)} vs {len(query)})")

  polyA_col = find_polyA_start(ref, args.min_run)
  if polyA_col < 0:
    sys.exit(f"error: no poly-A run of >= {args.min_run} found in reference")

  trimmed = query[:polyA_col]

  trimmed = trimmed.upper().replace("-", "")

  out = sys.stdout if args.output is None else open(args.output, "w")
  try:
    out.write(f">{headers[1]}\n")
    if args.wrap > 0:
      for i in range(0, len(trimmed), args.wrap):
        out.write(trimmed[i:i + args.wrap] + "\n")
    else:
      out.write(trimmed + "\n")
    sys.stderr.write(
      f"poly-A starts at alignment column {polyA_col} "
      f"(0-based); trimmed query from {len(query)} to {len(trimmed)} chars\n"
    )
  finally:
    if args.output is not None:
      out.close()


if __name__ == "__main__":
  main()