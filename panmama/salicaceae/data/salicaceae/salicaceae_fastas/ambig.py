#!/usr/bin/env python3
"""
ambig.py - identify and naively impute IUPAC ambiguous bases in FASTA.

Subcommands:
  idns    list ambiguous bases and their positions (stdout)
  navimp  naive imputation: N -> A, others -> lexicographically lowest
          canonical base they represent (e.g. Y -> C, R -> A)
"""

import argparse
import sys
import textwrap


def read_fasta(path):
  name, desc, chunks = None, "", []
  with open(path) as fh:
    for line in fh:
      line = line.rstrip("\n")
      if line.startswith(">"):
        if name is not None:
          yield name, desc, "".join(chunks)
        header = line[1:]
        parts = header.split(None, 1)
        name = parts[0] if parts else ""
        desc = header
        chunks = []
      else:
        chunks.append(line.strip())
    if name is not None:
      yield name, desc, "".join(chunks)


def write_fasta(path, records, wrap=80):
  with open(path, "w") as fh:
    for _, desc, seq in records:
      fh.write(f">{desc}\n")
      fh.write("\n".join(textwrap.wrap(seq, wrap)) + "\n")


IUPAC = {
  "R": ("A", "G"),
  "Y": ("C", "T"),
  "S": ("C", "G"),
  "W": ("A", "T"),
  "K": ("G", "T"),
  "M": ("A", "C"),
  "B": ("C", "G", "T"),
  "D": ("A", "G", "T"),
  "H": ("A", "C", "T"),
  "V": ("A", "C", "G"),
  "N": ("A", "C", "G", "T"),
}

NAIVE_MAP = {code: ("A" if code == "N" else bases[0])
              for code, bases in IUPAC.items()}


def cmd_idns(args):
  if args.concise:
    for name, _desc, seq in read_fasta(args.fasta):
      hits = [(i + 1, c) for i, c in enumerate(seq) if c.upper() in IUPAC]
      sites = ",".join(f"{p}:{b}" for p, b in hits)
      print(f"{name}\t{len(hits)}\t{sites}")
    return

  print("record\tposition_1based\tbase\trepresents")
  total = 0
  for name, _desc, seq in read_fasta(args.fasta):
    for i, c in enumerate(seq):
      u = c.upper()
      if u in IUPAC:
        represents = ",".join(IUPAC[u])
        print(f"{name}\t{i + 1}\t{c}\t{represents}")
        total += 1
  print(f"# total ambiguous bases: {total}", file=sys.stderr)


def cmd_navimp(args):
  out_records = []
  total = 0
  for name, desc, seq in read_fasta(args.fasta):
    chars = []
    for c in seq:
      u = c.upper()
      if u in IUPAC:
        repl = NAIVE_MAP[u]
        chars.append(repl if c.isupper() else repl.lower())
        total += 1
      else:
        chars.append(c)
    out_records.append((name, desc, "".join(chars)))
  write_fasta(args.out, out_records)
  print(f"[done] replaced {total} ambiguous bases; wrote {args.out}",
        file=sys.stderr)


def main():
  ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
  sub = ap.add_subparsers(dest="cmd", required=True)

  p1 = sub.add_parser("idns", help="identify ambiguous bases")
  p1.add_argument("fasta")
  p1.add_argument("-c", "--concise", action="store_true",
                  help="one line per record: name<TAB>count<TAB>pos:base,...")
  p1.set_defaults(func=cmd_idns)

  p2 = sub.add_parser("navimp", help="naive imputation to canonical bases")
  p2.add_argument("fasta")
  p2.add_argument("-o", "--out", required=True)
  p2.set_defaults(func=cmd_navimp)

  args = ap.parse_args()
  args.func(args)


if __name__ == "__main__":
  main()