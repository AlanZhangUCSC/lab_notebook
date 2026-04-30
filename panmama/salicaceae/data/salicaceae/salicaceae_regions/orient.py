# Usage:
#   python3 orient.py <alignment.sam> > orientations.tsv 2> warnings.log
#
# Input:
#   SAM file from minimap2, e.g.:
#     minimap2 -a seed_LSC.fasta queries_LSC.fasta > alignment.sam
#
# Output (stdout, TSV):
#   Sample  ref  rev
#   s1      ref  0
#   s2      ref  1
#
#   rev=1 indicates the sample should be reverse-complemented to match the seed.
#
# Warnings (stderr):
#   - Inconsistent orientations across primary mappings for a single query
#   - A query mapping to multiple references
#
# Notes:
#   - Unmapped (flag & 4), secondary (flag & 256), and supplementary (flag & 2048)
#     alignments are ignored. Only primary alignments are considered.

import sys
from collections import defaultdict

def parse_sam(sam_file):
  mappings = defaultdict(list)
  with open(sam_file) as f:
    for line in f:
      if line.startswith('@'):
        continue
      fields = line.rstrip('\n').split('\t')
      if len(fields) < 11:
        continue
      qname = fields[0]
      flag = int(fields[1])
      rname = fields[2]
      if flag & 4:
        continue
      if flag & 2304:
        continue
      is_reverse = 1 if (flag & 16) else 0
      mappings[qname].append((rname, is_reverse))
  return mappings

def main():
  if len(sys.argv) != 2:
    print("Usage: python3 orient.py <alignment.sam>", file=sys.stderr)
    sys.exit(1)

  mappings = parse_sam(sys.argv[1])

  print("Sample\tref\trev")
  for qname, hits in mappings.items():
    refs = {r for r, _ in hits}
    orientations = {o for _, o in hits}

    if len(orientations) > 1:
      print(f"WARNING: {qname} has inconsistent orientations across {len(hits)} primary mappings: {hits}",
            file=sys.stderr)

    if len(refs) > 1:
      print(f"WARNING: {qname} maps to multiple refs: {refs}", file=sys.stderr)

    ref = hits[0][0]
    rev = hits[0][1]
    print(f"{qname}\t{ref}\t{rev}")

if __name__ == "__main__":
  main()