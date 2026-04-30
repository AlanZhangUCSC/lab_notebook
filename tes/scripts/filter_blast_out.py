#!/usr/bin/env python3
"""
Filter BLAST outfmt 6 hits by resolving overlapping HSPs within each
(query, subject [, strand]) group using a greedy bitscore-first strategy.
Pre-filters on minimum length and minimum percent identity.

Expected input columns (17): qseqid sseqid staxid ssciname pident length
mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq
"""
import sys
import argparse
from dataclasses import dataclass
from typing import List, Iterator, Tuple, TextIO


@dataclass
class Hit:
  line: str
  sstart: int
  send: int
  bitscore: float

  @property
  def s_lo(self) -> int:
    return min(self.sstart, self.send)

  @property
  def s_hi(self) -> int:
    return max(self.sstart, self.send)

  @property
  def length(self) -> int:
    return self.s_hi - self.s_lo + 1


def parse_hit(line: str) -> Tuple[str, str, str, float, int, Hit]:
  fields = line.rstrip('\n').split('\t')
  qseqid = fields[0]
  sseqid = fields[1]
  pident = float(fields[4])
  aln_length = int(fields[5])
  sstart = int(fields[10])
  send = int(fields[11])
  bitscore = float(fields[13])
  sstrand = fields[14]
  return qseqid, sseqid, sstrand, pident, aln_length, Hit(line, sstart, send, bitscore)


def overlap_fraction_either(a: Hit, b: Hit) -> float:
  ov = min(a.s_hi, b.s_hi) - max(a.s_lo, b.s_lo) + 1
  if ov <= 0:
    return 0.0
  return ov / min(a.length, b.length)


def resolve_group(hits: List[Hit], overlap_thresh: float) -> List[Hit]:
  hits_sorted = sorted(hits, key=lambda h: h.bitscore, reverse=True)
  kept: List[Hit] = []
  for h in hits_sorted:
    conflict = False
    for k in kept:
      if overlap_fraction_either(h, k) >= overlap_thresh:
        conflict = True
        break
    if conflict:
      continue
    kept.append(h)
  return kept


def split_by_strand(
  hits_with_strand: List[Tuple[str, Hit]],
  strand_aware: bool,
) -> Iterator[List[Hit]]:
  if not strand_aware:
    yield [h for _, h in hits_with_strand]
    return
  plus_hits = [h for s, h in hits_with_strand if s == 'plus']
  minus_hits = [h for s, h in hits_with_strand if s == 'minus']
  if plus_hits:
    yield plus_hits
  if minus_hits:
    yield minus_hits


def group_stream(
  fh: TextIO,
  min_length: int,
  min_pident: float,
) -> Iterator[Tuple[str, str, List[Tuple[str, Hit]], int]]:
  """Yield (qseqid, sseqid, hits_with_strand, n_prefiltered_in_this_group)."""
  current_key = None
  buffer: List[Tuple[str, Hit]] = []
  n_prefiltered = 0
  for line in fh:
    if not line.strip():
      continue
    qseqid, sseqid, sstrand, pident, aln_length, hit = parse_hit(line)
    key = (qseqid, sseqid)
    if current_key is not None and key != current_key:
      yield current_key[0], current_key[1], buffer, n_prefiltered
      current_key = None
      buffer = []
      n_prefiltered = 0
    if aln_length < min_length or pident < min_pident:
      n_prefiltered += 1
      if current_key is None:
        current_key = key
      continue
    if current_key is None:
      current_key = key
      buffer = [(sstrand, hit)]
    else:
      buffer.append((sstrand, hit))
  if current_key is not None:
    yield current_key[0], current_key[1], buffer, n_prefiltered


def main():
  ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  ap.add_argument('input', help='BLAST outfmt 6 TSV file, or - for stdin')
  ap.add_argument('-o', '--output', default='-', help='output TSV')
  ap.add_argument(
    '-t', '--overlap-threshold', type=float, default=0.5,
    help='minimum overlap fraction (of shorter hit) at which the lower-'
         'bitscore hit is dropped; set to 1.0 to disable overlap resolution',
  )
  ap.add_argument(
    '--strand-agnostic', action='store_true',
    help='treat plus and minus strand hits as overlapping (default: treat '
         'them as independent, so a plus-strand hit at X cannot suppress a '
         'minus-strand hit at X)',
  )
  ap.add_argument(
    '--min-length', type=int, default=0,
    help='minimum alignment length (BLAST column 6) to retain a hit',
  )
  ap.add_argument(
    '--min-pident', type=float, default=70.0,
    help='minimum percent identity (BLAST column 5) to retain a hit',
  )
  ap.add_argument(
    '--stats', action='store_true',
    help='print counts to stderr',
  )
  args = ap.parse_args()

  fin = sys.stdin if args.input == '-' else open(args.input, 'r')
  fout = sys.stdout if args.output == '-' else open(args.output, 'w')

  n_groups = 0
  n_passed_prefilter = 0
  n_prefiltered_total = 0
  n_out = 0

  try:
    for qseqid, sseqid, hits_with_strand, n_prefiltered in group_stream(
      fin, args.min_length, args.min_pident,
    ):
      n_prefiltered_total += n_prefiltered
      if not hits_with_strand:
        continue
      n_groups += 1
      n_passed_prefilter += len(hits_with_strand)
      for hits in split_by_strand(hits_with_strand, not args.strand_agnostic):
        kept = resolve_group(hits, args.overlap_threshold)
        kept.sort(key=lambda h: h.s_lo)
        for h in kept:
          fout.write(h.line if h.line.endswith('\n') else h.line + '\n')
        n_out += len(kept)
      if n_groups % 100 == 0:
        n_in_so_far = n_passed_prefilter + n_prefiltered_total
        sys.stderr.write(
          f'\rprocessed {n_groups} groups, {n_in_so_far} lines'
        )
        sys.stderr.flush()
    if n_groups >= 100:
      sys.stderr.write('\n')
      sys.stderr.flush()
  finally:
    if fin is not sys.stdin:
      fin.close()
    if fout is not sys.stdout:
      fout.close()

  if args.stats:
    n_in_total = n_passed_prefilter + n_prefiltered_total
    dropped_overlap = n_passed_prefilter - n_out
    pct_pref = 100.0 * n_prefiltered_total / n_in_total if n_in_total else 0.0
    pct_ovl = 100.0 * dropped_overlap / n_passed_prefilter if n_passed_prefilter else 0.0
    pct_kept = 100.0 * n_out / n_in_total if n_in_total else 0.0
    sys.stderr.write(
      f'input_hits:            {n_in_total}\n'
      f'prefiltered_out:       {n_prefiltered_total} ({pct_pref:.1f}%)\n'
      f'passed_prefilter:      {n_passed_prefilter}\n'
      f'groups_processed:      {n_groups}\n'
      f'dropped_by_overlap:    {dropped_overlap} ({pct_ovl:.1f}% of prefiltered)\n'
      f'kept:                  {n_out} ({pct_kept:.1f}% of input)\n'
    )


if __name__ == '__main__':
  main()