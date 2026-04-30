#!/usr/bin/env python3
import argparse
import os
import sys
from collections import OrderedDict

def parse_fasta(path):
  seqs = OrderedDict()
  current = None
  buf = []
  with open(path) as f:
    for line in f:
      line = line.rstrip()
      if not line:
        continue
      if line.startswith('>'):
        if current is not None:
          seqs[current] = ''.join(buf).upper()
        current = line[1:].split()[0]
        buf = []
      else:
        buf.append(line)
    if current is not None:
      seqs[current] = ''.join(buf).upper()
  return seqs

def main():
  ap = argparse.ArgumentParser(
    description='Concatenate MSA FASTA files and produce a RAxML-style partition file for IQ-TREE.'
  )
  ap.add_argument('-i', '--input', nargs='+', required=True,
                  help='Input aligned FASTA files in desired concatenation order.')
  ap.add_argument('-o', '--out-fasta', required=True,
                  help='Output concatenated FASTA file.')
  ap.add_argument('-p', '--out-partition', required=True,
                  help='Output RAxML-style partition file.')
  ap.add_argument('-m', '--model', default='DNA',
                  help='Model label for partition file (e.g. DNA, AA, WAG, GTR). Default: DNA')
  ap.add_argument('--missing-char', default='-',
                  help='Character to fill for taxa missing from a partition. Default: -')
  ap.add_argument('--strict', action='store_true',
                  help='Fail if any taxon is missing from any alignment instead of padding with gaps.')
  args = ap.parse_args()

  alignments = []
  all_taxa = OrderedDict()
  for path in args.input:
    if not os.path.isfile(path):
      sys.exit(f'Error: file not found: {path}')
    seqs = parse_fasta(path)
    if not seqs:
      sys.exit(f'Error: no sequences in {path}')
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
      sys.exit(f'Error: sequences in {path} are not all the same length '
               f'(found lengths: {sorted(lengths)}). Is this a true alignment?')
    aln_len = lengths.pop()
    name = os.path.splitext(os.path.basename(path))[0]
    alignments.append((name, aln_len, seqs))
    for t in seqs:
      all_taxa[t] = None

  missing_report = []
  for name, _, seqs in alignments:
    missing = [t for t in all_taxa if t not in seqs]
    if missing:
      missing_report.append((name, missing))

  if missing_report:
    total_taxa = len(all_taxa)
    print('', file=sys.stderr)
    print('WARNING: missing taxa detected', file=sys.stderr)
    print('=' * 60, file=sys.stderr)
    for name, missing in missing_report:
      print(f'  [{name}] missing {len(missing)}/{total_taxa} taxa:', file=sys.stderr)
      for t in missing:
        print(f'    - {t}', file=sys.stderr)
    print('=' * 60, file=sys.stderr)
    if args.strict:
      sys.exit('Exiting due to --strict. Remove --strict to pad missing taxa with gaps.')
    else:
      print(f'Missing taxa will be padded with "{args.missing_char}" across their missing partitions.',
            file=sys.stderr)
      print('', file=sys.stderr)

  concat = OrderedDict((t, []) for t in all_taxa)
  for name, aln_len, seqs in alignments:
    gap = args.missing_char * aln_len
    for t in all_taxa:
      concat[t].append(seqs.get(t, gap))

  with open(args.out_fasta, 'w') as f:
    for t, parts in concat.items():
      f.write(f'>{t}\n')
      f.write(''.join(parts))
      f.write('\n')

  with open(args.out_partition, 'w') as f:
    start = 1
    for name, aln_len, _ in alignments:
      end = start + aln_len - 1
      f.write(f'{args.model}, {name} = {start}-{end}\n')
      start = end + 1

  total_len = sum(a[1] for a in alignments)
  print(f'Wrote {args.out_fasta}: {len(all_taxa)} taxa x {total_len} sites '
        f'from {len(alignments)} partitions', file=sys.stderr)
  print(f'Wrote {args.out_partition}', file=sys.stderr)

if __name__ == '__main__':
  main()