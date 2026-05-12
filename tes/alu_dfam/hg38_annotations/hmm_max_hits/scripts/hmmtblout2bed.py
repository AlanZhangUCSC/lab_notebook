import sys
import re

pattern = re.compile(r'^(.+?)::([^:]+):(\d+)-(\d+)\(([+-])\)$')

with open(sys.argv[1], 'r') as f:
  for line in f:
    fields = line.strip().split('\t')
    copy_name, bitscore = fields[2], fields[13]
    m = pattern.match(copy_name)
    if not m:
      print(f'Error: {copy_name} does not match the pattern', file=sys.stderr)
      exit(1)
    name, chrom, start, end, strand = m.groups()
    alignment_start, alignment_end = int(fields[6]), int(fields[7])
    if alignment_start > alignment_end:
      if strand == '+':
        strand = '-'
      elif strand == '-':
        strand = '+'
      else:
        raise ValueError(f'Invalid strand: {strand}')
        exit(1)
    print(f'{chrom}\t{start}\t{end}\t{name}\t{bitscore}\t{strand}')
