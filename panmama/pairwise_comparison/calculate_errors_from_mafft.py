import sys
from read_fasta import read_record

seqs = []
for idn, seq_str in read_record('-'):
  seqs.append((idn, seq_str))

idn1, seq1 = seqs[0]
idn2, seq2 = seqs[1]

print(f"Sequence lengths are {len(seq1)} and {len(seq2)}", file=sys.stderr)

match_string = ''

canonical = {
  'a': True, 't': True, 'c': True, 'g': True,
  'A': True, 'T': True, 'C': True, 'G': True
}

num_errors = 0
num_snps = 0
num_gaps = 0
num_gaps_from_ends = 0
num_base_to_ambiguous = 0
for i in range(len(seq1)):
  if seq1[i] != seq2[i]:
    match_string += ':'
    if seq1[i] == '-' or seq2[i] == '-':
      if seq1[i] in canonical or seq2[i] in canonical:
        num_errors += 1
        num_gaps += 1
    else:
      if seq1[i] in canonical and seq2[i] in canonical:
        num_snps += 1
        num_errors += 1
      else:
        num_base_to_ambiguous += 1
  else:
    match_string += '|'
  if i % 1000 == 0:
    print(f'\r{i} positions processed', end='', file=sys.stderr)

for i in range(len(seq1)):
  if seq1[i] == '-' or seq2[i] == '-':
    if seq1[i] in canonical or seq2[i] in canonical:
      num_gaps_from_ends += 1
  else:
    break
for i in range(len(seq1)-1, -1, -1):
  if seq1[i] == '-' or seq2[i] == '-':
    if seq1[i] in canonical or seq2[i] in canonical:
      num_gaps_from_ends += 1
  else:
    break

line_size = 80
for i in range(0, len(match_string), line_size):
    print(seq1[i:i+line_size])
    print(match_string[i:i+line_size])
    print(seq2[i:i+line_size])
    print()


print(f'{num_errors} errors out of {len(seq1)} positions')
print(f'{num_snps} snps')
print(f'{num_base_to_ambiguous} ambiguous snps')
print(f'{num_gaps} gaps')
print(f'{num_gaps - num_gaps_from_ends} gaps edge corrected')