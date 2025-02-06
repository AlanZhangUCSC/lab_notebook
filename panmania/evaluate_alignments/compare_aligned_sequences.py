import sys

import argparse

parser = argparse.ArgumentParser(description='Compare two aligned sequences.')
parser.add_argument('file1', help='First file containing the aligned sequence')
parser.add_argument('file2', nargs='?', help='Second file containing the aligned sequence (optional)', default=None)
parser.add_argument('--print_ranges', action='store_true', help='Print the ranges of mismatches')
args = parser.parse_args()


def read_two_sequences(fasta_file):
    sequences = []
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:  # If we have collected sequence data
                    sequences.append(''.join(current_seq))
                    current_seq = []  # Reset for next sequence
                    if len(sequences) >= 2:  # We already have two sequences
                        break
            else:
                current_seq.append(line)
        
        # Don't forget to append the last sequence
        if current_seq:
            sequences.append(''.join(current_seq))
    
    if len(sequences) < 2:
        raise ValueError("File must contain at least two sequences")
    
    return sequences[0], sequences[1]

def merge_ranges(numbers):
  if not numbers:
    return []
  
  numbers = sorted(numbers)
  ranges = []
  start = numbers[0]
  prev = numbers[0]
  
  for num in numbers[1:]:
    if num > prev + 1:
      ranges.append((start, prev))
      start = num
    prev = num
  
  ranges.append((start, prev))
  
  return ranges

if args.file2:
  file1=args.file1
  file2=args.file2
  with open(file1) as fh:
    next(fh)
    seq1 = fh.read().strip().replace('\n', '').upper()

  with open(file2) as fh:
    next(fh)
    seq2 = fh.read().strip().replace('\n', '').upper()
else:
  seq1, seq2 = read_two_sequences(args.file1)

assert(len(seq1) == len(seq2))
mismatch_positions = []
mismatch_count = 0
for i, (a, b) in enumerate(zip(seq1, seq2)):
  if a != b:
    mismatch_count += 1
    if args.print_ranges:
      mismatch_positions.append(i)

mismatch_ranges = merge_ranges(mismatch_positions)

print(mismatch_count)
if args.print_ranges:
  for cur_range in mismatch_ranges:
    print(f'{cur_range[1] - cur_range[0] + 1}\t{cur_range[0]}-{cur_range[1]}')
