import sys
import os
import pickle
import argparse
from collections import defaultdict
from bisect import bisect_left, bisect_right

parser = argparse.ArgumentParser(description='Compare two aligned sequences.')
parser.add_argument('file1', help='First file containing the aligned sequence')
parser.add_argument('file2', nargs='?', help='Second file containing the aligned sequence (optional)', default=None)
parser.add_argument('--block_ranges', help='File containing the block ranges')
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

def compare_sequences(seq1, seq2):
  mismatch_count = 0
  mismatch_positions = []
  nuc_positions = []
  mismatched_nuc_positions = []
  for i, (a, b) in enumerate(zip(seq1, seq2)):
    if a != '-' or b != '-':
      nuc_positions.append(i)
      if a != b:
        mismatched_nuc_positions.append(i)
    if a != b:
      mismatch_count += 1
      mismatch_positions.append(i)

  return mismatch_count, mismatch_positions, nuc_positions, mismatched_nuc_positions

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
all_mismatch_count, all_mismatch_positions, all_nuc_positions, all_mismatched_nuc_positions = compare_sequences(seq1, seq2)

all_mismatch_ranges = merge_ranges(all_mismatch_positions)

class FastIntervalMap:
  def __init__(self):
    self.starts = []
    self.ends = []    
    self.ids = []
      
  def add_interval(self, id, start, end):
    pos = bisect_left(self.starts, start)
    self.starts.insert(pos, start)
    self.ends.insert(pos, end)
    self.ids.insert(pos, id)
      
  def find_overlapping(self, start, end):
    left = bisect_right(self.ends, start - 1)
    right = bisect_left(self.starts, end + 1)
    return sorted([(self.ids[i], self.starts[i], self.ends[i], self.ends[i] - self.starts[i] + 1) 
                  for i in range(left, right)],
                  key=lambda x: (x[0], x[1]))

block_ranges_file = args.block_ranges
block_ranges_file_dir = os.path.dirname(block_ranges_file)
block_ranges_file_name = os.path.basename(block_ranges_file)
block_ranges_file_prefix = block_ranges_file_name.split('.')[0]
block_ranges_pickle_file = os.path.join(block_ranges_file_dir, f'{block_ranges_file_prefix}.pickle')

block_ranges = FastIntervalMap()
block_ranges_by_id = {}
# if os.path.exists(block_ranges_pickle_file):
#   print(f'Loading block ranges from {block_ranges_pickle_file}', file=sys.stderr)
#   with open(block_ranges_pickle_file, 'rb') as fh:
#     block_ranges = pickle.load(fh)
# else:
# print(f'Loading block ranges from {block_ranges_file}', file=sys.stderr)
with open(block_ranges_file) as fh:
  for i, line in enumerate(fh):
    start, end = map(int, line.strip().split('\t'))
    block_ranges.add_interval(i, start, end)
    block_ranges_by_id[i] = (start, end)
with open(block_ranges_pickle_file, 'wb') as fh:
  pickle.dump(block_ranges, fh)

gap_to_nuc_block_groups = []
seen_blocks = set()
corresponding_block_ids = []
for cur_range in all_mismatch_ranges:
  overlapping_blocks = block_ranges.find_overlapping(cur_range[0], cur_range[1])
  cur_corresponding_block_ids = set()
  for overlapping_block in overlapping_blocks:
    cur_corresponding_block_ids.add(overlapping_block[0])
  corresponding_block_ids.append(cur_corresponding_block_ids)
  cur_gap_to_nuc_block_group = []
  for overlapping_block in overlapping_blocks:
    block_id, start, end, length = overlapping_block
    if block_id in seen_blocks: continue
    seen_blocks.add(block_id)
    
    block_seq1 = seq1[start:end+1]
    block_seq2 = seq2[start:end+1]
    block1_is_gap = all(char == '-' for char in block_seq1)
    block2_is_gap = all(char == '-' for char in block_seq2)

    assert(not (block1_is_gap and block2_is_gap))
    if block1_is_gap and not block2_is_gap:
      cur_gap_to_nuc_block_group.append((overlapping_block, block_seq2, True))
    elif block2_is_gap and not block1_is_gap:
      cur_gap_to_nuc_block_group.append((overlapping_block, block_seq1, False))

  if cur_gap_to_nuc_block_group:
    gap_to_nuc_block_groups.append(cur_gap_to_nuc_block_group)

potential_ma_blocks_by_length = defaultdict(list)
for gap_to_nuc_block_group in gap_to_nuc_block_groups:
  for gap_to_nuc_block in gap_to_nuc_block_group:
    block_id, start, end, length = gap_to_nuc_block[0]
    block_seq = gap_to_nuc_block[1]
    is_seq2 = gap_to_nuc_block[2]
    potential_ma_blocks_by_length[length].append((block_id, start, end, block_seq, is_seq2))

ma_blocks = set()
ma_block_pairs = []
for length, blocks in potential_ma_blocks_by_length.items():
  if len(blocks) < 2: continue
  paired_blocks = []
  for i in range(len(blocks)):
    block1_id, block1_start, block1_end, block1_seq, block1_is_seq2 = blocks[i]
    for j in range(i+1, len(blocks)):
      block2_id, block2_start, block2_end, block2_seq, block2_is_seq2 = blocks[j]
      if block1_is_seq2 == block2_is_seq2: continue
      mismatch_count, mismatch_positions, nuc_positions, mismatched_nuc_positions = compare_sequences(block1_seq, block2_seq)
      similarity = 1 - (len(mismatched_nuc_positions) / len(nuc_positions))
      if similarity > 0.8:
        paired_blocks.append((block1_id, block2_id, mismatch_count, len(nuc_positions), len(mismatched_nuc_positions)))

  paired_blocks.sort(key=lambda x: x[4])

  seen_blocks = set()
  for paired_block in paired_blocks:
    block1_id, block2_id, mismatch_count, num_nuc_pairs, num_mismatched_nuc_pairs = paired_block
    if block1_id in seen_blocks or block2_id in seen_blocks: continue
    seen_blocks.add(block1_id)
    seen_blocks.add(block2_id)
    ma_block_pairs.append(paired_block)
    ma_blocks.add(block1_id)
    ma_blocks.add(block2_id)

mismatch_count_in_ma_blocks = 0
test = 0
for i,cur_range in enumerate(all_mismatch_ranges):
  cur_beg, cur_end = cur_range
  cur_len = cur_end - cur_beg + 1
  overlapping_block_ids = corresponding_block_ids[i]
  if len(overlapping_block_ids) == 1:
    overlapping_block_id = list(overlapping_block_ids)[0]
    if overlapping_block_id in ma_blocks:
      mismatch_count_in_ma_blocks += cur_len
  elif len(overlapping_block_ids) > 1:
    for cur_block_id in overlapping_block_ids:
      if cur_block_id in ma_blocks:
        cur_block_beg, cur_block_end = block_ranges_by_id[cur_block_id]
        for j in range(cur_block_beg, cur_block_end+1):
          if j >= cur_beg and j <= cur_end:
            mismatch_count_in_ma_blocks += 1
  else:
    raise ValueError(f'No overlapping blocks found for {cur_range}')

mismatch_count_in_corrected_blocks = 0

for ma_block_pair in ma_block_pairs:
  block1_id, block2_id, mismatch_count, num_nuc_pairs, num_mismatched_nuc_pairs = ma_block_pair
  mismatch_count_in_corrected_blocks += mismatch_count



print(f'#mismatch_count_in_ma_blocks: {mismatch_count_in_ma_blocks}')
print(f'#mismatch_count_in_corrected_blocks: {mismatch_count_in_corrected_blocks}')
print(f'#all_mismatch_count: {all_mismatch_count}')
print(f'#expected_mismatch_count_after_correction: {all_mismatch_count - mismatch_count_in_ma_blocks + mismatch_count_in_corrected_blocks}')
print(f'#block1\tblock2\tmismatch_count\tnum_nuc_pairs')
for ma_block_pair in ma_block_pairs:
  block1_id, block2_id, mismatch_count, num_nuc_pairs, num_mismatched_nuc_pairs = ma_block_pair
  print(f'{block1_id}\t{block2_id}\t{mismatch_count}\t{num_nuc_pairs}')







