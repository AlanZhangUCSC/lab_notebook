import sys
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Compare two fastags files.')
parser.add_argument('fastags_file_1', help='First fastags file')
parser.add_argument('fastags_file_2', help='Second fastags file')
args = parser.parse_args()

fastags_file_1 = args.fastags_file_1
fastags_file_2 = args.fastags_file_2


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

blocks_1 = []
blocks_2 = []

with open(fastags_file_1, 'r') as f:
  next(f)
  for line in f:
    line = line.strip()
    if line == '': continue
    blocks_1.append(line)

with open(fastags_file_2, 'r') as f:
  next(f)
  for line in f:
    line = line.strip()
    if line == '': continue
    blocks_2.append(line)

assert len(blocks_1) == len(blocks_2)

diffs_by_block = defaultdict(int)
potential_ma_blocks_1 = defaultdict(list)
potential_ma_blocks_2 = defaultdict(list)
for i, (block_1, block_2) in enumerate(zip(blocks_1, blocks_2)):
  block_1_on = block_1[0] != '.'
  block_2_on = block_2[0] != '.'

  if block_1_on and block_2_on:
    # both blocks are on, check diffs
    for c1, c2 in zip(block_1, block_2):
      if c1 != c2:
        diffs_by_block[i] += 1
  elif block_1_on and not block_2_on:
    # block 1 is on and block 2 is off, num_diff == num_nucs
    potential_ma_blocks_1[len(block_1)].append(i)
    for c in block_1:
      if c != '-':
        diffs_by_block[i] += 1
  elif not block_1_on and block_2_on:
    # block 1 is off and block 2 is on, num_diff == num_nucs
    potential_ma_blocks_2[len(block_2)].append(i)
    for c in block_2:
      if c != '-':
        diffs_by_block[i] += 1
  else:
    # both blocks are off, do nothing
    continue

ma_blocks = set()
ma_block_pairs = []
for block_len in sorted(potential_ma_blocks_1.keys()):
  paired_blocks = []
  block_1_ids = potential_ma_blocks_1[block_len]
  if block_len in potential_ma_blocks_2:
    block_2_ids = potential_ma_blocks_2[block_len]
    for block_1_id in block_1_ids:
      for block_2_id in block_2_ids:
        assert(block_1_id != block_2_id)
        mismatch_count, mismatch_positions, nuc_positions, mismatched_nuc_positions = compare_sequences(blocks_1[block_1_id], blocks_2[block_2_id])
        similarity = 1 - (len(mismatched_nuc_positions) / len(nuc_positions))
        if similarity > 0.75:
          paired_blocks.append((block_1_id, block_2_id, mismatch_count, len(nuc_positions), len(mismatched_nuc_positions)))

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

diffs_in_corrected_blocks = 0
for ma_block_pair in ma_block_pairs:
  block1_id, block2_id, mismatch_count, num_nuc_pairs, num_mismatched_nuc_pairs = ma_block_pair
  diffs_in_corrected_blocks += mismatch_count

total_diffs = 0
corrected_diffs = 0

for block_id, diffs in diffs_by_block.items():
  total_diffs += diffs
  if block_id not in ma_blocks:
    corrected_diffs += diffs

print(f'#total_diffs: {total_diffs}')
print(f'#corrected_diffs: {corrected_diffs}')
print(f'#diffs_in_corrected_blocks: {diffs_in_corrected_blocks}')

print(f'#fastags_file_1:{fastags_file_1}')
print(f'#fastags_file_2:{fastags_file_2}')
print('#ma_block_1\tma_block_2\tnum_diffs\tnum_nuc_pairs\tnum_mismatched_nuc_pairs')
for ma_block_pair in ma_block_pairs:
  block1_id, block2_id, mismatch_count, num_nuc_pairs, num_mismatched_nuc_pairs = ma_block_pair
  print(f'{block1_id}\t{block2_id}\t{mismatch_count}\t{num_nuc_pairs}\t{num_mismatched_nuc_pairs}')



