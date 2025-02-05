import sys
import random
import numpy as np
fasta_file = sys.argv[1]


random.seed(69)
np.random.seed(69)

def mutate_seq(seq, num_mutations):
  mutate_positions = np.random.choice(len(seq), num_mutations, replace=False)
  for pos in mutate_positions:
    current_nuc = seq[pos]
    new_nuc = random.choice([nuc for nuc in 'ACGT' if nuc != current_nuc])
    seq[pos] = new_nuc
  return seq

seq = []
with open(fasta_file, 'r') as f:
  for line in f:
    if line.startswith('>'): continue
    seq += list(line.strip())


for i in range(0, len(seq), 15):
  substring = seq[i:i+150]
  if len(substring) < 150: break
  for j in range(0, 50):
    substring_copy = substring.copy()
    mutated_substring = mutate_seq(substring_copy, j)
    seq_id = f'read_{i}_{j}'
    print(f'@{seq_id}')
    print(''.join(mutated_substring))
    print('+')
    print(''.join(['I'] * len(mutated_substring)))
