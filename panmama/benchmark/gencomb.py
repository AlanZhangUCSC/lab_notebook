import argparse
import itertools

#!/usr/bin/env python

parser = argparse.ArgumentParser(description='Process some inputs.')
parser.add_argument('--snps', nargs='+', help='List of SNPs')
parser.add_argument('--haplotypes', nargs='+', help='List of haplotypes')
parser.add_argument('--percent-mutated', nargs='+', type=float, help='List of percent mutated values')
parser.add_argument('--seq-types', nargs='+', help='List of sequence types')
parser.add_argument('--num-reads', nargs='+', type=int, help='List of number of reads')
parser.add_argument('--num-rep', type=int, help='Number of replicates')

args = parser.parse_args()

snps = list(map(int, args.snps)) if args.snps else []
haplotypes = list(map(int, args.haplotypes)) if args.haplotypes else []
percent_mutated = list(map(float, args.percent_mutated)) if args.percent_mutated else []
seq_types = list(args.seq_types) if args.seq_types else []
num_reads = list(map(int, args.num_reads)) if args.num_reads else []


# Validate percent_mutated
if any(pm == 0 for pm in percent_mutated):
  raise ValueError("percent-mutated cannot have 0.")

# Generate combinations
combinations = []
for seq_type in seq_types:
  for haplotype in haplotypes:
    for snp in snps:
      if snp == 0:
        for num_read in num_reads:
          for rep in range(args.num_rep):
            combinations.append((seq_type, haplotype, snp, 0, num_read, rep))  # Ignore percent-mutated 
      else:
        for pm in percent_mutated:
          for num_read in num_reads:
            if (pm * haplotype).is_integer():
              for rep in range(args.num_rep):
                combinations.append((seq_type, haplotype, snp, pm, num_read, rep))
              

# Output the combinations
print("seq_type\thaplotype\tsnp\tpercent_mutated\tnum_reads\trep")
for combo in combinations:
  print('\t'.join(map(str, combo)))