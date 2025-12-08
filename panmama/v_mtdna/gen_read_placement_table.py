# Given a <mgsr.abundance.out> file and a <mgsr.assignedReads.out> file, output a tsv of the reads that parsimoniously
# map to each node with format

# nodes r1  r2  r3  r4
# node1 1   1   1   0
# node2 1   0   0   0
# ...

# where 1 indicates the read places to that node parsimoniously, and 0 indicates it does not.

import os
import sys
import pandas as pd
import subprocess
from collections import defaultdict

abundance_path = sys.argv[1]
assigned_reads_path = sys.argv[2]
fastq_path = sys.argv[3]
output = sys.argv[4]


fastq_filtered = f'{fastq_path}.filtered.fq'
filter_cmd = f'~/tools/BBTools/bbduk.sh in={fastq_path} out={fastq_filtered} entropy=0.7'
subprocess.run(filter_cmd, shell=True, check=True)

fastq_filtered_sorted = f'{fastq_path}.filtered.sorted.fq'
sort_cmd = f"cat {fastq_filtered} | paste - - - - | sort -k1,1 -S 3G | tr '\\t' '\\n' > {fastq_filtered_sorted}"
subprocess.run(sort_cmd, shell=True, check=True)

fastq_path = fastq_filtered_sorted

representative_to_group = defaultdict(list)
with open(abundance_path, 'r') as fh:
  for line in fh.readlines():
    line = line.strip()
    haplotypes = line.split()[0].split(',')
    for hap in haplotypes:
      representative_to_group[haplotypes[0]].append(hap)

all_reads = set()
node_to_reads = defaultdict(set)
with open(assigned_reads_path, 'r') as fh:
  for line in fh.readlines():
    fields = line.strip().split()
    representative = fields[0]
    reads = list(map(int,fields[-1].split(',')))
    all_reads.update(reads)
    for hap in representative_to_group[representative]:
      node_to_reads[hap].update(reads)

read_ids = []
with open(fastq_path, 'r') as fh:
  for i, line in enumerate(fh.readlines()):
    if i % 4 == 0:
      read_id = line.strip().split()[0][1:]
      read_ids.append(read_id)

headers = ['SampleId'] + [read_ids[index] for index in sorted(list(all_reads))]
headers = list(map(str, headers))

data = []
for node in node_to_reads:
  reads = node_to_reads[node]
  row = [node]
  for i, read in enumerate(sorted(list(all_reads))):
    if int(read) in reads:
      row.append('1')
    else:
      row.append('0')
  data.append(row)

df = pd.DataFrame(data, columns=headers)

df = df.transpose().reset_index()
df.to_csv(output, sep='\t', index=False, header=False)

os.remove(fastq_filtered)
os.remove(fastq_filtered_sorted)
