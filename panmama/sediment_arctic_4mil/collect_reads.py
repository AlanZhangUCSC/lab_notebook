from collections import defaultdict
import argparse
import os

parser = argparse.ArgumentParser(
  description='Split MGSR assigned reads into per-taxon FASTQ files.'
)
parser.add_argument('input_dir', help='Directory containing *.mgsr.assignedReadsLCANode.annotated.out and matching *.mgsr.assignedReads.fastq files')
parser.add_argument('target_taxon', help='Taxon name to extract, or "all" to split by every taxon')
parser.add_argument('out_dir', help='Output directory for per-taxon FASTQ files')
args = parser.parse_args()

input_dir = args.input_dir
target_taxon = args.target_taxon
out_dir = args.out_dir

os.makedirs(out_dir, exist_ok=True)

out_files = {}

def get_out_handle(taxon):
  if taxon not in out_files:
    out_files[taxon] = open(os.path.join(out_dir, f'{taxon}.fastq'), 'w')
  return out_files[taxon]

for file in os.listdir(input_dir):
  if not file.endswith('.mgsr.assignedReadsLCANode.annotated.out'):
    continue

  prefix = '.'.join(file.split('.')[:-4])
  fastq_path = os.path.join(input_dir, f'{prefix}.mgsr.assignedReads.fastq')
  ann_path = os.path.join(input_dir, file)

  index_to_taxon = {}
  with open(ann_path, 'r') as f:
    for line in f:
      parts = line.rstrip('\n').split('\t')
      if len(parts) < 4:
        continue
      taxon = parts[1]
      if target_taxon != 'all' and taxon != target_taxon:
        continue
      for i in parts[3].split(','):
        index_to_taxon[int(i)] = taxon

  if not index_to_taxon:
    continue

  with open(fastq_path, 'r') as f:
    idx = 0
    while True:
      record = [f.readline() for _ in range(4)]
      if not record[0]:
        break
      taxon = index_to_taxon.get(idx)
      if taxon is not None:
        get_out_handle(taxon).writelines(record)
      idx += 1

for out_file in out_files.values():
  out_file.close()