import argparse
import sys

parser = argparse.ArgumentParser(description="Annotate reads with metadata.")
parser.add_argument("reads_file", help="Input file with reads")
parser.add_argument("meta_file", help="Metadata file")
parser.add_argument("column", nargs='?', default="Families", help="Metadata column to annotate (default: 'Families')")
args = parser.parse_args()

reads_file = args.reads_file
meta_file = args.meta_file
column = args.column

meta_dict = {}
with open(meta_file, 'r') as f:
  header = f.readline().strip().split('\t')
  col_index = -1
  for i, col in enumerate(header):
    if col == column:
      col_index = i
      break
  if col_index == -1:
    raise ValueError(f"Column {column} not found in meta file")
  for line in f:
    fields = line.strip().split('\t')
    if len(fields) >= col_index + 1:
      meta_dict[fields[0]] = fields[col_index]
    else:
      meta_dict[fields[0]] = '.'


with open(reads_file, 'r') as f:
  for line in f:
    nodes, num_reads, reads = line.strip().split('\t')

    meta_info = set()
    for node in nodes.split(','):
      meta_info.add(meta_dict[node])

    print(nodes, ','.join(meta_info), num_reads, reads, sep='\t')