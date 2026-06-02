import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Gather assignedReads outputs")
parser.add_argument('-i', '--input', nargs='+', required=True, help='Input file(s), 1 or more')
parser.add_argument('-o', '--output', required=True, help='Output file path')
args = parser.parse_args()

input_files = args.input

node_read_counts = defaultdict(int)
for input_file in input_files:
  with open(input_file, 'r') as f:
    for line in f:
      fields = line.strip().split('\t')
      nodes_id, read_count = fields[0], int(fields[2])
      node_read_counts[nodes_id] += read_count

with open(args.output, 'w') as f:
  for nodes_id, read_count in node_read_counts.items():
    f.write(f"{nodes_id}\t.\t{read_count}\t.\n")
