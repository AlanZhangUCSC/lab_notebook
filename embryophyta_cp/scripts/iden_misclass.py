import argparse
from ete3 import Tree
from collections import defaultdict
from misc import sister_leaf_names, closest_relative_names

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree', required=True, help='input newick file')
parser.add_argument('-m', '--metadata', required=True, help='input metadata file')
args = parser.parse_args()

samples = []
sample_to_order = {}
sample_to_class = {}
orders = defaultdict(list)
classes = defaultdict(list)

with open(args.metadata, 'r') as f:
  for line in f:
    fields = line.strip().split('\t')
    sample = fields[0]
    sample_class = fields[5]
    sample_order = fields[6]
    samples.append(sample)
    sample_to_order[sample] = sample_order
    sample_to_class[sample] = sample_class
    orders[sample_order].append(sample)
    classes[sample_class].append(sample)

tree = Tree(args.tree)

for leaf in tree.get_leaves():
  sample = leaf.name
  if sample not in sample_to_order:
    raise ValueError(f"Sample {sample} not found in metadata file")
    exit(1)
  order = sample_to_order[sample]
  members = orders[order]
  num_members = len(members)

  if num_members == 1:
    continue
  elif num_members == 2:
    other = members[0] if members[1] == sample else members[1]
    sisters = sister_leaf_names(leaf)
    if other not in sisters:
      print(f'{sample}\t{order}\t{num_members}\t{",".join(members)}\t{",".join(sisters)}')
  elif num_members >= 3:
    relatives = closest_relative_names(leaf)
    if len(relatives) < 2:
      raise ValueError(f"Invalid number of relatives: {len(relatives)}")
      exit(1)
    relatives_orders = {sample_to_order[r] for r in relatives}
    if len(relatives_orders) == 1:
      rel_order = next(iter(relatives_orders))
      if rel_order != order:
        print(f'{sample}\t{order}\t{num_members}\t{",".join(members)}\t{rel_order}\t{",".join(relatives)}')
    else:
      continue
  else:
    raise ValueError(f"Invalid number of members: {num_members}")
    exit(1)