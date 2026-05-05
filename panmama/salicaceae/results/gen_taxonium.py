import argparse
import dendropy
import json
import math

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--nwk', type=str, required=True, help='Newick tree file')
parser.add_argument('-n', '--node_read_counts_file', type=str, required=True, help='Node read counts file')
parser.add_argument('-l', '--node_read_counts_lca_file', type=str, required=True, help='Node read counts LCA file')
parser.add_argument('-e', '--template_file', type=str, required=True, help='Template file')
parser.add_argument('-g', '--gamma', type=float, default=0.5, help='Gamma value for power scale')
args = parser.parse_args()

nwk_string = ''
with open(args.nwk, 'r') as f:
  nwk_string = f.read().strip()

if args.node_read_counts_file == '' and args.node_read_counts_lca_file == '':
  print("Error: Either node read counts file or node read counts LCA file must be provided")
  exit(1)

tree = dendropy.Tree.get(data=nwk_string, schema="newick")
leaf_nodes = [node for node in tree.leaf_node_iter()]

node_read_counts = {}
max_read_count = 0
min_read_count = float('inf')
if args.node_read_counts_file != '':
  with open(args.node_read_counts_file, 'r') as f:
    for line in f:
      node_label, count, _ = line.strip().split('\t')
      for node_label_clean in node_label.split(','):
        node_label_clean = node_label_clean.replace('_', ' ')
        node_read_counts[node_label_clean] = int(count)
        max_read_count = max(max_read_count, int(count))
        min_read_count = min(min_read_count, int(count))

node_read_counts_lca = {}
if args.node_read_counts_lca_file != '':
  with open(args.node_read_counts_lca_file, 'r') as f:
    for line in f:
      node_label, count_lca, _ = line.strip().split('\t')
      for node_label_clean in node_label.split(','):
        node_label_clean = node_label_clean.replace('_', ' ')
        count_lca = int(count_lca)
        node_read_counts_lca[node_label_clean] = count_lca

tree_labels = set()
for node in tree.postorder_node_iter():
  tree_labels.add(node.taxon.label if node.taxon else node.label)

lca_labels = set(node_read_counts_lca.keys())

node_lca_count = {}
min_node_lca_count = float('inf')
max_node_lca_count = 0
for node in tree.postorder_node_iter():
  node_label = node.taxon.label if node.taxon else node.label
  node_lca_count[node_label] = node_read_counts_lca.get(node_label, 0)
  min_node_lca_count = min(min_node_lca_count, node_lca_count[node_label])
  max_node_lca_count = max(max_node_lca_count, node_lca_count[node_label])

node_subtree_lca_count = {}
min_node_subtree_lca_count = float('inf')
max_node_subtree_lca_count = 0
for node in tree.postorder_node_iter():
  node_label = node.taxon.label if node.taxon else node.label
  node_subtree_lca_count[node_label] = node_lca_count.get(node_label, 0)
  for child in node.child_nodes():
    child_label = child.taxon.label if child.taxon else child.label
    node_subtree_lca_count[node_label] += node_subtree_lca_count.get(child_label, 0)
  min_node_subtree_lca_count = min(min_node_subtree_lca_count, node_subtree_lca_count[node_label])
  max_node_subtree_lca_count = max(max_node_subtree_lca_count, node_subtree_lca_count[node_label])


def power_scale_breakpoints(max_val, n_breaks, gamma=0.5):
  if max_val <= 0:
    return [0] * n_breaks
  return [int(max_val * (i / (n_breaks - 1)) ** (1 / gamma)) for i in range(n_breaks)]

with open(args.template_file, "r") as f:
  for line in f:
    obj = json.loads(line)
    if "version" in obj:
      obj["config"]["colorRamps"] = {}
      color_ramps = obj["config"]["colorRamps"]

      base_scale = [[0,"#F5F5F5"],[0,"#F0B5B5"],[0,"#EC7676"],[0,"#E50000"],[0,"#B30000"]]
      color_ramps["read_counts"] = {"scale": [list(s) for s in base_scale]}
      color_ramps["lca_counts"] = {"scale": [list(s) for s in base_scale]}
      color_ramps["lca_subtree_counts"] = {"scale": [list(s) for s in base_scale]}

      for ramp_name, max_val in [
        ("read_counts", max_read_count),
        ("lca_counts", max_node_lca_count),
        ("lca_subtree_counts", max_node_subtree_lca_count),
      ]:
        breaks = power_scale_breakpoints(max_val, len(color_ramps[ramp_name]["scale"]))
        for i, color_scale in enumerate(color_ramps[ramp_name]["scale"]):
          color_scale[0] = breaks[i]
    elif "name" in obj:
      sample_ids = obj["name"]
      if sample_ids.startswith('node_'):
        sample_ids = '_'.join(sample_ids.split('_')[:-1])
      sample_ids = sample_ids.replace('_', ' ')
      obj['read_counts'] = node_read_counts.get(sample_ids, 0)
      obj['lca_counts'] = node_lca_count.get(sample_ids, 0)
      obj['lca_subtree_counts'] = node_subtree_lca_count.get(sample_ids, 0)
    print(json.dumps(obj))