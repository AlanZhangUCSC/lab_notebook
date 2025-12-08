import os
import sys
import pandas as pd
from ete3 import Tree

node_conversion_file = 'final_analysis/node_conversion.txt'
pathphynder_to_panman_node = {}
panman_to_pathphynder_node = {}
with open(node_conversion_file, 'r') as fh:
  node_map = {}
  for line in fh.readlines():
    fields = line.strip().split()
    pathphynder_node = fields[0]
    panman_node = fields[1]
    pathphynder_to_panman_node[pathphynder_node] = panman_node
    panman_to_pathphynder_node[panman_node] = pathphynder_node

panmama_out_dir = 'ancient_mammoth_mito_out'
pathphynder_out_dir = 'final_analysis/data/pathphynder/results_folder'

panmama_out_files = {}
pathphynder_out_files = {}

samples_prefix = set()

print('Loading Tree...')
tree = Tree('wooly_mammoth.panman.newick', format=1)

print(f'Loading panmama outputs from {panmama_out_dir}')
for filename in os.listdir(panmama_out_dir):
  if filename.endswith('.read_scores.mgsr.abundance.out'):
    cur_sample_prefix = filename.split('.')[0]
    samples_prefix.add(cur_sample_prefix)
    panmama_out_files[cur_sample_prefix] = os.path.join(panmama_out_dir, filename)

print(f'Loading pathphynder outputs from {pathphynder_out_dir}')
for filename in os.listdir(pathphynder_out_dir):
  if filename.endswith('.best_node_info.txt'):
    cur_sample_prefix = filename.split('.')[0]
    samples_prefix.add(cur_sample_prefix)
    pathphynder_out_files[cur_sample_prefix] = os.path.join(pathphynder_out_dir, filename)

print(f'Loaded {len(panmama_out_files)} panmama output files and {len(pathphynder_out_files)} pathphynder output files.')

data = []
headers = ['sample', 'panmama_mixed', 'panmama_node', 'pathphynder_node', 'distance']
for sample in samples_prefix:
  panmama_mixed = False
  panmama_node = 'NA'
  pathphynder_node = 'NA'
  distance = 'NA'

  print(f'Processing sample {sample}...', file=sys.stderr)
  if sample in panmama_out_files:
    nodes = []
    print(panmama_out_files[sample])
    with open(panmama_out_files[sample], 'r') as f:
      for line in f:
        nodes.append(line.strip().split()[0].split(','))
    if panmama_mixed:
      panmama_mixed = True
    
    nodes = nodes[0]
    lca = tree.get_common_ancestor(nodes)

    panmama_node = lca.name

  if sample in pathphynder_out_files:  
    with open(pathphynder_out_files[sample], 'r') as f:
      pathphynder_node = 'Node' + f.readlines()[1].strip().split()[1]
      pathphynder_node = pathphynder_to_panman_node[pathphynder_node]
  
  data.append([sample, panmama_mixed, panmama_node, pathphynder_node])
  if (panmama_node != 'NA') and (pathphynder_node != 'NA'):
    print(f'Calculating distance between {panmama_node} and {pathphynder_node}...', file=sys.stderr)
    distance = lca.get_distance(pathphynder_node)
    data[-1].append(distance)
  else:
    data[-1].append('NA')

df = pd.DataFrame(data, columns=headers)
df.to_csv('compare_panmama_pathphynder_mammoth.tsv', sep='\t', index=False)



