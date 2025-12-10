import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
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
    lca = tree.get_common_ancestor(nodes) if len(nodes) > 1 else tree.search_nodes(name=nodes[0])[0]
    panmama_node = lca.name if len(nodes) > 1 else nodes[0]
  
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

print('\nGenerating plot...')
df_filtered = df[df['distance'] != 'NA'].copy()
df_filtered['distance'] = pd.to_numeric(df_filtered['distance'])
distances = df_filtered['distance'].values

fig, ax = plt.subplots(figsize=(10, 6))

unique_distances, counts = np.unique(distances, return_counts=True)
distance_dict = dict(zip(unique_distances, counts))

all_distances = np.arange(int(np.min(distances)), int(np.max(distances)) + 1)
all_counts = [distance_dict.get(d, 0) for d in all_distances]

bars = ax.bar(all_distances, all_counts, edgecolor='black', alpha=0.7, width=0.8, color='steelblue')

for bar, count in zip(bars, all_counts):
  if count > 0:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2., height,
            f'{int(count)}',
            ha='center', va='bottom', fontsize=10)

ax.set_xlabel('Distance', fontsize=12)
ax.set_ylabel('Count', fontsize=12)
ax.set_title('Distance Counts Between Pathphynder and Panmama Nodes', fontsize=14)
ax.grid(axis='y', alpha=0.3)
ax.set_xticks(all_distances)

plt.tight_layout()
plt.savefig('compare_panmama_pathphynder_mammoth.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"\nTotal samples: {len(df)}")
print(f"Samples with valid distance: {len(df_filtered)}")
print(f"Samples with NA distance: {len(df) - len(df_filtered)}")
print(f"\nDistance statistics:")
print(f"Mean: {np.mean(distances):.2f}")
print(f"Median: {np.median(distances):.2f}")
print(f"Std: {np.std(distances):.2f}")
print(f"Min: {np.min(distances):.2f}")
print(f"Max: {np.max(distances):.2f}")
print(f"\nDistance value counts:")
for dist, count in zip(unique_distances, counts):
  print(f"Distance {dist:.1f}: {count} samples ({count/len(df_filtered)*100:.1f}%)")