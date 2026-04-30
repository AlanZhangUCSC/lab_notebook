import matplotlib.patches as mplpatches
from collections import defaultdict
import matplotlib.pyplot as plt
import dendropy
import sys
import argparse
import numpy as np

def get_depth(node, current_depth=0):
  if not node.child_nodes():
    return current_depth
  return max(get_depth(child, current_depth + 1) for child in node.child_nodes())

def get_node_level(node):
  level = 0
  while node.parent_node is not None:
    node = node.parent_node
    level += 1
  return level

def get_distance_to_deepest_leaf(node):
  if not node.child_nodes():
    return 0
  return max(1 + get_distance_to_deepest_leaf(child) for child in node.child_nodes())

def draw_branch(child_node, parent_node, node_coord, color, line_width=0.004):
  horizontal_line_len = abs(node_coord[child_node]['x_coord'] - node_coord[parent_node]['x_coord'])
  vertical_line_len = abs(node_coord[child_node]['y_coord'] - node_coord[parent_node]['y_coord'])

  horizontal_line = mplpatches.Rectangle(
    (node_coord[parent_node]['x_coord'] - line_width / 2,
     node_coord[child_node]['y_coord'] - line_width / 2),
    horizontal_line_len,
    line_width,
    facecolor=color,
    linewidth=0)

  if node_coord[child_node]['y_coord'] < node_coord[parent_node]['y_coord']:
    vertical_line = mplpatches.Rectangle(
      (node_coord[parent_node]['x_coord'] - line_width / 2,
       node_coord[child_node]['y_coord'] - line_width / 2),
      line_width,
      vertical_line_len,
      facecolor=color,
      linewidth=0)
  else:
    vertical_line = mplpatches.Rectangle(
      (node_coord[parent_node]['x_coord'] - line_width / 2,
       node_coord[parent_node]['y_coord'] - line_width / 2),
      line_width,
      vertical_line_len,
      facecolor=color,
      linewidth=0)

  return horizontal_line, vertical_line

def hex_to_rgb(hex_color):
  hex_color = hex_color.lstrip('#')
  return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def interpolate_color(t, color_stops, stop_positions):
  for i in range(len(stop_positions) - 1):
    if t <= stop_positions[i + 1]:
      local_t = (t - stop_positions[i]) / (stop_positions[i + 1] - stop_positions[i])
      c0 = hex_to_rgb(color_stops[i])
      c1 = hex_to_rgb(color_stops[i + 1])
      return tuple(c0[j] + local_t * (c1[j] - c0[j]) for j in range(3))
  return hex_to_rgb(color_stops[-1])

def get_node_color(read_count, max_read_count, min_read_count=0):
  color_stops = ["#D3D3D3", "#DE9E9E", "#E96969", "#F43434", "#FF0000"]
  stop_positions = [0.0, 0.25, 0.5, 0.75, 1.0]
  t = (read_count - min_read_count) / (max_read_count - min_read_count) if (max_read_count - min_read_count) > 0 else 0
  t = max(0.0, min(1.0, t))
  return interpolate_color(t, color_stops, stop_positions)

def ladderize(tree):
  num_tips = {}
  for node in tree.seed_node.postorder_iter():
    if node.is_leaf():
      num_tips[node] = 1
    else:
      num_tips[node] = sum(num_tips[child] for child in node.child_nodes())

  stack = [tree.seed_node]
  while stack:
    node = stack.pop()
    if not node.is_leaf():
      children = node.child_nodes()
      children.sort(key=lambda n: num_tips[n], reverse=True)
      node.set_child_nodes(children)
      stack.extend(children)

  return [node.taxon.label if node.taxon else node.label for node in reversed(list(tree.seed_node.preorder_iter())) if node.is_leaf()]


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--nwk', type=str, required=True, help='Newick tree file')
parser.add_argument('-n', '--node_read_counts_file', type=str, required=True, help='Node read counts file')
parser.add_argument('-l', '--node_read_counts_lca_file', type=str, required=True, help='Node read counts LCA file')
parser.add_argument('-m', '--metadata_file', type=str, required=False, default='', help='Metadata file')
parser.add_argument('--color-node-labels', type=str, required=False, default='', help='Color node labels by metadata column')
parser.add_argument('-o', '--out_prefix', type=str, required=True, help='Output prefix')
parser.add_argument('-c', '--color-by', type=str, default='lca_count', choices=['read_count', 'lca_count', 'lca_subtree_count', 'None'], help='Color by read count or LCA count')
parser.add_argument('-s', '--size-by', type=str, default='lca_subtree_count', choices=['read_count', 'lca_count', 'lca_subtree_count', 'None'], help='Size by read count or LCA count')
parser.add_argument('-a', '--above-branch-label', type=str, default='lca_count', choices=['read_count', 'lca_count', 'lca_subtree_count', None], help='Above branch label')
parser.add_argument('-b', '--below-branch-label', type=str, default='lca_subtree_count', choices=['read_count', 'lca_count', 'lca_subtree_count', None], help='Below branch label')
args = parser.parse_args()

nwk_string = ''
with open(args.nwk, 'r') as f:
  nwk_string = f.read().strip()

tree = dendropy.Tree.get(data=nwk_string, schema="newick")
leaf_nodes = [node for node in tree.leaf_node_iter()]
num_leaves = len(leaf_nodes)
max_depth = get_depth(tree.seed_node)
leaves_order = ladderize(tree)

print(f"Number of leaf nodes: {num_leaves}")
print(f"Max tree depth (edges): {max_depth}")

distance_to_root = {}
distance_to_deepest_leaf = {}
for node in tree.postorder_node_iter():
  label = node.taxon.label if node.taxon else node.label
  distance_to_root[label] = get_node_level(node)
  distance_to_deepest_leaf[label] = get_distance_to_deepest_leaf(node)

if args.metadata_file != '':
  metadata = {}
  color_by_column = -1
  with open(args.metadata_file, 'r') as f:
    for i, line in enumerate(f):
      if i == 0:
        for j, column in enumerate(line.split('\t')):
          column = column.strip()
          if column == '': continue
          if column == args.color_node_labels:
            color_by_column = j
            break
      else:
        if color_by_column == -1:
          print(f"Error: {args.color_node_labels} not found in metadata file")
          exit(1)
        columns = line.split('\t')
        metadata[columns[0].replace('_', ' ')] = columns[color_by_column]

if args.metadata_file != '' and args.color_node_labels != '':
  unique_categories = sorted(set(metadata.values()))
  n_cats = len(unique_categories)
  if n_cats <= 20:
    cmap = plt.cm.get_cmap('tab20', 20)
    category_colors = {cat: cmap(i) for i, cat in enumerate(unique_categories)}
  else:
    category_colors = {cat: plt.cm.hsv(i / n_cats) for i, cat in enumerate(unique_categories)}

node_read_counts = {}
max_read_count = 0
min_read_count = float('inf')
with open(args.node_read_counts_file, 'r') as f:
  for line in f:
    node_label, family, count = line.strip().split('\t')
    for node_label_clean in node_label.split(','):
      node_label_clean = node_label_clean.replace('_', ' ')
      node_read_counts[node_label_clean] = int(count)
      max_read_count = max(max_read_count, int(count))
      min_read_count = min(min_read_count, int(count))

node_read_counts_lca = {}
with open(args.node_read_counts_lca_file, 'r') as f:
  for line in f:
    node_label, count_lca = line.strip().split()
    for node_label_clean in node_label.split(','):
      node_label_clean = node_label_clean.replace('_', ' ')
      count_lca = int(count_lca)
      node_read_counts_lca[node_label_clean] = count_lca

tree_labels = set()
for node in tree.postorder_node_iter():
  tree_labels.add(node.taxon.label if node.taxon else node.label)

lca_labels = set(node_read_counts_lca.keys())
print("In tree but not in LCA file:", tree_labels - lca_labels)
print("In LCA file but not in tree:", lca_labels - tree_labels)

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

# --- Adaptive sizing parameters ---
# Scale text and markers based on how many leaves we have.
# The idea: inter_node_y_coord is the vertical gap between adjacent leaves
# in normalized [0,1] coords. We use fractions of that gap to place text
# so nothing overlaps regardless of tree size.

node_buffer = 0.02
inter_node_x_coord = (1 - node_buffer * 2) / max_depth
inter_node_y_coord = (1 - node_buffer * 2) / max(1, num_leaves - 1)

# Text size: scale with vertical gap, clamped to a reasonable range.
# At ~20 leaves on a 7.5in figure, 4pt works well. Scale linearly from there.
text_size = max(2.0, min(6.0, inter_node_y_coord * 80))
count_text_size = text_size * 0.7

# Node marker sizes: wider range for better visual distinction, sqrt-scaled
# so area is proportional to count.
min_node_size = max(1.5, inter_node_y_coord * 30)
max_node_size = min_node_size * 2.5

branch_width = max(0.001, inter_node_y_coord * 0.08)

leaf_label_pad = inter_node_x_coord * 0.08

# --- Compute node coordinates ---
node_coord = {}
for i, node_label in enumerate(leaves_order):
  node_coord[node_label] = {
    'x_coord': (1 - node_buffer) - (distance_to_deepest_leaf[node_label] * (1 - node_buffer * 2) / max_depth),
    'y_coord': (1 - node_buffer) - (i * inter_node_y_coord),
  }

for node in tree.postorder_node_iter():
  node_label = node.taxon.label if node.taxon else node.label
  if node_label in node_coord:
    continue
  children_y = [node_coord[child.taxon.label if child.taxon else child.label]['y_coord']
                 for child in node.child_nodes()]
  node_coord[node_label] = {
    'x_coord': (1 - node_buffer) - (distance_to_deepest_leaf[node_label] * (1 - node_buffer * 2) / max_depth),
    'y_coord': sum(children_y) / len(children_y),
  }

# --- Plot ---
plt.style.use('BME163')
plt.rcParams['svg.fonttype'] = 'none'

figure_size = (4.7, 7.5)
plt.figure(figsize=figure_size)

panel_buffer = 0.1
panel_size = (figure_size[0] - 2 * panel_buffer, figure_size[1] - 2 * panel_buffer)
panel_coor = (panel_buffer, panel_buffer)
panel = plt.axes([
  panel_coor[0] / figure_size[0],
  panel_coor[1] / figure_size[1],
  panel_size[0] / figure_size[0],
  panel_size[1] / figure_size[1]])

# Draw branches
for node in tree.postorder_node_iter():
  node_label = node.taxon.label if node.taxon else node.label
  if args.color_by == 'read_count':
    branch_color = get_node_color(node_read_counts.get(node_label, 0), max_read_count, min_read_count)
  elif args.color_by == 'lca_count':
    branch_color = get_node_color(node_lca_count.get(node_label, 0), max_node_lca_count, min_node_lca_count)
  elif args.color_by == 'lca_subtree_count':
    branch_color = get_node_color(node_subtree_lca_count.get(node_label, 0), max_node_subtree_lca_count, min_node_subtree_lca_count)
  else:
    branch_color = '#D3D3D3'

  if node.parent_node is None:
    root_branch_len = inter_node_x_coord * 0.5
    panel.add_patch(mplpatches.Rectangle(
      (node_coord[node_label]['x_coord'] - root_branch_len,
       node_coord[node_label]['y_coord'] - branch_width / 2),
      root_branch_len,
      branch_width,
      facecolor=branch_color,
      linewidth=0))
    continue

  parent_label = node.parent_node.taxon.label if node.parent_node.taxon else node.parent_node.label
  horizontal_line, vertical_line = draw_branch(node_label, parent_label, node_coord, branch_color, branch_width)
  panel.add_patch(horizontal_line)
  panel.add_patch(vertical_line)

# Draw text labels and node markers.
# To place text snugly against the node circle, we convert the marker radius
# from points (display units) into data-coordinate units. This way the gap
# between marker edge and text is consistent regardless of figure/panel size.

fig = plt.gcf()
renderer = fig.canvas.get_renderer()
ax_bbox = panel.get_window_extent(renderer)
ax_width_px = ax_bbox.width
ax_height_px = ax_bbox.height
xlim = panel.get_xlim()
ylim = panel.get_ylim()
pts_per_data_x = ax_width_px / (xlim[1] - xlim[0])
pts_per_data_y = ax_height_px / (ylim[1] - ylim[0])

text_pad_data_x = 2.0 / pts_per_data_x
text_pad_data_y = 0.2 / pts_per_data_y

leaf_set = set(leaves_order)

for node_label in node_coord:
  read_count = node_read_counts.get(node_label, 0)
  lca_count = node_lca_count.get(node_label, 0)
  lca_subtree_count = node_subtree_lca_count.get(node_label, 0)


  if args.color_by == 'read_count':
    color = get_node_color(read_count, max_read_count, min_read_count)
  elif args.color_by == 'lca_count':
    color = get_node_color(lca_count, max_node_lca_count, min_node_lca_count)
  elif args.color_by == 'lca_subtree_count':
    color = get_node_color(lca_subtree_count, max_node_subtree_lca_count, min_node_subtree_lca_count)
  else:
    color = '#D3D3D3'

  if args.size_by == 'lca_count':
    if max_node_lca_count > 0:
      frac = lca_count / max_node_lca_count
    else:
      frac = 0
    marker_size = min_node_size + np.sqrt(frac) * (max_node_size - min_node_size)
  elif args.size_by == 'read_count':
    if max_read_count > 0:
      frac = read_count / max_read_count
    else:
      frac = 0
    marker_size = min_node_size + np.sqrt(frac) * (max_node_size - min_node_size)
  elif args.size_by == 'lca_subtree_count':
    if max_node_subtree_lca_count > 0:
      frac = lca_subtree_count / max_node_subtree_lca_count
    else:
      frac = 0
    marker_size = min_node_size + np.sqrt(frac) * (max_node_size - min_node_size)
  else:
    marker_size = min_node_size

  nx = node_coord[node_label]['x_coord']
  ny = node_coord[node_label]['y_coord']

  panel.plot(nx, ny,
    marker='o',
    markersize=marker_size,
    color=color,
    markeredgewidth=0.1,
    markeredgecolor="#716f6e",
    linewidth=0)

  marker_radius_pts = marker_size / 2.0
  marker_radius_x = marker_radius_pts / pts_per_data_x
  marker_radius_y = marker_radius_pts / pts_per_data_y

  above_branch_label = ''
  if args.above_branch_label == 'read_count':
    above_branch_label = f'{read_count}'
  elif args.above_branch_label == 'lca_count':
    above_branch_label = f'{lca_count}'
  elif args.above_branch_label == 'lca_subtree_count':
    above_branch_label = f'{lca_subtree_count}'

  below_branch_label = ''
  if args.below_branch_label == 'read_count':
    below_branch_label = f'({read_count})'
  elif args.below_branch_label == 'lca_count':
    below_branch_label = f'({lca_count})'
  elif args.below_branch_label == 'lca_subtree_count':
    below_branch_label = f'({lca_subtree_count})'


  panel.text(
    nx - marker_radius_x - text_pad_data_x,
    ny + text_pad_data_y,
    f'{above_branch_label}',
    fontfamily='sans',
    fontsize=count_text_size,
    horizontalalignment='right',
    verticalalignment='bottom')

  panel.text(
    nx - marker_radius_x - text_pad_data_x,
    ny - text_pad_data_y,
    f'{below_branch_label}',
    fontfamily='sans',
    style='italic',
    fontsize=count_text_size,
    horizontalalignment='right',
    verticalalignment='top')


  if node_label in leaf_set:
    label_color = 'black'
    if args.metadata_file != '' and args.color_node_labels != '':
      category = metadata.get(node_label, None)
      if category is not None and category in category_colors:
        label_color = category_colors[category]

    panel.text(
      nx + marker_radius_x + leaf_label_pad,
      ny,
      node_label,
      fontfamily='sans',
      fontsize=text_size,
      fontweight='bold',
      verticalalignment='center',
      horizontalalignment='left',
      color=label_color)

panel.set_axis_off()


if args.metadata_file != '' and args.color_node_labels != '':
  from matplotlib.lines import Line2D
  legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=category_colors[cat],
           markersize=max(3, text_size * 0.8), label=cat, linewidth=0)
    for cat in unique_categories
  ]
  panel.legend(
    handles=legend_handles,
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),
    fontsize=max(2.0, text_size * 0.8),
    frameon=False,
    handletextpad=0.3,
    labelspacing=0.4,
    borderpad=0.2)
    
plt.savefig(f'{args.out_prefix}.png', dpi=1000, bbox_inches='tight')
plt.savefig(f'{args.out_prefix}.svg', dpi=1000, bbox_inches='tight')