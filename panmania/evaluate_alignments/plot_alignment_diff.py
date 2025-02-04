import matplotlib.patches as mplpatches
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import argparse
import glob
import sys
import re

##########
## Data ##
##########

diff_file = sys.argv[1]
pairs_file = sys.argv[2]
title = sys.argv[3]
output_file = sys.argv[4]

pairs_dict = defaultdict(lambda: defaultdict(tuple))
with open(pairs_file, 'r') as f:
  for line in f:
    parent, child = line.strip().split(' ')
    pairs_dict[parent][child] = None

with open(diff_file, 'r') as f:
  next(f)
  for line in f:
    parent, child, panman_diff, mafft_diff = line.strip().split('\t')
    pairs_dict[parent][child] = (int(panman_diff), int(mafft_diff))


internal_to_internal_diffs = []
internal_to_internal_scores = []
internal_to_leaf_diffs = []
internal_to_leaf_scores = []
for parent in pairs_dict.keys():
  for child in pairs_dict[parent].keys():
    # print(parent, child, pairs_dict[parent][child])
    try:
      if child.startswith('node_'):
        internal_to_internal_diffs.append(pairs_dict[parent][child][0] - pairs_dict[parent][child][1])
        internal_to_internal_scores.append((pairs_dict[parent][child][0], pairs_dict[parent][child][1]))
      else:
        internal_to_leaf_diffs.append(pairs_dict[parent][child][0] - pairs_dict[parent][child][1])
        internal_to_leaf_scores.append((pairs_dict[parent][child][0], pairs_dict[parent][child][1]))
    except:
      continue

internal_to_internal_diffs, internal_to_internal_scores = zip(*sorted(zip(internal_to_internal_diffs, internal_to_internal_scores)))
internal_to_internal_first_zero_index = next((i for i, diff in enumerate(internal_to_internal_diffs) if diff == 0), None)
internal_to_internal_first_positive_index = next((i for i, diff in enumerate(internal_to_internal_diffs) if diff > 0), None)

internal_to_leaf_diffs, internal_to_leaf_scores = zip(*sorted(zip(internal_to_leaf_diffs, internal_to_leaf_scores)))
internal_to_leaf_first_zero_index = next((i for i, diff in enumerate(internal_to_leaf_diffs) if diff == 0), None)
internal_to_leaf_first_positive_index = next((i for i, diff in enumerate(internal_to_leaf_diffs) if diff > 0), None)



############
## Figure ##
############

plt.style.use('BME163')
figure_size = (12, 6)
panel_size = (6 - 0.5, figure_size[1] / 2 - 0.5)

internal_to_internal_panel_coor = (0.4, figure_size[1] / 2)
internal_to_leaf_panel_coor = (0.4, 0.5 - 0.25)

plt.figure(figsize=figure_size)
internal_to_internal_panel = plt.axes(
  [internal_to_internal_panel_coor[0] / figure_size[0],
   internal_to_internal_panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

internal_to_leaf_panel = plt.axes(
  [internal_to_leaf_panel_coor[0] / figure_size[0],
   internal_to_leaf_panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

for panel, diffs, cutoffs, xlabel in zip([internal_to_internal_panel, internal_to_leaf_panel], [internal_to_internal_diffs, internal_to_leaf_diffs], [(internal_to_internal_first_zero_index, internal_to_internal_first_positive_index), (internal_to_leaf_first_zero_index, internal_to_leaf_first_positive_index)], ['internal-to-internal', 'internal-to-leaf']):
  ylim_kb = 5
  panel.set_xlim(-100, len(diffs) + 100)
  panel.set_ylim(-ylim_kb * 1000, ylim_kb * 1000)
  panel.set_xticks([])
  panel.set_yticks([1000 * i for i in range(-ylim_kb, ylim_kb + 1)])
  panel.set_yticklabels([f'{i}' for i in range(-ylim_kb, ylim_kb + 1)])
  panel.set_yticks([1000 * i for i in np.arange(-ylim_kb, ylim_kb, 0.5)], minor=True)
  panel.set_ylabel('Alignment difference (kbp)')
  panel.set_xlabel(f'{xlabel} ({len(diffs)} pairs)', fontsize=6)

  panel.grid(True, which='both')
  panel.grid(which='minor', alpha=0.2, linestyle='--')
  panel.grid(which='major', alpha=0.5, linestyle='--')
  first_zero_index, first_positive_index = cutoffs
  panel.plot([first_zero_index, first_zero_index], [-ylim_kb * 1000, ylim_kb * 1000], color='black', linewidth=0.5, linestyle='--')
  panel.plot([first_positive_index, first_positive_index], [-ylim_kb * 1000, ylim_kb * 1000], color='black', linewidth=0.5, linestyle='--')


internal_to_internal_colors = ['skyblue' if diff == 0 else 'yellowgreen' if diff < 0 else 'salmon' for diff in internal_to_internal_diffs]
internal_to_leaf_colors = ['skyblue' if diff == 0 else 'yellowgreen' if diff < 0 else 'salmon' for diff in internal_to_leaf_diffs]

internal_to_internal_panel.set_title(f'{title} (panMAN alignment - MAFFT alignment)', fontsize=14)

internal_to_internal_panel.scatter(
  range(len(internal_to_internal_diffs)),
  internal_to_internal_diffs,
  c=internal_to_internal_colors,
  s=0.5
)


internal_to_leaf_panel.scatter(
  range(len(internal_to_leaf_diffs)),
  internal_to_leaf_diffs,
  c=internal_to_leaf_colors,
  s=0.5
)

internal_to_internal_diff_diffs = internal_to_internal_diffs[internal_to_internal_first_positive_index:]
internal_to_internal_diff_scores = internal_to_internal_scores[internal_to_internal_first_positive_index:]
internal_to_leaf_diff_diffs = internal_to_leaf_diffs[internal_to_leaf_first_positive_index:]
internal_to_leaf_diff_scores = internal_to_leaf_scores[internal_to_leaf_first_positive_index:]

internal_to_internal_diff_panel_coor = (0.4 + panel_size[0] + 0.5, figure_size[1] / 2)
internal_to_leaf_diff_panel_coor = (0.4 + panel_size[0] + 0.5, 0.5 - 0.25)
internal_to_internal_diff_panel = plt.axes(
  [internal_to_internal_diff_panel_coor[0] / figure_size[0],
   internal_to_internal_diff_panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

internal_to_leaf_diff_panel = plt.axes(
  [internal_to_leaf_diff_panel_coor[0] / figure_size[0],
   internal_to_leaf_diff_panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

for panel, diffs, xlabel in zip([internal_to_internal_diff_panel, internal_to_leaf_diff_panel], [internal_to_internal_diff_diffs, internal_to_leaf_diff_diffs], ['internal-to-internal', 'internal-to-leaf']):
  ylim_kb = 5
  panel.set_xlim(-100, len(diffs) + 100)
  panel.set_ylim(-ylim_kb * 1000, ylim_kb * 1000)
  panel.set_xticks([])
  panel.set_yticks([1000 * i for i in range(-ylim_kb, ylim_kb + 1)])
  panel.set_yticklabels([f'{i}' for i in range(-ylim_kb, ylim_kb + 1)])
  panel.set_yticks([1000 * i for i in np.arange(-ylim_kb, ylim_kb, 0.5)], minor=True)
  panel.set_ylabel('Alignment difference (kbp)')
  panel.set_xlabel(f'{xlabel} ({len(diffs)} pairs)', fontsize=6)

  panel.grid(True, which='both')
  panel.grid(which='minor', alpha=0.2, linestyle='--')
  panel.grid(which='major', alpha=0.5, linestyle='--')
  panel.plot([-100, len(diffs) + 100], [0, 0], color='black', linewidth=0.5, linestyle='--')

internal_to_internal_diff_panel.set_title(f'{title} (panMAN alignment and MAFFT alignment for diffs > 0)', fontsize=14)

internal_to_internal_diff_panel.scatter(
  range(len(internal_to_internal_diff_diffs)),
  internal_to_internal_diff_diffs,
  color='salmon',
  s=0.5
)

internal_to_internal_diff_panel.scatter(
  range(len(internal_to_internal_diff_scores)),
  [diff[0] for diff in internal_to_internal_diff_scores],
  color='slateblue',
  s=0.5
)

internal_to_internal_diff_panel.scatter(
  range(len(internal_to_internal_diff_scores)),
  [diff[1] for diff in internal_to_internal_diff_scores],
  color='darkseagreen',
  s=0.5
)

internal_to_leaf_diff_panel.scatter(
  range(len(internal_to_leaf_diff_diffs)),
  internal_to_leaf_diff_diffs,
  color='salmon',
  s=0.5
)

internal_to_leaf_diff_panel.scatter(
  range(len(internal_to_leaf_diff_scores)),
  [diff[0] for diff in internal_to_leaf_diff_scores],
  color='slateblue',
  s=0.5
)

internal_to_leaf_diff_panel.scatter(
  range(len(internal_to_leaf_diff_scores)),
  [diff[1] for diff in internal_to_leaf_diff_scores],
  color='darkseagreen',
  s=0.5
)

plt.savefig(output_file, dpi=600)
