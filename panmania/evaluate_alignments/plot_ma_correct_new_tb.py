import matplotlib.patches as mplpatches
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np
import argparse
import glob
import sys
import re

##########
## Data ##
##########

diff_file = sys.argv[1]
title = sys.argv[2]
output_file = sys.argv[3]

pairs_dict = defaultdict(lambda: defaultdict(tuple))
with open(diff_file, 'r') as f:
  next(f)
  for line in f:
    parent, child, panman_diff, minimap_diff, expected_diff_after_correction, mismatch_count_in_ma_blocks= line.strip().split('\t')
    try:
      pairs_dict[parent][child] = (int(panman_diff), int(minimap_diff), int(expected_diff_after_correction), int(mismatch_count_in_ma_blocks))
    except:
      continue

total_mismatches_all = 0
total_mismatches_after_correction = 0

internal_to_internal_diffs = []
internal_to_internal_diffs_after_correction = []
internal_to_leaf_diffs = []
internal_to_leaf_diffs_after_correction = []
diff_bins = [0 for i in range(11)]
panman_better_than_minimap = 0
diff_bins_after_correction = [0 for i in range(11)]
panman_better_than_minimap_after_correction = 0


for parent in pairs_dict.keys():
  for child in pairs_dict[parent].keys():
    panman_diff, minimap_diff, expected_diff_after_correction, mismatch_count_in_ma_blocks = pairs_dict[parent][child]
    mismatches_all = panman_diff - minimap_diff
    mismatches_after_correction = expected_diff_after_correction + mismatch_count_in_ma_blocks - minimap_diff
    total_mismatches_all += mismatches_all
    total_mismatches_after_correction += mismatches_after_correction
    if mismatches_all > 0:
      try: diff_bins[mismatches_all // 100] += 1
      except: diff_bins[-1] += 1
    else:
      panman_better_than_minimap += 1
    if mismatches_after_correction > 0:
      try: diff_bins_after_correction[mismatches_after_correction // 100] += 1
      except: diff_bins_after_correction[-1] += 1
    else:
      panman_better_than_minimap_after_correction += 1
    try:
      if child.startswith('node_'):
        internal_to_internal_diffs.append(mismatches_all)
        internal_to_internal_diffs_after_correction.append(mismatches_after_correction)
      else:
        internal_to_leaf_diffs.append(mismatches_all)
        internal_to_leaf_diffs_after_correction.append(mismatches_after_correction)
    except:
      continue

internal_to_internal_diffs = sorted(internal_to_internal_diffs)
internal_to_internal_diffs_after_correction = sorted(internal_to_internal_diffs_after_correction)
internal_to_leaf_diffs = sorted(internal_to_leaf_diffs)
internal_to_leaf_diffs_after_correction = sorted(internal_to_leaf_diffs_after_correction)

diff_bins = [panman_better_than_minimap] + diff_bins
diff_bins_after_correction = [panman_better_than_minimap_after_correction] + diff_bins_after_correction

############
## Figure ##
############

plt.style.use('BME163')
figure_size = (12, 7)
diff_distribution_panel_size = (6 - 0.5, 6 - 1 + 0.25)
diff_distribution_panel_coor = (0.4, 1 + 0.25)
diff_panel_size = (6 - 0.5, 6 / 2 - 0.5)
legend_panel_size = (1, 0.5)
legend_panel_coor = (diff_distribution_panel_coor[0] + diff_distribution_panel_size[0] - 3, diff_distribution_panel_coor[1] + diff_distribution_panel_size[1] - 0.5)
internal_to_internal_diff_panel_coor = (0.4 + diff_panel_size[0] + 0.5, 1 + 6 / 2)
internal_to_leaf_diff_panel_coor = (0.4 + diff_panel_size[0] + 0.5, 1 + 0.5 - 0.25)

plt.figure(figsize=figure_size)
diff_distribution_panel = plt.axes(
  [diff_distribution_panel_coor[0] / figure_size[0], diff_distribution_panel_coor[1] / figure_size[1],
   diff_distribution_panel_size[0] / figure_size[0], diff_distribution_panel_size[1] / figure_size[1]]
)

legend_panel = plt.axes(
  [legend_panel_coor[0] / figure_size[0], legend_panel_coor[1] / figure_size[1],
   legend_panel_size[0] / figure_size[0], legend_panel_size[1] / figure_size[1]]
)

internal_to_internal_diff_panel = plt.axes(
  [internal_to_internal_diff_panel_coor[0] / figure_size[0], internal_to_internal_diff_panel_coor[1] / figure_size[1],
   diff_panel_size[0] / figure_size[0], diff_panel_size[1] / figure_size[1]]
)

internal_to_leaf_diff_panel = plt.axes(
  [internal_to_leaf_diff_panel_coor[0] / figure_size[0], internal_to_leaf_diff_panel_coor[1] / figure_size[1],
   diff_panel_size[0] / figure_size[0], diff_panel_size[1] / figure_size[1]]
)

########################
## distribution panel ##
########################

diff_distribution_panel.set_xlim(0, 12-0.2)
diff_distribution_panel.set_ylim(0, (max(max(diff_bins), max(diff_bins_after_correction), panman_better_than_minimap, panman_better_than_minimap_after_correction) // 100 + 1) * 100)
for i in range(12):
  for offset, diff_bin, color in zip([0, 0.4], [diff_bins[i], diff_bins_after_correction[i]], ['salmon', 'teal']):
    diff_distribution_panel.add_patch(
      mplpatches.Rectangle((i + offset, 0), 0.4, diff_bin, color=color)
    )
    if diff_bin <= 100:
      diff_distribution_panel.text(i + offset + 0.2, diff_bin, f'{diff_bin}', ha='center', va='bottom')
x_ticks = []
x_ticklabels = []
for i in range(12):
  x_ticks.append(0.4 + i)
  if i == 0:
    x_ticklabels.append('<0')
  elif i < 10:
    x_ticklabels.append(f'{i*100}-{i*100+100}')
  else:
    x_ticklabels.append(f'>{i*100}')
diff_distribution_panel.set_xticks(x_ticks)
diff_distribution_panel.set_xticklabels(x_ticklabels, rotation=90)
if diff_distribution_panel.get_ylim()[1] > 5000:
  diff_distribution_panel.set_yticks([i * 500 for i in range(int(diff_distribution_panel.get_ylim()[1] // 500 + 1))])
else:
  diff_distribution_panel.set_yticks([i * 100 for i in range(int(diff_distribution_panel.get_ylim()[1] // 100 + 1))])
diff_distribution_panel.set_ylabel('Number of pairs')
diff_distribution_panel.set_xlabel('panMAN and minimap2 alignment difference')
diff_distribution_panel.set_title(f'Distribution of alignment differences\n({title})', fontsize=12)

##################
## legend panel ##
##################

legend_panel.set_xlim(0, 1)
legend_panel.set_ylim(0, 0.5)
legend_panel.set_xticks([])
legend_panel.set_yticks([])
legend_panel.set_xticklabels([])
legend_panel.set_yticklabels([])

legend_panel.add_patch(
  mplpatches.Rectangle((0.1, 0.3), 0.1, 0.1, color='salmon')
)
legend_panel.add_patch(
  mplpatches.Rectangle((0.1, 0.1), 0.1, 0.1, color='teal')
)
legend_panel.text(0.3, 0.3, '(panMAN - minimap2) alignment before correction', ha='left', va='baseline')
legend_panel.text(0.3, 0.1, '(panMAN - minimap2) alignment after correction', ha='left', va='baseline')

legend_panel.axis('off')

#################
## diff panels ##
#################

for panel, diff_sets, xlabel in zip([internal_to_internal_diff_panel, internal_to_leaf_diff_panel], [[internal_to_internal_diffs, internal_to_internal_diffs_after_correction], [internal_to_leaf_diffs, internal_to_leaf_diffs_after_correction]], ['internal-to-internal', 'internal-to-leaf']):
  diffs, diffs_after_correction = diff_sets
  assert(len(diffs) == len(diffs_after_correction))

  panel.set_xlim(-100, len(diffs) + 100)
  panel.set_ylim(-5, 5)
  panel.set_xticks([])
  panel.set_yticks([i * 1000 for i in range(-5, 6)])
  panel.set_yticklabels([f'{i}' for i in range(-5, 6)])
  panel.set_yticks([i * 1000 for i in np.arange(-5, 6, 0.5)], minor=True)
  panel.set_ylabel('Alignment difference (kbp)')
  panel.set_xlabel(f'{xlabel} ({len(diffs)} pairs)', fontsize=8)

  panel.grid(True, which='both')
  panel.grid(which='minor', alpha=0.2, linestyle='--')
  panel.grid(which='major', alpha=0.5, linestyle='--')
  panel.plot([-100, len(diffs) + 100], [0, 0], color='black', linewidth=0.5, linestyle='--')

internal_to_internal_diff_panel.set_title(f'panMAN alignment and minimap2 alignment for diffs > 0\n({title})', fontsize=12)

internal_to_internal_diff_panel.scatter(
  range(len(internal_to_internal_diffs)),
  internal_to_internal_diffs,
  color='salmon',
  s=0.5
)

internal_to_internal_diff_panel.scatter(
  range(len(internal_to_internal_diffs)),
  internal_to_internal_diffs_after_correction,
  color='teal',
  s=0.5
)

internal_to_leaf_diff_panel.scatter(
  range(len(internal_to_leaf_diffs)),
  internal_to_leaf_diffs,
  color='salmon',
  s=0.5
)

internal_to_leaf_diff_panel.scatter(
  range(len(internal_to_leaf_diffs)),
  internal_to_leaf_diffs_after_correction,
  color='teal',
  s=0.5
)


plt.savefig(output_file, dpi=600)
print(diff_bins)
print(diff_bins_after_correction)