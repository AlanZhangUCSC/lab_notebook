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


alignment_diffs = []
for parent in pairs_dict.keys():
  for child in pairs_dict[parent].keys():
    # print(parent, child, pairs_dict[parent][child])
    try:
      alignment_diffs.append(pairs_dict[parent][child][0] - pairs_dict[parent][child][1])
    except:
      continue


############
## Figure ##
############

plt.style.use('BME163')
figure_size = (6, 3)
panel_size = (figure_size[0] - 0.5, figure_size[1] - 0.5)
panel_coor = (0.4,0.25)
plt.figure(figsize=figure_size)
panel = plt.axes(
  [panel_coor[0] / figure_size[0],
   panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

panel.set_xlim(-100, len(alignment_diffs) + 100)
panel.set_ylim(-5000, 5000)
panel.set_xticks([])
panel.set_yticks([1000 * i for i in range(-5, 6)])
panel.set_yticklabels([f'{i}' for i in range(-5, 6)])
panel.set_ylabel('Alignment difference (kbp)')
panel.set_title(title)

sorted_alignment_diffs = sorted(alignment_diffs)
colors = ['skyblue' if diff == 0 else 'salmon' if diff > 0 else 'mediumseagreen' for diff in sorted_alignment_diffs]
panel.scatter(
  range(len(sorted_alignment_diffs)),
  sorted_alignment_diffs,
  c=colors,
  s=0.5
)

plt.savefig(output_file, dpi=600)
print(max(alignment_diffs))
print(min(alignment_diffs))
