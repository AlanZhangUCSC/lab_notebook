import matplotlib.patches as mplpatches
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np
import argparse
import glob
import sys
import re

plt.style.use('BME163')

##########
## Data ##
##########

parser = argparse.ArgumentParser()
parser.add_argument('missing_data_summary', type=str, help='Input file containing missing data summary')
parser.add_argument('-t', '--title', type=str, required=True, help='Title of the plot')
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

missing_data_summary = args.missing_data_summary
title = args.title
output = args.output


connected_component_size_counts = []

read_cc_data = False
for line in open(missing_data_summary):
  if (not read_cc_data):
    if (line.startswith('Connected component size')):
      read_cc_data = True
    continue
  else:
    cc_size, count = map(int, line.strip().split())
    connected_component_size_counts.append((cc_size, count))

bin_size = 10
total_bins = 10
bins = [0 for _ in range(total_bins + 2)]
bin_labels = []


bins[0] = connected_component_size_counts[0][1]
bin_labels.append(f"1")
current_bin_index = 1
bin_labels.append(f"(1, {bin_size}]")
for cc_size, count in connected_component_size_counts[1:]:
  if (current_bin_index < total_bins + 1 and cc_size > bin_size * current_bin_index):
    current_bin_index += 1
    if current_bin_index != total_bins + 1:
      bin_labels.append(f"({bin_size * (current_bin_index - 1)}, {bin_size * current_bin_index}]")
    else:
      bin_labels.append(f">{bin_size * (current_bin_index - 1)}")
  bins[current_bin_index] += count

print(bins)
print(bin_labels)
print(sum(bins))

bins_in_thousands = [curbin / 1000 for curbin in bins]

print(bins_in_thousands)


############
## Figure ##
############
figure_size = (9, 8)
bottom_panel_size = (7, 4)
bottom_panel_coor = (1, 1.5)
top_panel_size = (7, 1)
top_panel_coor = (1, 1.5 + bottom_panel_size[1] + 0.5)

bottom_panel = plt.axes([bottom_panel_coor[0] / figure_size[0], bottom_panel_coor[1] / figure_size[1], bottom_panel_size[0] / figure_size[0], bottom_panel_size[1] / figure_size[1]])
top_panel = plt.axes([top_panel_coor[0] / figure_size[0], top_panel_coor[1] / figure_size[1], top_panel_size[0] / figure_size[0], top_panel_size[1] / figure_size[1]])


bottom_panel_y_limit = (max(bins_in_thousands[1:]) // 20 + 1) * 20
top_panel_y_limit = (max(bins_in_thousands) // 100 + 1) * 100
bottom_panel.set_ylim(0, bottom_panel_y_limit)
top_panel.set_ylim(bottom_panel_y_limit + 1, top_panel_y_limit)

bar_padding = 0.05
inter_bar_space = 0.02
x_range_per_bar = (1 - bar_padding * 2 - inter_bar_space * (len(bin_labels) - 1)) / len(bin_labels)
cur_x_start = bar_padding
x_ticks = []
for i, name in enumerate(bin_labels):
  x_ticks.append(cur_x_start + x_range_per_bar / 2)
  cur_x_start += x_range_per_bar + inter_bar_space

for x_start, count in zip(x_ticks, bins_in_thousands):
  bottom_panel.bar(x_start, count, width=x_range_per_bar, label=name, color='gray')
  top_panel.bar(x_start, count, width=x_range_per_bar, label=name, color='gray')

bottom_panel.set_xticks(x_ticks)
bottom_panel.set_xticklabels(bin_labels, rotation=45, ha='right')
bottom_panel.set_ylabel('Frequency (k)')
bottom_panel.set_xlabel('Neighborhood Size (Number of Connected Nodes)')
top_panel.set_yticks([top_panel_y_limit])
top_panel.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
top_panel.set_xticks([])
bottom_panel.set_yticks([i for i in range(0, int(bottom_panel_y_limit), 50)])

top_panel.spines['bottom'].set_visible(False)
bottom_panel.spines['top'].set_visible(False)

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
top_panel.plot([0, 1], [0, 0], transform=top_panel.transAxes, **kwargs)
bottom_panel.plot([0, 1], [1, 1], transform=bottom_panel.transAxes, **kwargs)

top_panel.set_title(title, fontsize=14, pad=15)

plt.savefig(output, dpi=600)
