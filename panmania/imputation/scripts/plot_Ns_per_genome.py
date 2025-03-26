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
parser.add_argument('-m', '--missing-data-summary', type=str, nargs='+', required=True, help='Input file containing missing data summary')
parser.add_argument('-s', '--summary', type=str, nargs='+', required=True, help='Input file containing normal data summary')
parser.add_argument('-n', '--names', type=str, nargs='+', required=True, help='Names of the species')
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

missing_data_summary = args.missing_data_summary
summary = args.summary
names = args.names
output = args.output

assert len(missing_data_summary) == len(summary) == len(names)

average_length_dict = {
  'RSV_4000': 15200,
  'SARS_8M': 30000,
  'SARS_20000': 30000,
  'HIV_20000': 9750,
  'ECOLI_1000': 5000000,
  'KLEBS_1000': 5491870,
  'TB_400': 4400000,
}



data = {}
for name, missing_data_summary_file, summary_file in zip(names, missing_data_summary, summary):
  data[name] = {}
  num_nodes = 0
  num_Ns = 0

  normal_data_summary_fh = open(summary_file, 'r')
  num_nodes = int(normal_data_summary_fh.readline().split(' ')[-1])
  normal_data_summary_fh.close()

  missing_data_summary_lines = open(missing_data_summary_file, 'r').readlines()
  num_Ns = int(missing_data_summary_lines[2].strip().split(' ')[-1])

  data[name]['N_per_base'] = num_Ns / (average_length_dict[name] * num_nodes)
  data[name]['N_per_genome'] = num_Ns / num_nodes

sorted_names = sorted(data.keys(), key=lambda x: data[x]['N_per_genome'], reverse=True)
data = {name: data[name] for name in sorted_names}
for name in sorted_names:
  print(f"{name}: {data[name]['N_per_base']} (average {data[name]['N_per_genome']:.2f} Ns per genomes)")


############
## Figure ##
############
figure_size = (8, 7)
bottom_panel_size = (7, 4)
bottom_panel_coor = (0.5, 0.75)
top_panel_size = (7, 1)
top_panel_coor = (0.5, 0.5 + bottom_panel_size[1] + 0.5)

bottom_panel = plt.axes([bottom_panel_coor[0] / figure_size[0], bottom_panel_coor[1] / figure_size[1], bottom_panel_size[0] / figure_size[0], bottom_panel_size[1] / figure_size[1]])
top_panel = plt.axes([top_panel_coor[0] / figure_size[0], top_panel_coor[1] / figure_size[1], top_panel_size[0] / figure_size[0], top_panel_size[1] / figure_size[1]])

bottom_panel.set_ylim(0, 100)
top_panel.set_ylim(100, 3000)

bar_padding = 0.05
inter_bar_space = 0.02
x_range_per_bar = (1 - bar_padding * 2 - inter_bar_space * (len(data.keys()) - 1)) / len(data.keys())
cur_x_start = bar_padding
x_ticks = []
for i, name in enumerate(data.keys()):
  data[name]['x_start'] = cur_x_start + x_range_per_bar / 2
  x_ticks.append(cur_x_start + x_range_per_bar / 2)
  cur_x_start += x_range_per_bar + inter_bar_space

for name in data.keys():
  bottom_panel.bar(data[name]['x_start'], data[name]['N_per_genome'], width=x_range_per_bar, label=name, color='gray')
  top_panel.bar(data[name]['x_start'], data[name]['N_per_genome'], width=x_range_per_bar, label=name, color='gray')
  if (data[name]['N_per_genome'] < 1):
    bottom_panel.text(data[name]['x_start'], data[name]['N_per_genome'], f"{data[name]['N_per_genome']:.2f}", ha='center', va='bottom')

print(x_ticks)
bottom_panel.set_xticks(x_ticks)
bottom_panel.set_xticklabels(data.keys(), rotation=15, ha='center')
top_panel.set_yticks([3000])
top_panel.set_xticks([])
bottom_panel.set_yticks([i for i in range(0, 100, 20)])

top_panel.spines['bottom'].set_visible(False)
bottom_panel.spines['top'].set_visible(False)

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
top_panel.plot([0, 1], [0, 0], transform=top_panel.transAxes, **kwargs)
bottom_panel.plot([0, 1], [1, 1], transform=bottom_panel.transAxes, **kwargs)

top_panel.set_title('Average number of Ns per genome', fontsize=14, pad=15)

plt.savefig(output, dpi=600)
