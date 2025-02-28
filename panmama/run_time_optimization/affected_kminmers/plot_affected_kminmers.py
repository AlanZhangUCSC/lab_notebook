import matplotlib.patches as mplpatches
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np
import argparse
import json
import glob
import sys
import re


########
## IO ##
########

input_file = sys.argv[1]
output_file = sys.argv[2]

##########
## Data ##
##########

class read_change:
  def __init__(self, node_id, read_index):
    self.node_id = node_id
    self.read_index = read_index
    self.changes = []
    self.kminmer_hits_before = []
    self.kminmer_hits_after = []
    self.kminmers_changes_added_to_chain = []
    self.kminmers_changes_removed_from_chain = []
    self.kminmers_changes_unchanged = []
  
  def add_kminmer_changes(self, kminmer_changes):
    self.changes.extend(kminmer_changes)

  def add_kminmer_hits(self, kminmer_hits, updated):
    if updated:
      self.kminmer_hits_after = kminmer_hits
    else:
      self.kminmer_hits_before = kminmer_hits

  def update_kminmers_change_types(self):
    if len(self.kminmer_hits_before) > len(self.kminmer_hits_after):
      self.kminmer_hits_after.extend([False] * (len(self.kminmer_hits_before) - len(self.kminmer_hits_after)))
    elif len(self.kminmer_hits_before) < len(self.kminmer_hits_after):
      self.kminmer_hits_before.extend([False] * (len(self.kminmer_hits_after) - len(self.kminmer_hits_before)))

    if len(self.kminmer_hits_before) == 0 and len(self.kminmer_hits_after) == 0:
      return

    for index in self.changes:
      try:
        if self.kminmer_hits_before[index] and not self.kminmer_hits_after[index]:
          self.kminmers_changes_removed_from_chain.append(index)
        elif not self.kminmer_hits_before[index] and self.kminmer_hits_after[index]:
          self.kminmers_changes_added_to_chain.append(index)
        else:
          self.kminmers_changes_unchanged.append(index)
      except:
        self.kminmers_changes_unchanged.append(index)

    self.kminmers_changes_added_to_chain.sort()
    self.kminmers_changes_removed_from_chain.sort()
    self.kminmers_changes_unchanged.sort()

read_changes = defaultdict(lambda: defaultdict(read_change))
cur_updated = None
cur_node_id = None
cur_read_index = None

processed_nodes = 0
with open(input_file, 'r') as f: lines = f.readlines()
for line in lines:
  line = line.strip()
  if line.startswith('>'):
    fields = line.split(' ')[1].split('_')
    cur_node_id = '_'.join(fields[:-1])
    if fields[-1] == 'before':
      processed_nodes += 1
      print(f'Processing {processed_nodes}th node, {cur_node_id}')
      cur_updated = False
    else:
      cur_updated = True
    continue

  if line.startswith('@'):
    cur_read_index = int(line[1:].split(',')[0])
    if cur_updated:
      assert cur_node_id in read_changes.keys() and cur_read_index in read_changes[cur_node_id].keys()
      continue
    read_changes[cur_node_id][cur_read_index] = read_change(cur_node_id, cur_read_index)
  elif line.startswith('['):
    minichains_str = line.split(' ')
    minichains = []
    for minichain_str in minichains_str:
      minichain = [None, None]
      minichain[0] = int(minichain_str.split(',')[0][1:])
      minichain[1] = int(minichain_str.split(',')[1][:-1])
      minichains.append(minichain)
    kminmer_hits = [False] * (minichains[-1][1] + 1)
    for minichain in minichains:
      for i in range(minichain[0], minichain[1]+1):
        kminmer_hits[i] = True
    read_changes[cur_node_id][cur_read_index].add_kminmer_hits(kminmer_hits, cur_updated)
  elif line != '':
    if cur_updated: continue
    read_changes[cur_node_id][cur_read_index].add_kminmer_changes(map(int, line.split(' ')))



  

change_types = defaultdict(lambda: defaultdict(int))
for node_id in read_changes.keys():
  for read_index in read_changes[node_id].keys():
    cur_read_change = read_changes[node_id][read_index]
    cur_read_change.update_kminmers_change_types()
    num_changes = len(cur_read_change.changes)
    num_changes_added = len(cur_read_change.kminmers_changes_added_to_chain)
    num_changes_removed = len(cur_read_change.kminmers_changes_removed_from_chain)
    num_changes_unchanged = len(cur_read_change.kminmers_changes_unchanged)
    change_types[num_changes]['added'] += num_changes_added
    change_types[num_changes]['removed'] += num_changes_removed
    change_types[num_changes]['unchanged'] += num_changes_unchanged


change_types = dict(sorted(change_types.items()))
print(json.dumps(change_types, indent=4))

max_total_changes = 0
for num_changes in change_types.keys():
  total_changes = sum(change_types[num_changes].values())
  if total_changes > max_total_changes:
    max_total_changes = total_changes

###########
# Figure ##
###########

plt.style.use('BME163')
figure_size = (12, 6)
panel_size = (figure_size[0] - 0.75, figure_size[1] - 0.75)

panel_coor = (0.5, 0.5)
plt.figure(figsize=figure_size)
panel = plt.axes(
  [panel_coor[0] / figure_size[0],
   panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

max_num_changes = max(change_types.keys())
num_bars = max_num_changes

offset = (1 / num_bars) * 0.1
bar_width = (1 - ((num_bars - 1) * offset)) / num_bars
x_axis_ticks = []

x_position = 0
for num_changes in range(1, max_num_changes + 1):
  panel.add_patch(
    mplpatches.Rectangle(
      (x_position, 0),
      bar_width,
      change_types[num_changes]['added'],
      color=(89/255, 150/255, 210/255)
    )
  )

  panel.add_patch(
    mplpatches.Rectangle(
      (x_position, change_types[num_changes]['added']),
      bar_width,
      change_types[num_changes]['removed'],
      color=(210/255, 183/255, 100/255)
    )
  )

  panel.add_patch(
    mplpatches.Rectangle(
      (x_position, change_types[num_changes]['added'] + change_types[num_changes]['removed']),
      bar_width,
      change_types[num_changes]['unchanged'],
      color=(118/255, 98/255, 180/255)
    )
  )
  x_axis_ticks.append(x_position + bar_width / 2)
  x_position += bar_width + offset

panel.set_xlim(0, 1)
panel.set_ylim(0, max_total_changes + 100)

panel.set_xticks(x_axis_ticks)
panel.set_xticklabels(range(1, max_num_changes + 1))


plt.savefig(output_file, dpi=600)