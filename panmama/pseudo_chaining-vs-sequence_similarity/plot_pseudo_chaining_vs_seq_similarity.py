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

sam_file = sys.argv[1]
chain_file = sys.argv[2]

def get_num_matches_from_md_tag(md_tag):
  numbers = re.findall(r'\d+', md_tag)
  return sum(map(int, numbers))

sam_matches = defaultdict(lambda: defaultdict(int))
with open(sam_file, 'r') as f:
  while True:
    line = f.readline()
    if not line: break
    fields = line.strip().split('\t')
    read_name = fields[0]
    idx, snps = map(int, read_name.split('_')[1:])
    ref_name = fields[2]
    md_tag = fields[-3]

    num_matches = 0 if ref_name == '*' else get_num_matches_from_md_tag(md_tag)
    if not sam_matches[idx][snps]:
      sam_matches[idx][snps] = num_matches
    else:
      if num_matches > sam_matches[idx][snps]:
        sam_matches[idx][snps] = num_matches

chain_matches = defaultdict(lambda: defaultdict(tuple))
with open(chain_file, 'r') as f:
  for line in f:
    fields = line.strip().split('\t')
    idx, snps = map(int, fields[0].split('_')[1:])
    chain_matches[idx][snps] = (int(fields[1]), int(fields[2]))

total_chain_entries = sum(len(snps) for snps in chain_matches.values())
total_sam_entries = sum(len(snps) for snps in sam_matches.values())

chain_scores = []
sam_scores = []
for idx in chain_matches.keys():
  for snps in chain_matches[idx].keys():
    if snps > 10: continue
    chain_scores.append(chain_matches[idx][snps][1] / chain_matches[idx][snps][0])
    sam_scores.append(sam_matches[idx][snps] / 150)

############
## Figure ##
############

plt.style.use('BME163')
figure_size = (6, 6)
panel_size = (figure_size[0] - 1, figure_size[1] - 1)
panel_coor = (0.5, 0.5)

plt.figure(figsize=figure_size)
panel = plt.axes(
  [panel_coor[0] / figure_size[0],
   panel_coor[1] / figure_size[1],
   panel_size[0] / figure_size[0],
   panel_size[1] / figure_size[1]]
)

panel.set_xlabel('Sequence similarity (matches / sequence length)')
panel.set_ylabel('Pseudo-chaining score / max possible chain score')
panel.scatter(sam_scores, chain_scores, s=0.5)

plt.savefig('pseudo_chaining_vs_seq_similarity.png', dpi=600)





    




