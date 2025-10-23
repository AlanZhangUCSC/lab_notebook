import sys
from collections import defaultdict
from itertools import combinations

node_scores_path = sys.argv[1]

identical_score_groups = defaultdict(lambda: defaultdict(list))

fh = open(node_scores_path, "r")
fh.readline()
for line in fh:
  fields = line.strip().split("\t")
  node_id = fields[0]
  sumReadScore = int(fields[2])
  sumMPReadScore = int(fields[7])
  identical_score_groups[sumReadScore][sumMPReadScore].append(node_id)
fh.close()

potential_identical_nodes = []
for sumReadScore in identical_score_groups:
  for sumMPReadScore in identical_score_groups[sumReadScore]:
    node_ids = identical_score_groups[sumReadScore][sumMPReadScore]
    if len(node_ids) > 1:
      potential_identical_nodes.append(node_ids)

potential_identical_nodes.sort(key=len, reverse=True)
for node_ids in potential_identical_nodes:
  for node_id1, node_id2 in combinations(node_ids, 2):
    print(f"{node_id1}\t{node_id2}")
  