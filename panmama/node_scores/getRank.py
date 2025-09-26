import sys
import math

true_haplotypes_file = sys.argv[1]
node_scores_file = sys.argv[2]
node_scores_sorted_by_node_scores_file = sys.argv[3]


class true_haplotype:
  def __init__(self, identifier, abundance):
    self.identifier = identifier
    self.abundance = abundance

true_haplotypes = []
with open(true_haplotypes_file, 'r') as f:
  for line in f:
    line = line.strip()
    identifier, abundance = line.split('\t')
    true_haplotypes.append(true_haplotype(identifier, abundance))

class selected_node:
  def __init__(self, node_id, kminmer_overlap, kminmer_overlap_rank, node_score, node_score_rank, node_score_rank_among_selected_nodes):
    self.node_id = node_id
    self.kminmer_overlap = kminmer_overlap
    self.kminmer_overlap_rank = kminmer_overlap_rank
    self.node_score = node_score
    self.node_score_rank = node_score_rank
    self.node_score_rank_among_selected_nodes = node_score_rank_among_selected_nodes


selected_nodes = {}
cur_kminmer_overlap = 999.0
cur_kminmer_overlap_rank = -1
minimum_kminmer_overlap_to_be_selected = -1
with open(node_scores_file, 'r') as f:
  for line in f:
    line = line.strip()
    node_id, kminmer_overlap, node_score = line.split(' ')
    kminmer_overlap = float(kminmer_overlap)
    node_score = float(node_score)
    if kminmer_overlap < cur_kminmer_overlap:
      cur_kminmer_overlap = kminmer_overlap
      cur_kminmer_overlap_rank += 1
      if cur_kminmer_overlap_rank == 1000: minimum_kminmer_overlap_to_be_selected = kminmer_overlap
    selected_nodes[node_id] = selected_node(node_id, kminmer_overlap, cur_kminmer_overlap_rank, node_score, -1, -1)
if minimum_kminmer_overlap_to_be_selected == -1: minimum_kminmer_overlap_to_be_selected = cur_kminmer_overlap

cur_node_score = float('inf')
cur_node_score_rank = -1
cur_node_score_rank_among_selected_nodes = -1
with open(node_scores_sorted_by_node_scores_file, 'r') as f:
  for line in f:
    line = line.strip()
    node_id, kminmer_overlap, node_score = line.split(' ')
    kminmer_overlap = float(kminmer_overlap)
    node_score = float(node_score)
    if node_score < cur_node_score:
      cur_node_score = node_score
      cur_node_score_rank += 1
      if kminmer_overlap >= minimum_kminmer_overlap_to_be_selected:
        cur_node_score_rank_among_selected_nodes += 1
    selected_nodes[node_id].node_score_rank = cur_node_score_rank
    selected_nodes[node_id].node_score_rank_among_selected_nodes = cur_node_score_rank_among_selected_nodes

for true_haplotype in true_haplotypes:
  if true_haplotype.identifier in selected_nodes:
    num_no_less_than_kminmer_overlap = sum(1 for node in selected_nodes.values() if node.kminmer_overlap >= selected_nodes[true_haplotype.identifier].kminmer_overlap)
    num_no_less_than_node_score = sum(1 for node in selected_nodes.values() if node.node_score >= selected_nodes[true_haplotype.identifier].node_score)
    num_no_less_than_node_score_among_selected_nodes = sum(1 for node in selected_nodes.values()
                                                          if node.node_score >= selected_nodes[true_haplotype.identifier].node_score
                                                          and node.kminmer_overlap >= minimum_kminmer_overlap_to_be_selected)
    print(true_haplotype.identifier,
          true_haplotype.abundance,
          selected_nodes[true_haplotype.identifier].kminmer_overlap,
          selected_nodes[true_haplotype.identifier].kminmer_overlap_rank,
          num_no_less_than_kminmer_overlap,
          selected_nodes[true_haplotype.identifier].node_score,
          selected_nodes[true_haplotype.identifier].node_score_rank,
          num_no_less_than_node_score,
          selected_nodes[true_haplotype.identifier].node_score_rank_among_selected_nodes,
          num_no_less_than_node_score_among_selected_nodes,
          sep='\t')
  else:
    print(true_haplotype.identifier,
          true_haplotype.abundance,
          "null",
          "null",
          "null",
          "null",
          "null",
          "null",
          "null",
          sep='\t')

