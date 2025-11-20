import os
import io
import sys
import bte
import math
import argparse
import subprocess
import calculate_distance_score_utils as utils
from collections import defaultdict
from read_fasta import read_record

parser = argparse.ArgumentParser(description='Calculate distance score.')
parser.add_argument('--tree-in', required=True, help='Input tree file')
parser.add_argument('--true-in', required=True, help='True input file')
parser.add_argument('--estm-in', required=True, help='Estimated input file')
parser.add_argument('--tool', required=True, help='Tool name (MAMA or WEPP)')
parser.add_argument('--prefix', required=True, help='Output prefix')
parser.add_argument('--fasta-dir', required=False, help='Directory containing FASTA files for MAMA')
parser.add_argument('--gaps', action='store_true', help='Include gaps in distance calculation')

args = parser.parse_args()


true_abundance = {}
for line in open(args.true_in):
  fields = line.strip().split("\t")
  true_abundance[fields[0]] = float(fields[1])

estm_abundance = {}
amb_to_rep = {}
rep_to_amb = defaultdict(list)
for line in open(args.estm_in):
  fields = line.strip().split("\t")
  nodes = fields[0].split(",")
  estm_abundance[nodes[0]] = float(fields[1])
  for node in nodes:
    amb_to_rep[node] = nodes[0]
    rep_to_amb[nodes[0]].append(node)

cached_diff = {}
whds = {}
wpds = {}
if (args.tool == 'WEPP'):
  tree = bte.MATree(args.tree_in)
  for true_node, abundance in true_abundance.items():
    if true_node in amb_to_rep:
      whds[true_node] = [true_node, abundance, estm_abundance[amb_to_rep[true_node]], 0, 0.0]
    else:
      true_hap = tree.get_haplotype(true_node)
      min_node = None
      min_dist = math.inf

      for amb_node in amb_to_rep.keys():
        amb_hap = tree.get_haplotype(amb_node)
        diff = len(true_hap.symmetric_difference(amb_hap))
        cached_diff[true_node + amb_node] = diff
        if diff < min_dist:
          min_dist = diff
          min_node = amb_node

      abundance_diff = abs(abundance - estm_abundance[amb_to_rep[min_node]])
      whds[true_node] = [amb_to_rep[min_node], abundance, estm_abundance[amb_to_rep[min_node]], min_dist, min_dist * abundance]

  for estm_node, abundance in estm_abundance.items():
    ambiguous_nodes = rep_to_amb[estm_node]
    is_true = False
    is_true_node = None
    for node in ambiguous_nodes:
      if node in true_abundance:
        is_true = True
        is_true_node = node
        break
    if is_true:
      wpds[estm_node] = [is_true_node, abundance, true_abundance[is_true_node], 0, 0.0]
    else:
      min_node = None
      min_dist = math.inf
      for ambiguous_node in ambiguous_nodes:
        for true_node in true_abundance.keys():
          diff = 0
          if (true_node + ambiguous_node) in cached_diff:
            diff = cached_diff[true_node + ambiguous_node]
          else:
            true_hap = tree.get_haplotype(true_node)
            estm_hap = tree.get_haplotype(ambiguous_node)
            diff = len(true_hap.symmetric_difference(estm_hap))
            cached_diff[true_node + ambiguous_node] = diff
          if diff < min_dist:
            min_dist = diff
            min_node = true_node
      abundance_diff = abs(abundance - true_abundance[min_node])
      wpds[estm_node] = [min_node, abundance, true_abundance[min_node], min_dist, min_dist * abundance]
elif args.tool == 'MAMA':
  if not args.fasta_dir:
    print('Error: --fasta-dir must be provided for MAMA tool', file=sys.stderr)
    exit(1)
  for true_node, abundance in true_abundance.items():
    if true_node in amb_to_rep:
      whds[true_node] = [true_node, abundance, estm_abundance[amb_to_rep[true_node]], 0, 0.0]
    else:
      true_hap_fasta_path = utils.get_fasta_path(true_node, args.fasta_dir)
      min_node = None
      min_dist = math.inf

      for amb_node in amb_to_rep.keys():
        amb_hap_fasta_path = utils.get_fasta_path(amb_node, args.fasta_dir)
        mafft_out = subprocess.run(['sh', '-c', f'cat {true_hap_fasta_path} {amb_hap_fasta_path} | mafft --quiet --auto -'], capture_output=True, text=True)
        seqs = []
        for idn, seq in read_record(io.StringIO(mafft_out.stdout)):
          seqs.append((idn, seq))
        num_errors, num_snps, num_gaps, num_gaps_from_ends, num_base_to_ambiguous, match_string = utils.get_distance(seqs[0][1], seqs[1][1], get_match_string=False)
        cached_diff[true_node + amb_node] = num_snps
        if (num_snps < min_dist):
          min_dist = num_snps
          min_node = amb_node
      abundance_diff = abs(abundance - estm_abundance[amb_to_rep[min_node]])
      whds[true_node] = [amb_to_rep[min_node], abundance, estm_abundance[amb_to_rep[min_node]], min_dist, min_dist * abundance]

  for estm_node, abundance in estm_abundance.items():
    ambiguous_nodes = rep_to_amb[estm_node]
    is_true = False
    is_true_node = None
    for node in ambiguous_nodes:
      if node in true_abundance:
        is_true = True
        is_true_node = node
        break
    if is_true:
      wpds[estm_node] = [is_true_node, abundance, true_abundance[is_true_node], 0, 0.0]
    else:
      min_node = None
      min_dist = math.inf
      for ambiguous_node in ambiguous_nodes:
        for true_node in true_abundance.keys():
          diff = 0
          if (true_node + ambiguous_node) in cached_diff:
            diff = cached_diff[true_node + ambiguous_node]
          else:
            true_hap_fasta_path = utils.get_fasta_path(true_node, args.fasta_dir)
            amb_hap_fasta_path = utils.get_fasta_path(ambiguous_node, args.fasta_dir)
            mafft_out = subprocess.run(['sh', '-c', f'cat {true_hap_fasta_path} {amb_hap_fasta_path} | mafft --quiet --auto -'], capture_output=True, text=True)
            seqs = []
            for idn, seq in read_record(io.StringIO(mafft_out.stdout)):
              seqs.append((idn, seq))
            num_errors, num_snps, num_gaps, num_gaps_from_ends, num_base_to_ambiguous, match_string = utils.get_distance(seqs[0][1], seqs[1][1], get_match_string=False)
            diff = num_snps
            cached_diff[true_node + ambiguous_node] = num_snps
          if diff < min_dist:
            min_dist = diff
            min_node = true_node
      abundance_diff = abs(abundance - true_abundance[min_node])
      wpds[estm_node] = [min_node, abundance, true_abundance[min_node], min_dist, min_dist * abundance]
else:
  print(f'Error: Unknown tool {args.tool}.. (must be either WEPP or MAMA)', file=sys.stderr)
  exit(1)

whds_out = args.prefix + '.wepp.whd.tsv' if args.tool == 'WEPP' else args.prefix + '.mama.whd.tsv'
with open(whds_out, 'w') as whd_fh:
  whd_fh.write(f'True_node\tClosest_estimate\tTrue_abundance\tEstimated_abundance\tDistance_to_closest\tWeighted_haplotype_dist\n')
  for true_node, distance in whds.items():
    whd_fh.write(f'{true_node}\t{distance[0]}\t{distance[1]}\t{distance[2]}\t{distance[3]}\t{distance[4]}\n')

# whds_diffs_out = args.prefix + '.wepp.whd_diffs.tsv' if args.tool == 'WEPP' else args.prefix + '.mama.whd_diffs.tsv'
# with open(whds_diffs_out, 'w') as whd_diff_fh:
#   whd_diff_fh.write(f'True_node\tClosest_estimate\tDiffs\n')
#   for true_node, distance in whds.items():
#     estm_node = distance[0]
#     whd_diff_fh.write(f'{true_node}\t{estm_node}\t')
#     if (distance[3] > 0):
#       true_hap = tree.get_haplotype(true_node)
#       estm_hap = tree.get_haplotype(estm_node)
#       diffs = list(true_hap.symmetric_difference(estm_hap))
#       diffs.sort(key=lambda x: int(x[1:-1]))
#       for i in range(len(diffs)):
#         whd_diff_fh.write(f'{diffs[i]}')
#         if i < len(diffs) - 1: whd_diff_fh.write(',')
#         else: whd_diff_fh.write('\n')
#     else:
#       whd_diff_fh.write('.\n')

wpds_out = args.prefix + '.wepp.wpd.tsv' if args.tool == 'WEPP' else args.prefix + '.mama.wpd.tsv'
with open(wpds_out, 'w') as wpd_fh:
  wpd_fh.write(f'Estimated_node\tClosest_true\tEstimated_abundance\tTrue_abundance\tDistance_to_closest\tWeighted_peak_dist\n')
  for estm_node, distance in wpds.items():
    wpd_fh.write(f'{estm_node}\t{distance[0]}\t{distance[1]}\t{distance[2]}\t{distance[3]}\t{distance[4]}\n')

sum_distance_out = args.prefix + '.wepp.distance_sum.tsv' if args.tool == 'WEPP' else args.prefix + '.mama.distance_sum.tsv'
with open(sum_distance_out, 'w') as sum_fh:
  total_whd = sum([distance[4] for distance in whds.values()])
  total_wpd = sum([distance[4] for distance in wpds.values()])
  sum_fh.write(f'Tool\tTotal_Weighted_Haplotype_Distance\tTotal_Weighted_Peak_Distance\n')
  sum_fh.write(f'WEPP\t{total_whd}\t{total_wpd}\n')

