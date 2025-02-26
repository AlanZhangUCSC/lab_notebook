from collections import defaultdict
from itertools import product
import numpy as np
import scipy.stats as stats
import subprocess
import argparse
import gzip
import math
import sys
import os
import heuristic_v5_utils as utils
import json

parser = argparse.ArgumentParser()
parser.add_argument('--vcfs', type=str, nargs='+', help='Paths to the allvar VCF files')
parser.add_argument('--bams', type=str, nargs='+', help='Paths to the BAM files')
parser.add_argument('--pcfs', type=str, nargs='+', help='Paths to the PCF files')
parser.add_argument('--outputs', type=str, nargs='+', help='Paths to the output files')
parser.add_argument('--abundances', type=str, help='Path to the haplotype abundance file')
parser.add_argument('--assignments', type=str, help='Path to reference assignments')
parser.add_argument('--fasta', type=str, help='Path to the aligned FASTA file')
parser.add_argument('--min_coverage', type=int, default=10, help='Minimum coverage to call a variant')
parser.add_argument('--min_allele_depth', type=int, default=5, help='Minimum allele depth to call a variant')
parser.add_argument('--min_coverage_per_ref_assignment_distribution', type=int, default=10, help='Minimum coverage to call a variant per reference assignment distribution')
parser.add_argument('--min_allele_depth_per_ref_assignment_distribution', type=int, default=2, help='Minimum allele depth to call a variant per reference assignment distribution')
parser.add_argument('--error_rate', type=float, default=0.05, help='Sequencing error rate')
parser.add_argument('--confidence_threshold', type=float, default=2, help='Confidence threshold for using statistical significance')
parser.add_argument('--debug', type=str, nargs='+', help='Debug mode')
args = parser.parse_args()

vcf_paths = args.vcfs
bam_paths = args.bams
pcf_paths = args.pcfs
output_paths = args.outputs
abundance_file_path = args.abundances
aligned_fasta_path = args.fasta
reference_assignments_path = args.assignments
error_rate = args.error_rate
confidence_threshold = args.confidence_threshold
debug = False
debug_positions = set()
debug_haplotype = None
if args.debug:
  debug = True
  debug_haplotype = args.debug[0]
  debug_positions = set(int(pos) for pos in args.debug[1:])


# Load haplotype abundance
haplotype_abundance = {}
with utils.open_file(abundance_file_path) as fh:
  for line in fh:
    line = line.strip()
    if not line: continue
    if line.startswith('#'): continue
    fields = line.split('\t')
    haplotype = fields[0].split(',')[0]
    abundance = float(fields[1])
    haplotype_abundance[haplotype] = abundance

# Sort haplotypes by abundance in descending order and assign binary codes
haplotype_to_binary = {}
binary_to_haplotype = {}
sorted_haplotypes = sorted(haplotype_abundance.items(), key=lambda x: x[1])
for i, (haplotype, _) in enumerate(sorted_haplotypes):
  haplotype_to_binary[haplotype] = 2**i
  binary_to_haplotype[2**i] = haplotype

# Load reference assignments
reference_assignments = {}
with utils.open_file(reference_assignments_path) as fh:
  for i, line in enumerate(fh):
    line = line.strip()
    if not line: continue
    fields = line.split('\t')
    query_name = fields[0]
    reference_info_list = fields[1].strip().split(' ')
    reference_assignments[query_name] = []
    for reference_info in reference_info_list:
      reference, score = reference_info.split(',')
      if reference in haplotype_to_binary:
        reference_assignments[query_name].append(utils.assigned_reference(haplotype_to_binary[reference], float(score)))

# Load pcf info
pcf_info = defaultdict(lambda: defaultdict(list))
for pcf_path in pcf_paths:
  with utils.open_file(pcf_path) as fh:
    for line in fh:
      line = line.strip()
      if not line: continue
      fields = line.split('\t')
      haplotype, pos = fields[0], int(fields[1])
      pcf_info[haplotype][pos] = line


resolved_alleles_by_haplotype = defaultdict(lambda: defaultdict(str))
unresolved_sites_by_haplotype = defaultdict(list)
out_str_by_haplotype = defaultdict(list)
position_info_by_haplotype_dict = {}
files_by_haplotype = defaultdict(list)
ref_seqs = {}
fhs = []
for vcf_path in args.vcfs: fhs.append(utils.open_file(vcf_path))

for i in range(len(fhs)):
  vcf_fh, bam_path, output_path = fhs[i], bam_paths[i], output_paths[i]
  print(f'Processing {vcf_fh.name}', file=sys.stderr)
  cur_header_lines = []
  for vcf_line in vcf_fh:
    vcf_line = vcf_line.strip()
    if vcf_line.startswith('#'):
      cur_header_lines.append(vcf_line)
      continue
    
    fields = vcf_line.split('\t')
    haplotype, pos, indel = fields[0], int(fields[1]), fields[7].startswith('INDEL')
    if out_str_by_haplotype[haplotype] == []:
      out_str_by_haplotype[haplotype] += cur_header_lines
    if haplotype not in files_by_haplotype:
      files_by_haplotype[haplotype] = [vcf_fh, bam_path, output_path]
    if haplotype not in position_info_by_haplotype_dict:
      position_info_by_haplotype, ref_seqs = utils.process_aligned_fasta(aligned_fasta_path, haplotype)
      position_info_by_haplotype_dict[haplotype] = position_info_by_haplotype

    if pos not in pcf_info[haplotype]: continue
    haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, candidate_alleles, ref_assignment_distributions = utils.parse_pcf_line(pcf_info[haplotype][pos])
    alleles = [ref_allele] + alt_alleles

    if debug and haplotype == debug_haplotype and pos in debug_positions: print('-' * 100)
    
    if len(candidate_alleles) == 1:
      gt = alleles.index(list(candidate_alleles)[0])
      assert(gt != 0)
      if debug and haplotype == debug_haplotype:
        print('\t'.join(fields[:-1]), end='')
        print(f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      out_str_by_haplotype[haplotype].append('\t'.join(fields[:-1]) + f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      continue
    elif indel:
      if debug and haplotype == debug_haplotype:
        print(vcf_line)
      out_str_by_haplotype[haplotype].append(vcf_line)
      continue
    
    debug_line = debug and haplotype == debug_haplotype and pos in debug_positions
    cross_referenced_alleles_by_ref_assignment_distributions = utils.cross_reference_alleles_2(
      pos, haplotype, position_info_by_haplotype_dict[haplotype], ref_seqs, ref_assignment_distributions, 
      pcf_info, candidate_alleles, haplotype_abundance, binary_to_haplotype, haplotype_to_binary, 
      args.min_coverage_per_ref_assignment_distribution, args.min_allele_depth_per_ref_assignment_distribution, error_rate,
      debug_line
    )

    
    min_max_score_diff = min(ref_assignment_distributions.keys())
    alleles_support_count = defaultdict(int)
    for haplotypes_binary, cross_referenced_alleles in cross_referenced_alleles_by_ref_assignment_distributions[min_max_score_diff].items():
      top_allele = cross_referenced_alleles[0][0]
      alleles_support_count[top_allele] += 1
    alleles_support_count = dict(sorted(alleles_support_count.items(), key=lambda x: x[1], reverse=True))

    log_likelihood_ratio_by_num_haplotypes_assigned = {}
    for haplotypes_binary, cross_referenced_alleles in cross_referenced_alleles_by_ref_assignment_distributions[min_max_score_diff].items():
      num_haplotypes_assigned = bin(haplotypes_binary).count('1')
      alleles_total_count = ref_assignment_distributions[min_max_score_diff][haplotypes_binary].get_passed_total_count()
      top_allele = cross_referenced_alleles[0][0]
      if len(cross_referenced_alleles) == 1:
        if num_haplotypes_assigned not in log_likelihood_ratio_by_num_haplotypes_assigned:
          log_likelihood_ratio_by_num_haplotypes_assigned[num_haplotypes_assigned] = (top_allele, 0)
      else:
        log_likelihood_ratio_normalized = (cross_referenced_alleles[0][1] - cross_referenced_alleles[1][1]) / math.log(alleles_total_count)
        assert(log_likelihood_ratio_normalized >= 0)
        if num_haplotypes_assigned not in log_likelihood_ratio_by_num_haplotypes_assigned:
          log_likelihood_ratio_by_num_haplotypes_assigned[num_haplotypes_assigned] = (top_allele, log_likelihood_ratio_normalized)
        else:
          if log_likelihood_ratio_normalized > log_likelihood_ratio_by_num_haplotypes_assigned[num_haplotypes_assigned][1]:
            log_likelihood_ratio_by_num_haplotypes_assigned[num_haplotypes_assigned] = (top_allele, log_likelihood_ratio_normalized)


    if debug and haplotype == debug_haplotype and pos in debug_positions:
      print(vcf_line)
      print(haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, candidate_alleles)
      print('cross_referenced_alleles_by_ref_assignment_distributions')
      for max_score_diff in cross_referenced_alleles_by_ref_assignment_distributions:
        for haplotypes_binary, cross_referenced_alleles in cross_referenced_alleles_by_ref_assignment_distributions[max_score_diff].items():
          print(f"{bin(haplotypes_binary).count('1')} {haplotypes_binary:0{len(haplotype_to_binary)}b}", cross_referenced_alleles)
      print()

      print('log_likelihood_ratio_by_num_haplotypes_assigned')
      for num_haplotypes_assigned in sorted(log_likelihood_ratio_by_num_haplotypes_assigned.keys()):
        print(f'{num_haplotypes_assigned}: {log_likelihood_ratio_by_num_haplotypes_assigned[num_haplotypes_assigned]}')
      print()

      print('ref_assignment_distributions')
      for max_score_diff in ref_assignment_distributions:
        print(f'max_score_diff: {max_score_diff}')
        for haplotypes_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff].items():
          print(f"{haplotypes_binary:0{len(haplotype_to_binary)}b}", 
            dict(ref_assignment_distribution.alleles_count),
            dict(ref_assignment_distribution.unpassed_alleles_count),
            dict(ref_assignment_distribution.all_alleles_count), sep='\t', end='')
          print()
        print()
      
      print('alleles_support_count')
      print(alleles_support_count)
      print()

    if len(alleles_support_count) == 1:
      allele = list(alleles_support_count.keys())[0]
      resolved_alleles_by_haplotype[haplotype][pos] = allele
      gt = alleles.index(allele)
      if gt == 0: continue
      if debug and haplotype == debug_haplotype:
        print('\t'.join(fields[:-1]), end='')
        print(f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      out_str_by_haplotype[haplotype].append('\t'.join(fields[:-1]) + f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      continue
    elif len(alleles_support_count) >= 2:
      total_likelihood_ratio_by_allele = defaultdict(float)
      for num_haplotypes_assigned, (top_allele, log_likelihood_ratio_normalized) in log_likelihood_ratio_by_num_haplotypes_assigned.items():
        total_likelihood_ratio_by_allele[top_allele] += log_likelihood_ratio_normalized
      sorted_total_likelihood_ratio_by_allele = sorted(total_likelihood_ratio_by_allele.items(), key=lambda x: x[1], reverse=True)
      if len(sorted_total_likelihood_ratio_by_allele) == 1 or sorted_total_likelihood_ratio_by_allele[0][1] - sorted_total_likelihood_ratio_by_allele[1][1] >= confidence_threshold:
        allele = sorted_total_likelihood_ratio_by_allele[0][0]
        resolved_alleles_by_haplotype[haplotype][pos] = allele
        gt = alleles.index(allele)
        if gt == 0: continue
        if debug and haplotype == debug_haplotype:
          print('\t'.join(fields[:-1]), end='')
          print(f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
        out_str_by_haplotype[haplotype].append('\t'.join(fields[:-1]) + f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
        continue
    
    unresolved_sites_by_haplotype[haplotype].append((pos, set(alleles_support_count.keys())))

    

# for fh in fhs: fh.seek(0)
# for haplotype, unresolved_sites in unresolved_sites_by_haplotype.items():
#   vcf_fh, bam_path, output_path = files_by_haplotype[haplotype]
#   position_info_by_haplotype = position_info_by_haplotype_dict[haplotype]
#   for pos, candidate_alleles in unresolved_sites:
#     print(f'Processing unresolved site {haplotype} {pos}', file=sys.stderr)
#     candidate_alleles_by_haplotype = defaultdict(set)
#     haplotype_, pos_, ref_allele, alt_alleles, ref_depth, alt_depths, candidate_alleles, ref_assignment_distributions = utils.parse_pcf_line(pcf_info[haplotype][pos])
#     assert(haplotype == haplotype_ and pos == pos_)
#     alignable_haplotypes = utils.get_alignable_haplotypes(ref_assignment_distributions, binary_to_haplotype)
#     candidate_alleles_by_haplotype[haplotype] = set(candidate_alleles)
#     alleles = set(candidate_alleles)
#     base_by_query_name = defaultdict(str)
#     conflicting_base_query = set()
#     for i, cur_other_haplotype in enumerate(alignable_haplotypes):
#       curother_ref_allele, curother_pos = None, None
#       if cur_other_haplotype == haplotype:
#         curother_ref_allele, curother_pos = ref_allele, pos
#       else:
#         curother_ref_allele, curother_pos = position_info_by_haplotype[cur_other_haplotype][pos-1]

#       if cur_other_haplotype != haplotype:
#         if curother_ref_allele == '-': continue
#         if curother_pos in resolved_alleles_by_haplotype[cur_other_haplotype]:
#           candidate_alleles_by_haplotype[cur_other_haplotype] = set([resolved_alleles_by_haplotype[cur_other_haplotype][curother_pos]])
#           alleles.add(resolved_alleles_by_haplotype[cur_other_haplotype][curother_pos])
#         else:
#           if curother_pos in pcf_info[cur_other_haplotype]:
#             cur_other_haplotype, curother_pos, curother_ref_allele, curother_alt_alleles, curother_ref_depth, curother_alt_depths, curother_candidate_alleles, curother_ref_assignment_distributions = utils.parse_pcf_line(pcf_info[cur_other_haplotype][curother_pos])
#             candidate_alleles_by_haplotype[cur_other_haplotype] = set(curother_candidate_alleles)
#             for allele in curother_candidate_alleles:
#               alleles.add(allele)
#           else:
#             if curother_ref_allele in ('A', 'C', 'G', 'T'):
#               candidate_alleles_by_haplotype[cur_other_haplotype] = set([curother_ref_allele])
#               alleles.add(curother_ref_allele)
#             else:
#               candidate_alleles_by_haplotype[cur_other_haplotype] = set(candidate_alleles)
      
#       curother_vcf_fh, curother_bam_path = fhs[i], bam_paths[i]
#       sam_lines_at_pos = utils.get_sam_lines_at_pos(curother_bam_path, cur_other_haplotype, curother_pos)
#       for sam_line in sam_lines_at_pos:
#         qname, template_info, base, mapq, pass_filter = utils.process_sam_line(sam_line, curother_pos, False, alleles)
#         print(qname, template_info, base, mapq, pass_filter)
#         if base not in alleles or qname is None: continue
#         if qname not in conflicting_base_query:
#           if qname in base_by_query_name:
#             if pass_filter and mapq == 60 and base_by_query_name[qname] != base:
#               conflicting_base_query.add(qname)
#               base_by_query_name.pop(qname)
#           else:
#             if pass_filter and mapq == 60:
#               print('adding to base_by_query_name', qname, base)
#               base_by_query_name[qname] = base
      
#     all_alleles_count = defaultdict(int)
#     for query_name, base in base_by_query_name.items():
#       if base in ('A', 'T', 'C', 'G') and query_name not in conflicting_base_query:
#         all_alleles_count[base] += 1
#     total_alleles_count = sum(all_alleles_count.values())
#     all_observed_allele_frequency = {allele: count / total_alleles_count for allele, count in all_alleles_count.items()}

#     allele_combinations_dict, kl_divergences, allele_combinations_expected_allele_frequency = utils.calculate_allele_combination_kl_divergence(candidate_alleles_by_haplotype, haplotype_abundance, all_observed_allele_frequency)

#     if debug and haplotype == debug_haplotype and pos in debug_positions:
#       print()
#       for cur_haplotype, cur_candidate_alleles in candidate_alleles_by_haplotype.items():
#         print(f'{cur_haplotype} candidate alleles: {cur_candidate_alleles}')
#       print(list(candidate_alleles_by_haplotype.keys()))
#       print(f'all_alleles_count: {dict(sorted((k, v) for k, v in all_alleles_count.items()))}')
#       print(f'all_observed_allele_frequency: {dict(sorted((k, round(v, 3)) for k, v in all_observed_allele_frequency.items()))}')
#       for allele_combination_dict, kl_divergence, expected_allele_frequency in zip(allele_combinations_dict, kl_divergences, allele_combinations_expected_allele_frequency):
#         for _, allele in allele_combination_dict.items():
#           print(allele, end='')
#         print(f'\t{round(kl_divergence, 5)}\t{dict(sorted((k, round(v, 3)) for k, v in expected_allele_frequency.items()))}')
#       print()

#     minimum_kl_divergence_index = kl_divergences.index(min(kl_divergences))
#     minimum_kl_divergence_allele = allele_combinations_dict[minimum_kl_divergence_index][haplotype]
#     gt = ([ref_allele] + alt_alleles).index(minimum_kl_divergence_allele)
#     if gt == 0: continue
#     if debug and haplotype == debug_haplotype:
#       print('\t'.join(fields[:-1]), end='')
#       print(f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
#     out_str_by_haplotype[haplotype].append('\t'.join(fields[:-1]) + f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))

if not debug:
  for haplotype, out_str in out_str_by_haplotype.items():
    output_path = files_by_haplotype[haplotype][2]
    with open(output_path, 'w') as fh:
      for line in out_str:
        fh.write(line + '\n')
