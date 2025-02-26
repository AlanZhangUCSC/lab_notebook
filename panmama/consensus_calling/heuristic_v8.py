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
parser.add_argument('--min_allele_depth', type=str, default='5,0.02', help='Minimum allele raw depth and ratio to call a variant (the highest of the two is taken)')
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


resolved_alleles_by_haplotype = defaultdict(lambda: defaultdict(list))
unresolved_sites_by_haplotype = defaultdict(lambda: defaultdict(list))
unresolved_site_info_by_haplotype = defaultdict(lambda: defaultdict(list))
out_str_by_haplotype = defaultdict(lambda: defaultdict(str))
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
    if out_str_by_haplotype[haplotype]['headers'] == '':
      out_str_by_haplotype[haplotype]['headers'] = '\n'.join(cur_header_lines)
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
      resolved_alleles_by_haplotype[haplotype][pos] = [list(candidate_alleles)[0], []]
      if debug and haplotype == debug_haplotype:
        print('\t'.join(fields[:-1]), end='')
        print(f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      out_str_by_haplotype[haplotype][pos] = ('\t'.join(fields[:-1]) + f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      continue
    elif indel:
      if debug and haplotype == debug_haplotype:
        print(vcf_line)
      out_str_by_haplotype[haplotype][pos] = vcf_line
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
    log_likelihood_ratio_by_read_groups = {}
    for haplotypes_binary, cross_referenced_alleles in cross_referenced_alleles_by_ref_assignment_distributions[min_max_score_diff].items():
      num_haplotypes_assigned = bin(haplotypes_binary).count('1')
      alleles_total_count = ref_assignment_distributions[min_max_score_diff][haplotypes_binary].get_passed_total_count()
      top_allele = cross_referenced_alleles[0][0]
      if len(cross_referenced_alleles) == 1:
        log_likelihood_ratio_by_read_groups[haplotypes_binary] = (top_allele, 0)
        if num_haplotypes_assigned not in log_likelihood_ratio_by_num_haplotypes_assigned:
          log_likelihood_ratio_by_num_haplotypes_assigned[num_haplotypes_assigned] = (top_allele, 0)
      else:
        log_likelihood_ratio_normalized = (cross_referenced_alleles[0][1] - cross_referenced_alleles[1][1]) / math.log(alleles_total_count)
        assert(log_likelihood_ratio_normalized >= 0)
        log_likelihood_ratio_by_read_groups[haplotypes_binary] = (top_allele, log_likelihood_ratio_normalized)
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

      print('log_likelihood_ratio_by_read_groups')
      for max_score_diff in cross_referenced_alleles_by_ref_assignment_distributions:
        for haplotypes_binary, (top_allele, log_likelihood_ratio_normalized) in log_likelihood_ratio_by_read_groups.items():
          print(f"{bin(haplotypes_binary).count('1')} {haplotypes_binary:0{len(haplotype_to_binary)}b}", top_allele, log_likelihood_ratio_normalized)
      print()

      print('log_likelihood_ratio_by_num_haplotypes_assigned')
      for num_haplotypes_assigned, (top_allele, log_likelihood_ratio_normalized) in sorted(log_likelihood_ratio_by_num_haplotypes_assigned.items()):
        print(f'{num_haplotypes_assigned}: {top_allele} {log_likelihood_ratio_normalized}')
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

    
    total_likelihood_ratio_by_allele = defaultdict(float)
    for read_group, (top_allele, log_likelihood_ratio_normalized) in log_likelihood_ratio_by_read_groups.items():
      total_likelihood_ratio_by_allele[top_allele] += log_likelihood_ratio_normalized
    sorted_total_likelihood_ratio_by_allele = sorted(total_likelihood_ratio_by_allele.items(), key=lambda x: x[1], reverse=True)
    if debug and haplotype == debug_haplotype and pos in debug_positions:
      print(sorted_total_likelihood_ratio_by_allele)
      print()
    
    if len(sorted_total_likelihood_ratio_by_allele) <= 1:
      try:
        likelihood_ratio_difference = sorted_total_likelihood_ratio_by_allele[0][1] 
      except:
        likelihood_ratio_difference = -1
    else:
      likelihood_ratio_difference = sorted_total_likelihood_ratio_by_allele[0][1] - sorted_total_likelihood_ratio_by_allele[1][1]
    if likelihood_ratio_difference >= confidence_threshold:
      allele = sorted_total_likelihood_ratio_by_allele[0][0]
      resolved_alleles_by_haplotype[haplotype][pos] = [allele, sorted_total_likelihood_ratio_by_allele]
      gt = alleles.index(allele)
      if gt == 0: continue
      if debug and haplotype == debug_haplotype:
        print('\t'.join(fields[:-1]), end='')
        print(f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:]))
      out_str_by_haplotype[haplotype][pos] = '\t'.join(fields[:-1]) + f'\t{gt}:' + ':'.join(fields[-1].split(':')[1:])
      continue


    
    unresolved_sites_by_haplotype[haplotype][pos] = [sorted_total_likelihood_ratio_by_allele, utils.get_alignable_haplotypes(ref_assignment_distributions, binary_to_haplotype)]
    unresolved_site_info_by_haplotype[haplotype][pos] = [alleles, fields]


for cur_haplotype in unresolved_sites_by_haplotype:
  for cur_pos in unresolved_sites_by_haplotype[cur_haplotype]:
    # if already resolved -> continue
    if cur_pos in resolved_alleles_by_haplotype[cur_haplotype]:
      continue
    if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
      print('-' * 100)
      print(cur_haplotype, cur_pos)

    # check if there are any cryptic mutations (allele not seen in the references)
    ref_alleles = set()
    observed_alleles = set()
    ref_by_pcf = set()
    liklihood_ratio_by_haplotype = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    cur_likelihood_ratio_by_allele, alignable_haplotypes = unresolved_sites_by_haplotype[cur_haplotype][cur_pos]
    for alignable_haplotype in alignable_haplotypes:
      alignable_haplotype_ref, alignable_haplotype_pos = position_info_by_haplotype_dict[cur_haplotype][alignable_haplotype][cur_pos-1]
      alignable_haplotype_pos += 1
      ref_alleles.add(alignable_haplotype_ref)
      if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
        print(alignable_haplotype, alignable_haplotype_pos, alignable_haplotype_ref)
      if alignable_haplotype in resolved_alleles_by_haplotype and alignable_haplotype_pos in resolved_alleles_by_haplotype[alignable_haplotype]:
        observed_alleles.add(resolved_alleles_by_haplotype[alignable_haplotype][alignable_haplotype_pos][0])
        for allele, likelihood_ratio in resolved_alleles_by_haplotype[alignable_haplotype][alignable_haplotype_pos][1]:
          liklihood_ratio_by_haplotype[alignable_haplotype][allele] = likelihood_ratio
          observed_alleles.add(allele)
      elif alignable_haplotype in unresolved_sites_by_haplotype and alignable_haplotype_pos in unresolved_sites_by_haplotype[alignable_haplotype]:
        for allele, likelihood_ratio in unresolved_sites_by_haplotype[alignable_haplotype][alignable_haplotype_pos][0]:
          liklihood_ratio_by_haplotype[alignable_haplotype][allele] = likelihood_ratio
          observed_alleles.add(allele)
      else:
        ref_by_pcf.add(alignable_haplotype)
        observed_alleles.add(alignable_haplotype_ref)

    for haplotype in liklihood_ratio_by_haplotype:
      for allele, likelihood_ratio in liklihood_ratio_by_haplotype[haplotype].items():
        observed_alleles.add(allele)

    cryptic_allele = set()
    for allele in observed_alleles:
      if allele not in ref_alleles:
        cryptic_allele.add(allele)

    if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
      print()
      print('ref_alleles:', ref_alleles)
      print('observed_alleles:', observed_alleles)
      print('cryptic_allele:', cryptic_allele)
      print()
        
    # if there is no cryptic mutation -> assign reference allele
    if not cryptic_allele:
      if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
        print('no cryptic allele, assigning reference allele\n')
      ref_allele = position_info_by_haplotype_dict[cur_haplotype][cur_haplotype][cur_pos-1][0]
      resolved_alleles_by_haplotype[cur_haplotype][cur_pos] = [ref_allele, unresolved_sites_by_haplotype[cur_haplotype][cur_pos][0]]
      continue
        
    # if there are cryptic mutations
    # check if it's already assigned
    if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
      print('Contains cryptic alleles....')
    assigned_cryptic_alleles = set()
    for alignable_haplotype in alignable_haplotypes:
      alignable_haplotype_pos = position_info_by_haplotype_dict[cur_haplotype][alignable_haplotype][cur_pos-1][1] + 1
      if alignable_haplotype in resolved_alleles_by_haplotype and alignable_haplotype_pos in resolved_alleles_by_haplotype[alignable_haplotype]:
        resolved_allele = resolved_alleles_by_haplotype[alignable_haplotype][alignable_haplotype_pos][0]
        if resolved_allele in cryptic_allele:
          if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
            print(f'Cryptic allele {resolved_allele} already assigned to {alignable_haplotype} at {alignable_haplotype_pos}')
          assigned_cryptic_alleles.add(resolved_allele)
    
    cryptic_alleles_to_assign = set()
    for cryptic_allele in cryptic_allele:
      if cryptic_allele not in assigned_cryptic_alleles:
        cryptic_alleles_to_assign.add(cryptic_allele)

    if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
      print(f'cryptic_alleles_to_assign: {cryptic_alleles_to_assign}\n')

    # if already assigned -> assign to reference allele
    # else                -> assign to unresolved haplotype with the highest likelihood ratio
    if not cryptic_alleles_to_assign:
      if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
        print('no cryptic alleles to assign, assigning reference allele')
      ref_allele = position_info_by_haplotype_dict[cur_haplotype][cur_haplotype][cur_pos-1][0]
      resolved_alleles_by_haplotype[cur_haplotype][cur_pos] = [ref_allele, unresolved_sites_by_haplotype[cur_haplotype][cur_pos][0]]
      continue
    else:
      for cryptic_allele in cryptic_alleles_to_assign:
        max_likelihood_ratio = -1
        max_likelihood_haplotypes = set()
        for alignable_haplotype in alignable_haplotypes:
          if alignable_haplotype not in unresolved_sites_by_haplotype: continue
          alignable_haplotype_pos = position_info_by_haplotype_dict[cur_haplotype][alignable_haplotype][cur_pos-1][1] + 1
          if alignable_haplotype_pos not in unresolved_sites_by_haplotype[alignable_haplotype]: continue
          if not unresolved_sites_by_haplotype[alignable_haplotype][alignable_haplotype_pos]: continue
          likelihood_ratio_by_allele = dict(unresolved_sites_by_haplotype[alignable_haplotype][alignable_haplotype_pos][0])
          if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
            print(alignable_haplotype, list(dict(likelihood_ratio_by_allele).items()))
          if cryptic_allele in likelihood_ratio_by_allele:
            if round(likelihood_ratio_by_allele[cryptic_allele], 10) > max_likelihood_ratio:
              max_likelihood_ratio = round(likelihood_ratio_by_allele[cryptic_allele], 10)
              max_likelihood_haplotypes = {alignable_haplotype}
            elif round(likelihood_ratio_by_allele[cryptic_allele], 10) == max_likelihood_ratio:
              max_likelihood_haplotypes.add(alignable_haplotype)
        if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
          print(f'\nmax_likelihood_haplotypes: {max_likelihood_haplotypes}\nmax_likelihood_ratio: {max_likelihood_ratio}\n')
        if max_likelihood_haplotypes and cur_haplotype in max_likelihood_haplotypes:
          resolved_alleles_by_haplotype[cur_haplotype][cur_pos] = [cryptic_allele, unresolved_sites_by_haplotype[cur_haplotype][cur_pos][0]]
          alleles, fields = unresolved_site_info_by_haplotype[cur_haplotype][cur_pos]
          if cur_haplotype == debug_haplotype and cur_pos in debug_positions:
            print(f'assigning {cryptic_allele} to {cur_haplotype} at {cur_pos}\n')
            print('\t'.join(fields[:-1]), end='')
            print(f'\t{alleles.index(cryptic_allele)}:' + ':'.join(fields[-1].split(':')[1:]))
          out_str_by_haplotype[cur_haplotype][cur_pos] = '\t'.join(fields[:-1]) + f'\t{alleles.index(cryptic_allele)}:' + ':'.join(fields[-1].split(':')[1:])
          break
    




if not debug:
  for haplotype, out_str in out_str_by_haplotype.items():
    output_path = files_by_haplotype[haplotype][2]
    with open(output_path, 'w') as fh:
      fh.write(out_str['headers'] + '\n')
      out_str.pop('headers')
      for cur_pos in sorted(out_str.keys()):
        fh.write(out_str[cur_pos] + '\n')