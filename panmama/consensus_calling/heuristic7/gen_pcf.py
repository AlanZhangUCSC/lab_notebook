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
parser.add_argument('--vcf', type=str, help='Path to the allvar VCF file')
parser.add_argument('--abundances', type=str, help='Path to the haplotype abundance file')
parser.add_argument('--bam', type=str, help='Path to the BAM file')
parser.add_argument('--assignments', type=str, help='Path to reference assignments')
parser.add_argument('--fasta', type=str, help='Path to the aligned FASTA file')
parser.add_argument('--min_coverage', type=int, default=10, help='Minimum coverage to call a variant')
parser.add_argument('--min_allele_depth', type=int, default=5, help='Minimum allele depth to call a variant')
parser.add_argument('--min_coverage_per_ref_assignment_distribution', type=int, default=5, help='Minimum allele count to call a variant per reference assignment distribution')
parser.add_argument('--min_allele_depth_per_ref_assignment_distribution', type=int, default=2, help='Minimum allele depth to call a variant per reference assignment distribution')
parser.add_argument('--error_rate', type=float, default=0.05, help='Sequencing error rate')
args = parser.parse_args()

allvar_vcf_path = args.vcf
abundance_file_path = args.abundances
bam_file_path = args.bam
aligned_fasta_path = args.fasta
reference_assignments_path = args.assignments
haplotype_abundance = {}
reference_assignments = {}
error_rate = args.error_rate

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
aligned_sequences = {}
position_info_by_haplotype = None
sorted_haplotypes = sorted(haplotype_abundance.items(), key=lambda x: x[1])
for i, (haplotype, _) in enumerate(sorted_haplotypes):
  haplotype_to_binary[haplotype] = 2**i
  binary_to_haplotype[2**i] = haplotype


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

for vcf_line in utils.open_file(allvar_vcf_path):
  vcf_line = vcf_line.strip()
  if vcf_line.startswith('#'):
    continue
  fields = vcf_line.split('\t')
  haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, indel = utils.parse_vcf_line(vcf_line)

  # skip if no alternative alleles
  if alt_alleles is None: continue
  alleles = [ref_allele] + alt_alleles
  # skip if overall depth is too low
  if ref_depth + sum(alt_depths) < args.min_coverage: continue

  # skip if no alternative allele has abundance > estimated abundance of the target haplotype and 
  # low abundant alternative alleles are significantly lower than the target haplotype abundance
  possible_alt, possible_alleles = utils.assess_possible_alt(ref_depth, alt_depths, ref_allele, alt_alleles, haplotype_abundance[haplotype], error_rate, args.min_allele_depth)
  if not possible_alt: continue

  query_names = defaultdict(list)
  sam_lines_at_pos = utils.get_sam_lines_at_pos(bam_file_path, haplotype, pos)
  for sam_line in sam_lines_at_pos:
    qname, template_info, base, mapq, pass_filter = utils.process_sam_line(sam_line, pos, indel, alleles)
    if base not in alleles or qname is None: continue
    query_names[qname].append((base, mapq, template_info, pass_filter))

  # Group reads by the numebr of haplotypes they assign to and how different they are from the target haplotype
  ref_assignment_distributions = utils.group_reads_by_haplotype_assignment(query_names, haplotype_to_binary, haplotype, haplotype_abundance, reference_assignments)

  # if pos in [11565, 20562, 22658, 799, 2415, 28951, 27723]:
  #   print(vcf_line)
  #   print('ref_assignment_distributions')
  #   for max_score_diff in ref_assignment_distributions:
  #     print(f'max_score_diff: {max_score_diff}')
  #     for haplotypes_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff].items():
  #       print(f"{haplotypes_binary:0{len(haplotype_to_binary)}b}", 
  #         dict(ref_assignment_distribution.alleles_count),
  #         dict(ref_assignment_distribution.unpassed_alleles_count),
  #         dict(ref_assignment_distribution.all_alleles_count), sep='\t', end='')
  #       print()
  #     print()

  target_candidate_alleles = utils.select_alleles(haplotype, haplotype_abundance, ref_assignment_distributions, haplotype_to_binary, error_rate, possible_alleles, int(args.min_coverage_per_ref_assignment_distribution), int(args.min_allele_depth_per_ref_assignment_distribution), indel)


  if len(target_candidate_alleles) > 1 or (len(target_candidate_alleles) == 1 and list(target_candidate_alleles)[0] != ref_allele):
    utils.print_pcf_line(haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, target_candidate_alleles, ref_assignment_distributions)
  
  