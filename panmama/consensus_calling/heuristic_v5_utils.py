import os
import re
import sys
import gzip
import math
import subprocess
import numpy as np
from scipy import stats
from scipy.stats import chi2
from collections import defaultdict
from itertools import product
import warnings

def open_file(path):
  return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def process_aligned_fasta(input_file, reference):
  # Read input FASTA file
  sequences = {}
  with open(input_file, 'r') as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        current_id = line[1:]
        sequences[current_id] = ''
      else:
        sequences[current_id] += line

  ref_seq = None
  ref_id = None
  for seq_id, sequence in sequences.items():
    if reference in seq_id:
      ref_seq = sequence
      ref_id = seq_id
      break

  results = {seq_id: [] for seq_id in sequences}
  num_gaps = {seq_id: 0 for seq_id in sequences}

  ungapped_ref_seqs = defaultdict(str)
  for gapped_pos, ref_base in enumerate(ref_seq):
    if ref_base != '-':
      for seq_id, sequence in sequences.items():
        cur_base = sequence[gapped_pos]
        if cur_base == '-':
          results[seq_id].append((cur_base.upper(), -1))
        else:
          results[seq_id].append((cur_base.upper(), gapped_pos-num_gaps[seq_id]))
          ungapped_ref_seqs[seq_id] += cur_base.upper()
    else:
      for seq_id, sequence in sequences.items():
        if sequence[gapped_pos] != '-':
          ungapped_ref_seqs[seq_id] += sequence[gapped_pos].upper()

    for seq_id, sequence in sequences.items():
      if sequence[gapped_pos] == '-':
        num_gaps[seq_id] += 1

  return results, ungapped_ref_seqs

def parse_vcf_line(vcf_line):
  vcf_line = vcf_line.strip()
  if not vcf_line: return None
  assert not vcf_line.startswith('#')
  fields = vcf_line.split('\t')
  haplotype, pos, ref_allele, alt, sample_info = fields[0], int(fields[1]), fields[3], fields[4], fields[-1]
  indel = False
  if fields[7].startswith('INDEL'):
    indel = True
  alt_alleles = alt.split(',')
  allele_depths = list(map(int, sample_info.split(':')[-1].split(',')))
  ref_depth = allele_depths[0]
  alt_depths = allele_depths[1:]

  if alt == '.':
    return (haplotype, pos, ref_allele, None, ref_depth, None, indel)
  return (haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, indel)

def g_test(observed_allele_proportion, expected_allele_proportion, depth):
  observed_props = np.array([observed_allele_proportion, 1-observed_allele_proportion])
  expected_props = np.array([expected_allele_proportion, 1-expected_allele_proportion])
  observed_counts = observed_props * depth
  expected_counts = expected_props * depth

  g_components = 2 * observed_counts * np.log(observed_props / expected_props)
  g_stat = np.sum(g_components)
  df = len(observed_props) - 1
  p_value = 1 - stats.chi2.cdf(g_stat, df)
  return g_stat, p_value

def assess_possible_alt(ref_depth, alt_depths, ref_allele, alt_alleles, haplotype_abundance, error_rate, min_allele_depth):
  all_depth = ref_depth + sum(alt_depths)
  possible_alt = False
  possible_alleles = set()
  min_allele_depth, min_allele_ratio = map(float, min_allele_depth.split(','))
  if ref_depth >= max(min_allele_depth, min_allele_ratio * all_depth):
    if ref_depth >= all_depth * max(haplotype_abundance - 0.1, 0):
      possible_alleles.add(ref_allele)
    elif ref_depth > 0:
      g_stat, pvalue = g_test(ref_depth / all_depth, max(haplotype_abundance - 0.1, 0), all_depth)
      if pvalue >= 0.05:
        possible_alleles.add(ref_allele)


  for alt_allele, alt_depth in zip(alt_alleles, alt_depths):
    if alt_depth >= max(min_allele_depth, min_allele_ratio * all_depth):
      if alt_depth >= all_depth * max(haplotype_abundance - 0.1, 0):
        possible_alt = True
        possible_alleles.add(alt_allele)
      elif alt_depth > 0:
        g_stat, pvalue = g_test(alt_depth / all_depth, max(haplotype_abundance - 0.1, 0), all_depth)
        if pvalue >= 0.05:
          possible_alt = True
          possible_alleles.add(alt_allele)
  
  return possible_alt, possible_alleles

def get_sam_lines_at_pos(bam_file_path, haplotype, pos):
  if (not os.path.exists(bam_file_path+'.bai')):
    print(f"Indexing BAM file {bam_file_path}...", file=sys.stderr)
    subprocess.run(f"samtools index {bam_file_path}", shell=True, check=True)
  cmd = f"samtools view {bam_file_path} '{haplotype}:{pos}-{pos}'"
  proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
  reads_at_pos = proc.stdout
  return [line.strip() for line in reads_at_pos.split('\n') if line.strip()]

def parse_cigar(cigar_str):
    """Parse CIGAR string and return list of operations and lengths."""
    import re
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)]

def get_ref_to_read_pos(cigar, ref_pos, start_pos):
    """Convert reference position to read position accounting for CIGAR operations."""
    current_ref_pos = start_pos
    current_read_pos = 1  # 1-based position in read
    left_clip = 0
    
    # First handle any clips at the start
    if cigar and cigar[0][1] in 'HS':
        left_clip = cigar[0][0]
    
    for length, op in cigar:
        if op in 'MX=':  # Operations that consume both reference and read
            if current_ref_pos <= ref_pos <= current_ref_pos + length:
                return (ref_pos - current_ref_pos) + current_read_pos
            current_ref_pos += length
            current_read_pos += length
        elif op == 'I':    # Insertion (consumes read but not reference)
            current_read_pos += length
        elif op in 'DN':    # Deletion (consumes reference but not read)
            current_ref_pos += length
        elif op in 'HS':   # Clips (if not at start, must be at end)
            if current_read_pos == 1:  # Left clip
                current_read_pos += length
    
    return None  # Position not found in read


def parse_md_string(md_string):
  md_string = md_string.strip()[5:]
  matches = list(map(int, re.findall(r'\d+', md_string)))
  return sum(matches)

def process_sam_line(line, target_position, indel, alleles):
  # Split the line into fields
  fields = line.strip().split('\t')
  
  # Extract relevant fields
  qname = fields[0]
  flag = int(fields[1])
  start_pos = int(fields[3])
  mapq = int(fields[4])
  cigar = fields[5]
  seq = fields[9]
  md = fields[-1]
  pass_filter = True

  if cigar == '*':
    return None, None, None, None, None, None

  if not (flag & 0x2):
    pass_filter = False
  
  if flag & 0x100:
    return None, None, None, None, None, None

  
  # Parse CIGAR string
  cigar_ops = parse_cigar(cigar)
  
  # Convert reference position to read position
  read_pos = get_ref_to_read_pos(cigar_ops, target_position, start_pos)
  
  # Skip if position not found in read
  if read_pos is None or read_pos <= 0 or read_pos > len(seq):
      return None, None, None, None, None, None

  num_matches = parse_md_string(md)
  
  # Extract the base at the given position
  base = None
  if not indel:
    base = seq[read_pos - 1]  # Convert to 0-based for Python string indexing
    if base not in alleles:
      return None, None, None, None, None, None
  else:
    ref = alleles[0]
    indels = {}
    for allele in alleles[1:]:
      current_indel = (abs(len(allele) - len(ref)), 'I' if len(allele) > len(ref) else 'D')
      indels[''.join(map(str, current_indel))] = allele
    current_base_num = 0
    base = ref
    for i, (length, op) in enumerate(cigar_ops):
      current_base_num += length
      if current_base_num == read_pos:
        try:
          next_length, next_op = cigar_ops[i+1]
          indel_cigar = ''.join([str(next_length), str(next_op)])
          if indel_cigar in indels:
            base = indels[indel_cigar]
            break
        except:
          pass

  template_info = 0 if flag & 0x40 else 1  # Check first/second in pair flag
  assert not (flag & 0x40 and flag & 0x80), "Read cannot be both first and second in pair"
  return qname, template_info, base, mapq, pass_filter, num_matches

class assigned_reference:
  def __init__(self, binary_code, score):
    self.binary_code = binary_code
    self.score = score

  def __str__(self):
    return f"{self.binary_code},{self.score}"
  
  def __repr__(self):
    return f"{self.binary_code},{self.score}"


class ref_assignment_distribution:
  def __init__(self, haplotype_binary, max_score_diff):
    self.haplotype_binary = haplotype_binary
    self.max_score_diff = max_score_diff

    self.alleles_count = defaultdict(int)
    self.alleles_mapq = defaultdict(list)
    self.alleles_num_matches = defaultdict(lambda: defaultdict(int))

    self.unpassed_alleles_count = defaultdict(int)
    self.unpassed_alleles_mapq = defaultdict(list)
    self.unpassed_alleles_num_matches = defaultdict(lambda: defaultdict(int))

    self.all_alleles_count = defaultdict(int)

    self.candidate_allele = ''
  
  def add_allele(self, base, passed_filter, num_matches=None, mapq=None):
    if passed_filter:
      self.alleles_count[base] += 1
      if mapq is not None:
        self.alleles_mapq[base].append(mapq)
      if num_matches is not None:
        self.alleles_num_matches[base][num_matches] += 1
    else:
      self.unpassed_alleles_count[base] += 1
      if mapq is not None:
        self.unpassed_alleles_mapq[base].append(mapq)
      if num_matches is not None:
        self.unpassed_alleles_num_matches[base][num_matches] += 1
    self.all_alleles_count[base] += 1

  def get_total_allele_count(self):
    total_allele_count = defaultdict(int)
    for allele in self.alleles_count:
      total_allele_count[allele] += self.alleles_count[allele]
    for allele in self.unpassed_alleles_count:
      total_allele_count[allele] += self.unpassed_alleles_count[allele]
    return total_allele_count
  
  def get_passed_total_count(self):
    return sum(self.alleles_count.values())
  
  def get_unpassed_total_count(self):
    return sum(self.unpassed_alleles_count.values())

  def get_all_total_count(self):
    return self.get_passed_total_count() + self.get_unpassed_total_count()
  
  def get_alleles(self):
    return sorted(self.alleles_count.keys(), key=lambda x: self.alleles_count[x], reverse=True)
  
  def get_unpassed_alleles(self):
    return sorted(self.unpassed_alleles_count.keys(), key=lambda x: self.unpassed_alleles_count[x], reverse=True)

  def get_all_alleles(self):
    return sorted(self.all_alleles_count.keys(), key=lambda x: self.all_alleles_count[x], reverse=True)
  

def update_ref_assignment_distribution(haplotype, paired_qname, base, pass_filter, mapq, num_matches, haplotype_abundance, reference_assignments, ref_assignment_distributions, haplotype_to_binary):
  haplotypes_binary = 0
  target_score = -1
  highest_score = -1
  for ref in reference_assignments[paired_qname]:
    if ref.binary_code == haplotype_to_binary[haplotype]:
      target_score = ref.score
    if ref.score > highest_score:
      highest_score = ref.score

  for ref in reference_assignments[paired_qname]:
    if ref.score >= target_score:
      haplotypes_binary += ref.binary_code

  max_score_diff = highest_score - target_score
  if max_score_diff not in ref_assignment_distributions:
    ref_assignment_distributions[max_score_diff] = {}
  if haplotypes_binary not in ref_assignment_distributions[max_score_diff]:
    ref_assignment_distributions[max_score_diff][haplotypes_binary] = ref_assignment_distribution(haplotypes_binary, max_score_diff)
  ref_assignment_distributions[max_score_diff][haplotypes_binary].add_allele(base, pass_filter, num_matches, mapq)

def group_reads_by_haplotype_assignment(query_names, haplotype_to_binary, haplotype, haplotype_abundance, reference_assignments):
  ref_assignment_distributions = {}
  for qname, qinfo in query_names.items():
    if len(qinfo) > 2: continue
    if len(qinfo) == 2:
      if (qinfo[0][0] != qinfo[1][0]):
        continue
    for base, mapq, template_info, pass_filter, num_matches in qinfo:
      paired_qname = qname + '/1' if template_info == 0 else qname + '/2'
      update_ref_assignment_distribution(haplotype, paired_qname, base, pass_filter, mapq, num_matches, haplotype_abundance, reference_assignments, ref_assignment_distributions, haplotype_to_binary)
  
  # Sort ref assigment distributions by how close they are to the target haplotype
  ref_assignment_distributions = dict(sorted(ref_assignment_distributions.items(), key=lambda x: x[0]))
  # Sort read groups by the number of haplotypes they assign to
  for max_score_diff in ref_assignment_distributions:
    ref_assignment_distributions[max_score_diff] = dict(sorted(ref_assignment_distributions[max_score_diff].items(), key=lambda x: bin(x[0]).count('1')))
  return ref_assignment_distributions


def get_total_count_of_allele(ref_assignment_distributions, allele=None):
  total_count = 0
  for max_score_diff in ref_assignment_distributions:
    for haplotype_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff]:
      if allele is None:
        total_count += ref_assignment_distribution.get_total_count()
      else:
        total_count += ref_assignment_distribution.alleles_count[allele] + ref_assignment_distribution.unpassed_alleles_count[allele]
  return total_count

def get_passed_count_of_allele(ref_assignment_distributions, allele=None):
  total_count = 0
  for max_score_diff in ref_assignment_distributions:
    for haplotype_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff]:
      if allele is None:
        total_count += ref_assignment_distribution.get_passed_total_count()
      else:
        total_count += ref_assignment_distribution.alleles_count[allele]
  return total_count

def get_unpassed_count_of_allele(ref_assignment_distributions, allele=None):
  total_count = 0
  for max_score_diff in ref_assignment_distributions:
    for haplotype_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff]:
      if allele is None:
        total_count += ref_assignment_distribution.get_unpassed_total_count()
      else:
        total_count += ref_assignment_distribution.unpassed_alleles_count[allele]
  return total_count


def compare_full_alignment(query, start_pos, cigar_ops, ref_seqs, reference):
  cur_query_pos = 0
  cur_ref_pos = start_pos - 1

  query_ranges = []
  ref_ranges = []
  for length, op in cigar_ops:
    if op in 'HP':
      continue
    elif op in 'MX=':
      query_ranges.append((cur_query_pos, cur_query_pos + length - 1))
      ref_ranges.append((cur_ref_pos, cur_ref_pos + length - 1))
      cur_query_pos += length
      cur_ref_pos += length
    elif op in 'IS':
      cur_query_pos += length
    elif op in 'DN':
      cur_ref_pos += length
  
  num_diffs = 0
  for query_range, ref_range in zip(query_ranges, ref_ranges):
    for query_i, ref_i in zip(range(query_range[0], query_range[1]), range(ref_range[0], ref_range[1])):
      if query[query_i] != ref_seqs[reference][ref_i]:
        num_diffs += 1
  
  return num_diffs

def get_full_alignment_scores(target_position, sam_lines_at_pos, ref_seqs, reference, candidate_alleles, indel):
  alignment_scores_by_allele = defaultdict(lambda: defaultdict(int))
  for line in sam_lines_at_pos:
    fields = line.strip().split('\t')
    
    qname = fields[0]
    flag = int(fields[1])
    start_pos = int(fields[3])
    mapq = int(fields[4])
    cigar = fields[5]
    seq = fields[9]

    if cigar == '*' or not (flag & 0x2) or (flag & 0x100):
      continue
    
    # Parse CIGAR string
    cigar_ops = parse_cigar(cigar)
    
    # Convert reference position to read position
    read_pos = get_ref_to_read_pos(cigar_ops, target_position, start_pos)
    
    # Skip if position not found in read
    if read_pos is None or read_pos <= 0 or read_pos > len(seq):
      continue
    
    # Extract the base at the given position
    base = None
    if not indel:
      base = seq[read_pos - 1]  # Convert to 0-based for Python string indexing
      if base not in candidate_alleles:
        continue
    else:
      ref = candidate_alleles[0]
      indels = {}
      for allele in candidate_alleles[1:]:
        current_indel = (abs(len(allele) - len(ref)), 'I' if len(allele) > len(ref) else 'D')
        indels[''.join(map(str, current_indel))] = allele
      current_base_num = 0
      base = ref
      try:
        for i, (length, op) in enumerate(cigar_ops):
          current_base_num += length
          if current_base_num == read_pos:
              next_length, next_op = cigar_ops[i+1]
              indel_cigar = ''.join([str(next_length), str(next_op)])
              if indel_cigar in indels:
                base = indels[indel_cigar]
                break
      except:
        continue

    assert not (flag & 0x40 and flag & 0x80), "Read cannot be both first and second in pair"
    diffs = compare_full_alignment(seq, start_pos, cigar_ops, ref_seqs, reference)
    alignment_scores_by_allele[base][diffs] += 1
  
  for allele in alignment_scores_by_allele:
    alignment_scores_by_allele[allele] = dict(sorted(alignment_scores_by_allele[allele].items(), key=lambda x: x[0]))
  return alignment_scores_by_allele

def select_alleles(haplotype, haplotype_abundance, ref_assignment_distributions, all_matches_by_allele, haplotype_to_binary, error_rate, possible_alleles, min_coverage_per_ref_assignment_distribution, min_allele_depth_per_ref_assignment_distribution, indel, debug_line=False):
  target_haplotype_abundance = haplotype_abundance[haplotype]
  target_haplotype_binary = haplotype_to_binary[haplotype]
  candidate_alleles = set()
  min_max_score_diff = min(list(ref_assignment_distributions.keys()))
  cur_score_distribution = ref_assignment_distributions[min_max_score_diff]
  if target_haplotype_binary in cur_score_distribution:
    assert bin(target_haplotype_binary).count('1') == 1, "Target haplotype binary is not a single haplotype"
    cur_assignment_distribution = cur_score_distribution[target_haplotype_binary]
    alleles = cur_assignment_distribution.get_all_alleles()
    if (alleles[0] in possible_alleles and cur_assignment_distribution.all_alleles_count[alleles[0]] >= min_allele_depth_per_ref_assignment_distribution):
      if len(alleles) == 1:
        candidate_alleles.add(alleles[0])
        return candidate_alleles
      elif len(alleles) > 1:
        # if (cur_assignment_distribution.alleles_count[alleles[0]] > cur_assignment_distribution.alleles_count[alleles[1]]):
        #   candidate_alleles.add(alleles[0])
        #   return candidate_alleles
        max_score_1 = max(all_matches_by_allele[alleles[0]].keys())
        try:
          max_score_2 = max(all_matches_by_allele[alleles[1]].keys())
        except:
          max_score_2 = 0
        tolerance = 2
        if max_score_1 > max_score_2 + tolerance:
          candidate_alleles.add(alleles[0])
          return candidate_alleles
        elif max_score_2 > max_score_1 + tolerance:
          candidate_alleles.add(alleles[1])
          return candidate_alleles
        else:
          if (cur_assignment_distribution.all_alleles_count[alleles[0]] > cur_assignment_distribution.all_alleles_count[alleles[1]]): 
            candidate_alleles.add(alleles[0])
            return candidate_alleles

  for haplotypes_set_binary, cur_ref_assignment_distribution in ref_assignment_distributions[min_max_score_diff].items():
    sorted_pairs = sorted(haplotype_abundance.items(), key=lambda x: x[1])
    sorted_haplotypes = [h for h, _ in sorted_pairs]
    sorted_abundances = [a for _, a in sorted_pairs]

    sum_haplotype_abundance = 0
    for i in range(len(sorted_haplotypes)):
      if haplotypes_set_binary & (1 << i):
        sum_haplotype_abundance += sorted_abundances[i]
    relative_target_haplotype_abundance = max(target_haplotype_abundance - 0.1, 0) / sum_haplotype_abundance

    num_haplotypes_assigned = bin(haplotypes_set_binary).count('1')
    sorted_alleles = sorted(cur_ref_assignment_distribution.alleles_count.items(), key=lambda x: x[1], reverse=True)
    sorted_alleles = sorted_alleles[:min(num_haplotypes_assigned, len(sorted_alleles))]


    cur_total_allele_count = cur_ref_assignment_distribution.get_passed_total_count()

    # depth too low to consider alleles
    if cur_total_allele_count < min_coverage_per_ref_assignment_distribution: continue
    cur_candidate_alleles = set()
    if debug_line:
      print(f"{haplotypes_set_binary:0{len(haplotype_to_binary)}b}")
    for allele, count in sorted_alleles:
      if debug_line:
        print(allele, count, sep='\t')
      if allele in possible_alleles:
        if count < min_allele_depth_per_ref_assignment_distribution:
          if debug_line:
            print(f'skipping {allele} because count < min_allele_depth_per_ref_assignment_distribution: {count}')
          continue
        if debug_line:
          print(f'count / cur_total_allele_count: {count / cur_total_allele_count}')
          print(f'relative_target_haplotype_abundance * (1 - error_rate): {relative_target_haplotype_abundance * (1 - error_rate)}')
        if count / cur_total_allele_count >= relative_target_haplotype_abundance * (1 - error_rate):
          if debug_line:
            print(f'adding {allele} because count / cur_total_allele_count >= relative_target_haplotype_abundance * (1 - error_rate): {count / cur_total_allele_count} >= {relative_target_haplotype_abundance * (1 - error_rate)}')
          cur_candidate_alleles.add(allele)
        else:
          g_stat, pvalue = g_test(count / cur_total_allele_count, relative_target_haplotype_abundance, cur_total_allele_count)
          if pvalue >= 0.05:
            if debug_line:
              print(f'adding {allele} because pvalue >= 0.05: {pvalue}')
            cur_candidate_alleles.add(allele)
          else:
            if debug_line:
              print(f'skipping {allele} because pvalue < 0.05: {pvalue}')
    
    if not cur_candidate_alleles:
      all_high_depth = True
      for allele, count in sorted_alleles:
        if count < min_allele_depth_per_ref_assignment_distribution:
          all_high_depth = False
          break
      if all_high_depth:
        cur_candidate_alleles = set(allele for allele, _ in sorted_alleles)
    for allele in cur_candidate_alleles:
      candidate_alleles.add(allele)
  
  if not candidate_alleles: candidate_alleles = possible_alleles
  return candidate_alleles

def print_pcf_line(haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, target_candidate_alleles, ref_assignment_distributions):
  print(haplotype, pos, ref_allele, ','.join(alt_alleles), ref_depth, ','.join(map(str, alt_depths)), ','.join(target_candidate_alleles), sep='\t', end='\t')
  ref_assignment_distributions_str = ''
  all_alleles = [ref_allele] + alt_alleles
  for max_score_diff in ref_assignment_distributions:
    ref_assignment_distributions_str += f'{max_score_diff}:'
    for haplotypes_set_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff].items():
      ref_assignment_distributions_str += f'{haplotypes_set_binary},'
      for allele in all_alleles:
        if allele in ref_assignment_distribution.alleles_count:
          ref_assignment_distributions_str += f'{ref_assignment_distribution.alleles_count[allele]},'
        else:
          ref_assignment_distributions_str += '0,'
      for allele in all_alleles:
        if allele in ref_assignment_distribution.unpassed_alleles_count:
          ref_assignment_distributions_str += f'{ref_assignment_distribution.unpassed_alleles_count[allele]},'
        else:
          ref_assignment_distributions_str += '0,'
    ref_assignment_distributions_str = ref_assignment_distributions_str[:-1] + ';'
  
  print(ref_assignment_distributions_str)

def parse_ref_assignment_distribution(ref_assignment_distribution_str, alleles):
  ref_assignment_distributions = {}
  
  for by_max_score_diff in ref_assignment_distribution_str.split(';'):
    if by_max_score_diff == '': continue
    max_score_diff, distribution_str = by_max_score_diff.split(':')
    max_score_diff = float(max_score_diff)
    ref_assignment_distributions[max_score_diff] = {}
    distribution_fields = distribution_str.split(',')
    assert(len(distribution_fields) % (len(alleles) * 2 + 1) == 0)
    for i in range(0, len(distribution_fields), len(alleles) * 2 + 1):
      by_haplotypes_set_binary = distribution_fields[i:i+len(alleles)*2+1]
      haplotypes_set_binary = int(by_haplotypes_set_binary[0])
      ref_assignment_distributions[max_score_diff][haplotypes_set_binary] = ref_assignment_distribution(haplotypes_set_binary, max_score_diff)
      for j, (allele, count) in enumerate(zip(alleles * 2, by_haplotypes_set_binary[1:])):
        count = int(count)
        if j < len(alleles):
          for k in range(count):
            ref_assignment_distributions[max_score_diff][haplotypes_set_binary].add_allele(allele, True, None)
        else:
          for k in range(count):
            ref_assignment_distributions[max_score_diff][haplotypes_set_binary].add_allele(allele, False, None)
    

  return ref_assignment_distributions

def parse_pcf_line(pcf_line):
  fields = pcf_line.split('\t')
  haplotype, pos = fields[0], int(fields[1])
  ref_allele = fields[2]
  alt_alleles = fields[3].split(',')
  ref_depth = int(fields[4])
  alt_depths = list(map(int, fields[5].split(',')))
  candidate_alleles = fields[6].split(',')
  ref_assignment_distribution_str = fields[7]
  ref_assignment_distributions = parse_ref_assignment_distribution(ref_assignment_distribution_str, [ref_allele] + alt_alleles)
  return haplotype, pos, ref_allele, alt_alleles, ref_depth, alt_depths, candidate_alleles, ref_assignment_distributions

def calculate_mle(observed_allele_frequency, expected_allele_frequency):
  smoothed_expected_allele_frequency = {allele: (freq + 1e-10) for allele, freq in expected_allele_frequency.items()}
  return sum(observed_allele_frequency[allele] * math.log(smoothed_expected_allele_frequency[allele]) for allele in observed_allele_frequency)

def cross_reference_alleles(pos, target_haplotype, position_info_by_haplotype, ref_seqs, ref_assignment_distributions, pcf_info, possible_alleles, haplotype_abundance, binary_to_haplotype, haplotype_to_binary, min_coverage_per_ref_assignment_distribution, min_allele_depth_per_ref_assignment_distribution, error_rate, debug_line=False):
  target_haplotype_abundance = haplotype_abundance[target_haplotype]
  target_haplotype_binary = haplotype_to_binary[target_haplotype]
  candidate_alleles = set()
  min_max_score_diff = min(list(ref_assignment_distributions.keys()))
  cur_score_distribution = ref_assignment_distributions[min_max_score_diff]

  cross_referenced_alleles_by_ref_assignment_distribution = defaultdict(lambda: defaultdict(list))
  for haplotypes_set_binary, cur_ref_assignment_distribution in ref_assignment_distributions[min_max_score_diff].items():
    sorted_pairs = sorted(haplotype_abundance.items(), key=lambda x: x[1])
    sorted_haplotypes = [h for h, _ in sorted_pairs]
    sorted_abundances = [a for _, a in sorted_pairs]

    sum_haplotype_abundance = 0
    for i in range(len(sorted_haplotypes)):
      if haplotypes_set_binary & (1 << i):
        sum_haplotype_abundance += sorted_abundances[i]
    relative_target_haplotype_abundance = target_haplotype_abundance / sum_haplotype_abundance

    num_haplotypes_assigned = bin(haplotypes_set_binary).count('1')
    sorted_alleles = sorted(cur_ref_assignment_distribution.alleles_count.items(), key=lambda x: x[1], reverse=True)
    sorted_alleles = sorted_alleles[:min(num_haplotypes_assigned, len(sorted_alleles))]

    cur_total_allele_count = cur_ref_assignment_distribution.get_passed_total_count()

    # depth too low to consider alleles
    if cur_total_allele_count < 1 / relative_target_haplotype_abundance or cur_total_allele_count < min_coverage_per_ref_assignment_distribution: continue
    cur_candidate_alleles = set()

    for allele, count in sorted_alleles:
      if allele in possible_alleles:
        if count < min_allele_depth_per_ref_assignment_distribution: continue
        if count / cur_total_allele_count >= relative_target_haplotype_abundance * (1 - error_rate):
          cur_candidate_alleles.add(allele)
        else:
          g_stat, pvalue = g_test(count / cur_total_allele_count, relative_target_haplotype_abundance, cur_total_allele_count)
          if pvalue >= 0.05:
            cur_candidate_alleles.add(allele)

    if len(cur_candidate_alleles) == 1:
      cross_referenced_alleles_by_ref_assignment_distribution[min_max_score_diff][haplotypes_set_binary].append([list(cur_candidate_alleles)[0], 0, 0])
    else:
      other_haplotypes_info = []
      haplotypes = binary_set_to_haplotypes(haplotypes_set_binary, binary_to_haplotype)
      alleles_list = []
      for cur_haplotype in haplotypes:
        if cur_haplotype == target_haplotype:
          alleles_list.append(list(cur_candidate_alleles))
        else:
          cur_hap_ref_base, cur_hap_pos = position_info_by_haplotype[cur_haplotype][pos-1]
          other_haplotypes_info.append([cur_haplotype, cur_hap_ref_base, cur_hap_pos])
          cur_hap_pos += 1
          if cur_hap_pos in pcf_info[cur_haplotype]:
            alleles_list.append(pcf_info[cur_haplotype][cur_hap_pos].split('\t')[6].split(','))
          else:
            if cur_hap_ref_base in ('A', 'C', 'G', 'T'):
              alleles_list.append([cur_hap_ref_base])
            else:
              alleles_list.append(list(possible_alleles))

      if debug_line:
        print()
        print(f'haplotypes_set_binary: {haplotypes_set_binary:0{len(haplotype_to_binary)}b}')
        for _haplotype, _alleles in zip(haplotypes, alleles_list):
          print(_haplotype, _alleles)
        for _haplotype, _ref_base, _pos in other_haplotypes_info:
          print(_haplotype, _ref_base, _pos)
        print()

      allele_combinations = list(product(*alleles_list))
      allele_combinations_dict = []
      for allele_combination in allele_combinations:
        combination_dict = {}
        for haplotype, allele in zip(haplotypes, allele_combination):
          combination_dict[haplotype] = allele
        allele_combinations_dict.append(combination_dict)
      
      allele_combinations_expected_allele_frequency = []

      observed_allele_frequency = defaultdict(int)
      observed_allele_relative_frequency = defaultdict(int)
      observed_total_allele_count = 0
      for allele in cur_candidate_alleles:
        observed_allele_frequency[allele] = cur_ref_assignment_distribution.alleles_count[allele]
        observed_total_allele_count += cur_ref_assignment_distribution.alleles_count[allele]
      observed_allele_relative_frequency = {allele: observed_allele_frequency[allele] / observed_total_allele_count for allele in cur_candidate_alleles}

      kl_divergences = []
      mles = []
      total_abundance = sum(haplotype_abundance[h] for h in haplotypes)
      for i, allele_combination in enumerate(allele_combinations):
        expected_allele_frequency = defaultdict(int)
        for haplotype, allele in zip(haplotypes, allele_combination):
          expected_allele_frequency[allele] += haplotype_abundance[haplotype] / total_abundance
        allele_combinations_expected_allele_frequency.append(expected_allele_frequency)
        kl_divergences.append(calculate_kl_divergence(observed_allele_relative_frequency, expected_allele_frequency))
        mles.append(calculate_mle(observed_allele_frequency, expected_allele_frequency))
      
      # print(observed_allele_relative_frequency)
      # for kl_divergence, mle, allele_combination, expected_allele_frequency in zip(kl_divergences, mles, allele_combinations_dict, allele_combinations_expected_allele_frequency):
      #   print(kl_divergence, mle, allele_combination, dict(expected_allele_frequency))
      # print()

      max_mle_by_allele = {}
      for mle, kl_divergence, allele_combination in zip(mles, kl_divergences, allele_combinations_dict):
        allele = allele_combination[target_haplotype]
        if allele not in max_mle_by_allele or mle > max_mle_by_allele[allele][0]:
          max_mle_by_allele[allele] = [mle, kl_divergence]
      
      sorted_max_mle_by_allele = dict(sorted(max_mle_by_allele.items(), key=lambda item: item[1], reverse=True))
      
      for allele, (mle, kl_divergence) in sorted_max_mle_by_allele.items():
        cross_referenced_alleles_by_ref_assignment_distribution[min_max_score_diff][haplotypes_set_binary].append([allele, mle, kl_divergence])

  return cross_referenced_alleles_by_ref_assignment_distribution

def cross_reference_alleles_2(pos, target_haplotype, position_info_by_haplotype, ref_seqs, ref_assignment_distributions, pcf_info, possible_alleles, haplotype_abundance, binary_to_haplotype, haplotype_to_binary, min_coverage_per_ref_assignment_distribution, min_allele_depth_per_ref_assignment_distribution, error_rate, debug_line=False):
  target_haplotype_abundance = haplotype_abundance[target_haplotype]
  target_haplotype_binary = haplotype_to_binary[target_haplotype]
  candidate_alleles = set()
  min_max_score_diff = min(list(ref_assignment_distributions.keys()))
  cur_score_distribution = ref_assignment_distributions[min_max_score_diff]

  cross_referenced_alleles_by_ref_assignment_distribution = defaultdict(lambda: defaultdict(list))
  for haplotypes_set_binary, cur_ref_assignment_distribution in ref_assignment_distributions[min_max_score_diff].items():
    sorted_pairs = sorted(haplotype_abundance.items(), key=lambda x: x[1])
    sorted_haplotypes = [h for h, _ in sorted_pairs]
    sorted_abundances = [a for _, a in sorted_pairs]

    sum_haplotype_abundance = 0
    for i in range(len(sorted_haplotypes)):
      if haplotypes_set_binary & (1 << i):
        sum_haplotype_abundance += sorted_abundances[i]
    relative_target_haplotype_abundance = target_haplotype_abundance / sum_haplotype_abundance

    num_haplotypes_assigned = bin(haplotypes_set_binary).count('1')
    if num_haplotypes_assigned == 1: continue
    sorted_alleles = sorted(cur_ref_assignment_distribution.alleles_count.items(), key=lambda x: x[1], reverse=True)
    sorted_alleles = sorted_alleles[:min(num_haplotypes_assigned, len(sorted_alleles))]

    cur_total_allele_count = cur_ref_assignment_distribution.get_passed_total_count()

    # depth too low to consider alleles
    if cur_total_allele_count < 1 / relative_target_haplotype_abundance or cur_total_allele_count < min_coverage_per_ref_assignment_distribution: continue
    cur_candidate_alleles = set()

    for allele, count in sorted_alleles:
      if allele in possible_alleles:
        if count < min_allele_depth_per_ref_assignment_distribution: continue
        cur_candidate_alleles.add(allele)
        # if count / cur_total_allele_count >= relative_target_haplotype_abundance * (1 - error_rate):
        #   cur_candidate_alleles.add(allele)
        # else:
        #   g_stat, pvalue = g_test(count / cur_total_allele_count, relative_target_haplotype_abundance, cur_total_allele_count)
        #   if pvalue >= 0.05:
        #     cur_candidate_alleles.add(allele)

    if len(cur_candidate_alleles) == 1:
      cross_referenced_alleles_by_ref_assignment_distribution[min_max_score_diff][haplotypes_set_binary].append([list(cur_candidate_alleles)[0], 0, 0])
      if debug_line:
        other_haplotypes_info = []
        haplotypes = binary_set_to_haplotypes(haplotypes_set_binary, binary_to_haplotype)
        alleles_list = []
        for cur_haplotype in haplotypes:
          if cur_haplotype == target_haplotype:
            alleles_list.append(list(cur_candidate_alleles))
          else:
            cur_hap_ref_base, cur_hap_pos = position_info_by_haplotype[cur_haplotype][pos-1]
            other_haplotypes_info.append([cur_haplotype, cur_hap_ref_base, cur_hap_pos])
            cur_hap_pos += 1
            if cur_hap_pos in pcf_info[cur_haplotype]:
              alleles_list.append(pcf_info[cur_haplotype][cur_hap_pos].split('\t')[6].split(','))
            else:
              if cur_hap_ref_base in ('A', 'C', 'G', 'T'):
                alleles_list.append([cur_hap_ref_base])
              else:
                alleles_list.append(list(possible_alleles))

        if debug_line:
          print()
          print(f'haplotypes_set_binary: {haplotypes_set_binary:0{len(haplotype_to_binary)}b}')
          for _haplotype, _alleles in zip(haplotypes, alleles_list):
            print(_haplotype, _alleles)
          for _haplotype, _ref_base, _pos in other_haplotypes_info:
            print(_haplotype, _ref_base, _pos)
          print()
    else:
      other_haplotypes_info = []
      haplotypes = binary_set_to_haplotypes(haplotypes_set_binary, binary_to_haplotype)
      alleles_list = []
      for cur_haplotype in haplotypes:
        if cur_haplotype == target_haplotype:
          alleles_list.append(list(cur_candidate_alleles))
        else:
          cur_hap_ref_base, cur_hap_pos = position_info_by_haplotype[cur_haplotype][pos-1]
          other_haplotypes_info.append([cur_haplotype, cur_hap_ref_base, cur_hap_pos])
          cur_hap_pos += 1
          if cur_hap_pos in pcf_info[cur_haplotype]:
            alleles_list.append(pcf_info[cur_haplotype][cur_hap_pos].split('\t')[6].split(','))
          else:
            if cur_hap_ref_base in ('A', 'C', 'G', 'T'):
              alleles_list.append([cur_hap_ref_base])
            else:
              alleles_list.append(list(possible_alleles))

      if debug_line:
        print()
        print(f'haplotypes_set_binary: {haplotypes_set_binary:0{len(haplotype_to_binary)}b}')
        for _haplotype, _alleles in zip(haplotypes, alleles_list):
          print(_haplotype, _alleles)
        for _haplotype, _ref_base, _pos in other_haplotypes_info:
          print(_haplotype, _ref_base, _pos)
        print()

      allele_combinations = list(product(*alleles_list))
      allele_combinations_dict = []
      for allele_combination in allele_combinations:
        combination_dict = {}
        for haplotype, allele in zip(haplotypes, allele_combination):
          combination_dict[haplotype] = allele
        allele_combinations_dict.append(combination_dict)
      
      allele_combinations_expected_allele_frequency = []

      observed_allele_frequency = defaultdict(int)
      observed_allele_relative_frequency = defaultdict(int)
      observed_total_allele_count = 0
      for allele in cur_candidate_alleles:
        observed_allele_frequency[allele] = cur_ref_assignment_distribution.alleles_count[allele]
        observed_total_allele_count += cur_ref_assignment_distribution.alleles_count[allele]
      observed_allele_relative_frequency = {allele: observed_allele_frequency[allele] / observed_total_allele_count for allele in cur_candidate_alleles}

      kl_divergences = []
      mles = []
      total_abundance = sum(haplotype_abundance[h] for h in haplotypes)
      for i, allele_combination in enumerate(allele_combinations):
        expected_allele_frequency = defaultdict(int)
        for haplotype, allele in zip(haplotypes, allele_combination):
          expected_allele_frequency[allele] += haplotype_abundance[haplotype] / total_abundance
        allele_combinations_expected_allele_frequency.append(expected_allele_frequency)
        kl_divergences.append(calculate_kl_divergence(observed_allele_relative_frequency, expected_allele_frequency))
        mles.append(calculate_mle(observed_allele_frequency, expected_allele_frequency))
      
      # print(observed_allele_relative_frequency)
      # for kl_divergence, mle, allele_combination, expected_allele_frequency in zip(kl_divergences, mles, allele_combinations_dict, allele_combinations_expected_allele_frequency):
      #   print(kl_divergence, mle, allele_combination, dict(expected_allele_frequency))
      # print()

      max_mle_by_allele = {}
      for mle, kl_divergence, allele_combination in zip(mles, kl_divergences, allele_combinations_dict):
        allele = allele_combination[target_haplotype]
        if allele not in max_mle_by_allele or mle > max_mle_by_allele[allele][0]:
          max_mle_by_allele[allele] = [mle, kl_divergence]
      
      sorted_max_mle_by_allele = dict(sorted(max_mle_by_allele.items(), key=lambda item: item[1], reverse=True))
      
      for allele, (mle, kl_divergence) in sorted_max_mle_by_allele.items():
        cross_referenced_alleles_by_ref_assignment_distribution[min_max_score_diff][haplotypes_set_binary].append([allele, mle, kl_divergence])

  return cross_referenced_alleles_by_ref_assignment_distribution

def score_alleles(target_position, haplotype, haplotype_abundance, ref_assignment_distributions, haplotype_to_binary, binary_to_haplotype, error_rate, possible_alleles, ref_allele, sam_lines_at_pos, ref_seqs, indel, debug=False):
  # pass
  target_haplotype_abundance = haplotype_abundance[haplotype]
  target_haplotype_binary = haplotype_to_binary[haplotype]
  candidate_allele_by_ref_assignment_distribution = defaultdict(lambda: defaultdict(set))
  final_candidate_alleles = set()
  if 0 in ref_assignment_distributions:
    cur_score_distribution = dict(ref_assignment_distributions[0])
    if target_haplotype_binary in cur_score_distribution:
      assert bin(target_haplotype_binary).count('1') == 1, "Target haplotype binary is not a single haplotype"
      cur_assignment_distribution = cur_score_distribution[target_haplotype_binary]
      alleles = cur_assignment_distribution.get_alleles()
      if len(alleles) == 1 or (len(alleles) > 1 and cur_assignment_distribution.alleles_count[alleles[0]] > cur_assignment_distribution.alleles_count[alleles[1]]):
        if alleles[0] in possible_alleles:
          final_candidate_alleles.add(alleles[0])
          candidate_allele_by_ref_assignment_distribution[0][target_haplotype_binary].add(alleles[0])
          return final_candidate_alleles, candidate_allele_by_ref_assignment_distribution
    
    candidate_alleles = set()
    for haplotypes_set_binary, cur_ref_assignment_distribution in ref_assignment_distributions[0].items():
      sorted_pairs = sorted(haplotype_abundance.items(), key=lambda x: x[1])
      sorted_haplotypes = [h for h, _ in sorted_pairs]
      sorted_abundances = [a for _, a in sorted_pairs]

      sum_haplotype_abundance = 0
      for i in range(len(sorted_haplotypes)):
        if haplotypes_set_binary & (1 << i):
          sum_haplotype_abundance += sorted_abundances[i]
      relative_target_haplotype_abundance = target_haplotype_abundance / sum_haplotype_abundance

      num_haplotypes_assigned = bin(haplotypes_set_binary).count('1')
      sorted_alleles = sorted(cur_ref_assignment_distribution.alleles_count.items(), key=lambda x: x[1], reverse=True)
      sorted_alleles = sorted_alleles[:min(num_haplotypes_assigned, len(sorted_alleles))]


      cur_total_allele_count = cur_ref_assignment_distribution.get_passed_total_count()

      # depth too low to consider alleles
      if cur_total_allele_count < 1 / relative_target_haplotype_abundance: continue

      # if len(sorted_alleles) == 1:
      #   cur_candidate_allele = sorted_alleles[0][0]
      #   if cur_candidate_allele in possible_alleles:
      #     candidate_allele_by_ref_assignment_distribution[0][haplotypes_set_binary].add(cur_candidate_allele)
      # else:
      #   if num_haplotypes_assigned <= 3:
      #     haplotypes = binary_set_to_haplotypes(haplotypes_set_binary, binary_to_haplotype)


      for allele, count in sorted_alleles:
        if allele in possible_alleles and count / cur_total_allele_count >= relative_target_haplotype_abundance * (1-error_rate):
          candidate_allele_by_ref_assignment_distribution[0][haplotypes_set_binary].add(allele)
      

    if len(candidate_allele_by_ref_assignment_distribution[0]) == 0:
      candidate_alleles = possible_alleles
    elif len(candidate_allele_by_ref_assignment_distribution[0]) == 1:
      haplotypes_set_binary, alleles = list(candidate_allele_by_ref_assignment_distribution[0].items())[0]
      if len(alleles) == 1:
        candidate_alleles = alleles
      else:
        candidate_alleles = possible_alleles
    else:
      for haplotypes_set_binary, alleles in candidate_allele_by_ref_assignment_distribution[0].items():
        if len(alleles) == 1:
          candidate_alleles.add(list(alleles)[0])
          
      # for haplotypes_set_binary, alleles in candidate_allele_by_ref_assignment_distribution[0].items():
      #   if len(alleles) == 1:
      #     cur_candidate_allele = list(alleles)[0]
      #     negate_allele = False
      #     for max_score_diff in ref_assignment_distributions.keys():
      #       if max_score_diff == 0: continue
      #       if haplotypes_set_binary in ref_assignment_distributions[max_score_diff]:
      #         cur_score_ref_assignment_distribution = ref_assignment_distributions[max_score_diff][haplotypes_set_binary]
      #         if len(cur_score_ref_assignment_distribution.alleles_count) == 1:
      #           cur_score_allele = list(cur_score_ref_assignment_distribution.alleles_count.keys())[0]
      #           candidate_allele_by_ref_assignment_distribution[max_score_diff][haplotypes_set_binary].add(cur_score_allele)
      #           if cur_score_allele == cur_candidate_allele:
      #             negate_allele = True
      #             break
      #     if not negate_allele:
      #       candidate_alleles.add(cur_candidate_allele)

    if not candidate_alleles:
      candidate_alleles = possible_alleles
    
    if debug and len(candidate_alleles) > 1:
      alignment_scores_by_allele = get_full_alignment_scores(target_position, sam_lines_at_pos, ref_seqs, haplotype, candidate_alleles, indel)
      for allele in alignment_scores_by_allele:
        for diffs, count in alignment_scores_by_allele[allele].items():
            print(f'{allele}: {diffs} {count}')
        print()
    return candidate_alleles, candidate_allele_by_ref_assignment_distribution

  
  return possible_alleles, candidate_allele_by_ref_assignment_distribution

def calculate_kl_divergence(observed_allele_frequency, expected_allele_frequency):
  kl_divergence = 0
  for allele, observed_freq in observed_allele_frequency.items():
    kl_divergence += observed_freq * math.log(observed_freq / (expected_allele_frequency[allele] + 1e-10))
  return kl_divergence

def calculate_allele_combination_kl_divergence(candidate_alleles_by_haplotype, haplotype_abundance, observed_allele_frequency):  
  haplotypes = list(candidate_alleles_by_haplotype.keys())
  allele_lists = []
  for haplotype, alleles in candidate_alleles_by_haplotype.items():
    allele_lists.append(alleles)
  allele_combinations = list(product(*allele_lists))

  allele_combinations_dict = []
  for allele_combination in allele_combinations:
    combination_dict = {}
    for haplotype, allele in zip(haplotypes, allele_combination):
      combination_dict[haplotype] = allele
    allele_combinations_dict.append(combination_dict)
  
  allele_combinations_expected_allele_frequency = []

  kl_divergences = []
  total_abundance = sum(haplotype_abundance[h] for h in candidate_alleles_by_haplotype)
  for i, allele_combination in enumerate(allele_combinations):
    expected_allele_frequency = defaultdict(int)
    for haplotype, allele in zip(haplotypes, allele_combination):
      expected_allele_frequency[allele] += haplotype_abundance[haplotype] / total_abundance
    allele_combinations_expected_allele_frequency.append(expected_allele_frequency)
    kl_divergences.append(calculate_kl_divergence(observed_allele_frequency, expected_allele_frequency))



  return allele_combinations_dict, kl_divergences, allele_combinations_expected_allele_frequency

def get_alignable_haplotypes(ref_assignment_distributions, binary_to_haplotype):
  alignable_haplotypes = set()
  for max_score_diff in ref_assignment_distributions:
    for haplotypes_binary, ref_assignment_distribution in ref_assignment_distributions[max_score_diff].items():
      # Split into individual binary numbers with single 1
      binary_num = int(haplotypes_binary)
      while binary_num:
        # Extract rightmost 1 bit
        single_one = binary_num & (-binary_num)
        if single_one in binary_to_haplotype:
          alignable_haplotypes.add(binary_to_haplotype[single_one])
        binary_num &= (binary_num - 1)  # Clear rightmost 1 bit
  return alignable_haplotypes

def binary_set_to_haplotypes(haplotypes_set_binary, binary_to_haplotype):
  haplotypes = set()
  while haplotypes_set_binary:
    single_one = haplotypes_set_binary & (-haplotypes_set_binary)
    if single_one in binary_to_haplotype:
      haplotypes.add(binary_to_haplotype[single_one])
    haplotypes_set_binary &= (haplotypes_set_binary - 1)
  return haplotypes

