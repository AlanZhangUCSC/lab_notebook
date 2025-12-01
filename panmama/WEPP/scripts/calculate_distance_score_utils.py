import bte
import os
import sys
import subprocess

def place_and_get_diff(og_tree_path, ref_path, sample_path, tmp_dir, prefix):
  if not os.path.exists(tmp_dir): os.makedirs(tmp_dir, exist_ok=False)

  sample_fasta = os.path.join(tmp_dir, prefix + '.sample.fasta')
  subprocess.run(["cp", sample_path, sample_fasta], check=True)

  head_cmd = ["head", "-n", "1", sample_fasta]
  header = subprocess.check_output(head_cmd, text=True).strip()
  seq = subprocess.check_output(["tail", "-n", "+2", sample_fasta], text=True).strip()
  tmp_sample_id = header.split(sep=' ')[0][1:] + '_tmp_for_distance'
  with open(sample_fasta, 'w') as f: f.write(f'>{tmp_sample_id}\n{seq}\n')


  tmp_alignment = os.path.join(tmp_dir, prefix + '.aligned.tmp.fasta')
  with open(tmp_alignment, 'w') as f:
      subprocess.run(["mafft", "--quiet", "--auto", "--keeplength", "--addfragments", sample_fasta, ref_path], stdout=f, check=True)

  tmp_new_tree = os.path.join(tmp_dir, prefix + '.new.tmp.pb')
  tmp_vcf = os.path.join(tmp_dir, prefix + '.aligned.tmp.vcf')
  subprocess.run(["conda", "run", "-n", "usher", "faToVcf", tmp_alignment, tmp_vcf], check=True)
  subprocess.run(["conda", "run", "-n", "usher", "usher", "-i", og_tree_path, "-v", tmp_vcf, "-o" , tmp_new_tree, "-d", tmp_dir], check=True)

  tree = bte.MATree(tmp_new_tree)
  haplotype = tree.get_haplotype("0107.CN.2020.BD203AJ.OL401524")
  haplotype_2 = tree.get_haplotype("0107.CN.2020.BD203AJ.OL401524_tmp_for_distance")
  print(f'haplotype num mutations: {len(haplotype)}')
  print(f'haplotype_2 num mutations: {len(haplotype_2)}')
  
  mutations_dict = {}
  mutation_path_file = os.path.join(tmp_dir, 'mutation-paths.txt')
  with open(mutation_path_file, 'r') as f:
    line = f.readline().strip()
    mutation_path = line.split(sep='\t')[1].strip()
    for node_mutations in mutation_path.split(sep=' '):
      mutations = node_mutations.split(sep=':')[1].strip()
      for mutation in mutations.split(sep=','):
        ref = mutation[0]
        alt = mutation[-1]
        coord = int(mutation[1:-1])
        if coord in mutations_dict:
          mutations_dict[coord][1] = alt
        else:
          mutations_dict[coord] = [ref, alt]
  
  num_mutations_to_ref = 0
  for coord in mutations_dict:
    ref, alt = mutations_dict[coord]
    if ref != alt:
      num_mutations_to_ref += 1
  print(f'Number of mutations from sample to reference: {num_mutations_to_ref}')

def get_fasta_path(node, fasta_dir):
  fasta_file = ''.join(c if c.isalnum() else '_' for c in node) + '.unaligned.fasta'
  fasta_path = os.path.join(fasta_dir, fasta_file)
  return fasta_path

def get_sequence_from_fasta(node, fasta_dir):
  fasta_path = get_fasta_path(node, fasta_dir)
  seqs = []
  for idn, seq_str in read_record(fasta_path):
    seqs.append((idn, seq_str))
  return seqs[0][1]

def get_distance(seq1, seq2, get_match_string=True):
  canonical = {
    'a': True, 't': True, 'c': True, 'g': True,
    'A': True, 'T': True, 'C': True, 'G': True
  }

  match_string = ''

  num_errors = 0
  num_snps = 0
  num_gaps = 0
  num_gaps_from_ends = 0
  num_base_to_ambiguous = 0
  for i in range(len(seq1)):
    if seq1[i] != seq2[i]:
      if get_match_string: match_string += ':'
      if seq1[i] == '-' or seq2[i] == '-':
        if seq1[i] in canonical or seq2[i] in canonical:
          num_errors += 1
          num_gaps += 1
      else:
        if seq1[i] in canonical and seq2[i] in canonical:
          num_snps += 1
          num_errors += 1
        else:
          num_base_to_ambiguous += 1
    else:
      if get_match_string: match_string += '|'

  for i in range(len(seq1)):
    if seq1[i] == '-' or seq2[i] == '-':
      if seq1[i] in canonical or seq2[i] in canonical:
        num_gaps_from_ends += 1
    else:
      break
  for i in range(len(seq1)-1, -1, -1):
    if seq1[i] == '-' or seq2[i] == '-':
      if seq1[i] in canonical or seq2[i] in canonical:
        num_gaps_from_ends += 1
    else:
      break
  
  return num_errors, num_snps, num_gaps, num_gaps_from_ends, num_base_to_ambiguous, match_string


        


  
