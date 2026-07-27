import os
import sys
from utils import get_covered_ref_fragments, merge_ref_fragments, align_and_measure_divergence

import argparse

parser = argparse.ArgumentParser(
    description="Calculate covered reference distance from FASTQ and abundance/haplotype files"
)
parser.add_argument("abundance_file", help="Path to abundance file")
parser.add_argument("fragment_dir", help="Path to fragment directory")
parser.add_argument("manifest_file", help="Path to manifest file")
parser.add_argument("--total_fragments", required=True, type=int, help="Total number of fragments")
parser.add_argument("--fasta_dir", required=False, default="", help="Path to directory containing FASTA files")
args = parser.parse_args()

abundance_file = args.abundance_file
fragment_dir = args.fragment_dir
manifest_file = args.manifest_file
fasta_dir = args.fasta_dir
total_fragments = args.total_fragments

true_haplotypes = {}
with open(manifest_file, "r") as f:
  for line in f:
    line = line.strip()
    if line.startswith('n_total'): continue
    fields = line.split("\t")
    true_haplotype = fields[1]
    proportion = float(fields[2])
    fragment_file = os.path.join(fragment_dir, f'{true_haplotype}_{total_fragments}.og_fragment.fa.gz')
    true_haplotypes[true_haplotype] = (fragment_file, proportion)
    

with open(abundance_file, "r") as f:
  for line in f:
    nodes_string, abundance = line.strip().split("\t")
    nodes = nodes_string.split(",")
    # if true_haplotype in nodes:
    #   print(line.strip() + "\t0")
    #   continue

    rep_node = None
    for node in nodes:
      if not node.startswith("node_"):
        rep_node = node
        break
    if rep_node is None: rep_node = nodes[0]
    
    min_divergence = float('inf')
    min_divergence_haplotype = None
    for true_haplotype, info in true_haplotypes.items():
      fragment_file, proportion = info
      divergence = align_and_measure_divergence(os.path.join(fasta_dir, rep_node + ".fasta"), fragment_file)
      total_mismatches = divergence.substitutions + divergence.inserted_bases + divergence.deleted_bases
      if total_mismatches < min_divergence:
        min_divergence = total_mismatches
        min_divergence_haplotype = true_haplotype


    print(line.strip() + "\t" + min_divergence_haplotype + "\t" + str(min_divergence))
    
   


    