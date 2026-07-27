import os
import sys
from utils import get_covered_ref_fragments, merge_ref_fragments, align_and_measure_divergence, align_and_measure_divergence_whole_genome

import argparse

parser = argparse.ArgumentParser(
    description="Calculate covered reference distance from FASTQ and abundance/haplotype files"
)
parser.add_argument("abundance_file", help="Path to abundance file")
parser.add_argument("fragment_file", help="Path to fragment file")
parser.add_argument("true_haplotype", help="Name of the true haplotype")
parser.add_argument("--fasta_dir", required=False, default="", help="Path to directory containing FASTA files")
parser.add_argument("-o", "--output_file", required=False, default="", help="Path to output file")
args = parser.parse_args()

abundance_file = args.abundance_file
fragment_file = args.fragment_file
true_haplotype = args.true_haplotype
fasta_dir = args.fasta_dir

fh = None
if args.output_file != "":
  fh = open(args.output_file, "w")

with open(abundance_file, "r") as f:
  for line in f:
    nodes_string, abundance = line.strip().split("\t")
    nodes = nodes_string.split(",")
    if true_haplotype in nodes:
      if fh is not None:
        fh.write(line.strip() + "\t0\t0\n")
      else:
        print(line.strip() + "\t0\t0")
      continue

    rep_node = None
    for node in nodes:
      if not node.startswith("node_"):
        rep_node = node
        break
    if rep_node is None: rep_node = nodes[0]

    divergence = align_and_measure_divergence(os.path.join(fasta_dir, rep_node + ".fasta"), fragment_file)
    divergence_whole_genome = align_and_measure_divergence_whole_genome(os.path.join(fasta_dir, true_haplotype + ".fasta"), os.path.join(fasta_dir, rep_node + ".fasta"))
    total_mismatches = divergence.substitutions + divergence.inserted_bases + divergence.deleted_bases
    total_mismatches_whole_genome = divergence_whole_genome.substitutions + divergence_whole_genome.inserted_bases + divergence_whole_genome.deleted_bases
    if fh is not None:
      fh.write(line.strip() + "\t" + str(total_mismatches) + "\t" + str(total_mismatches_whole_genome) + "\n")
    else:
      print(line.strip() + "\t" + str(total_mismatches) + "\t" + str(total_mismatches_whole_genome))
   


    