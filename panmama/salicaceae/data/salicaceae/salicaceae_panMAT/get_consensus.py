from Bio import AlignIO
from Bio.Align import AlignInfo
import sys

import argparse

parser = argparse.ArgumentParser(description="Generate a gap consensus sequence from a FASTA alignment.")
parser.add_argument("input_file", help="Input FASTA alignment file")
parser.add_argument("output_file", help="Output file for consensus sequence")
parser.add_argument("--header", default="gap_consensus", help="Sequence header for consensus (default: gap_consensus)")
parser.add_argument("--threshold", type=float, default=0.5, help="Threshold for consensus (default: 0.5)")
parser.add_argument("--ambiguous", default="N", help="Ambiguous characters to include in consensus (default: N)")

args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file
header = args.header
threshold = args.threshold

aln = AlignIO.read(input_file, "fasta")
summary = AlignInfo.SummaryInfo(aln)
gap_consensus = str(summary.gap_consensus(threshold=threshold, ambiguous=ambiguous))
with open(output_file, "w") as f:
    f.write(f">{header}\n")
    f.write(gap_consensus)