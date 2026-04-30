#!/usr/bin/env python3
import sys
from Bio import SeqIO
from collections import Counter

def majority_consensus(aln_file, threshold=0.5, gap_threshold=0.5):
  seqs = [str(rec.seq).upper() for rec in SeqIO.parse(aln_file, "fasta")]
  if not seqs:
    return ""
  
  L = len(seqs[0])
  N = len(seqs)
  consensus = []
  
  for i in range(L):
    column = [s[i] for s in seqs]
    counts = Counter(column)
    
    gap_frac = counts.get('-', 0) / N
    if gap_frac > gap_threshold:
      continue
    
    non_gap = {b: c for b, c in counts.items() if b != '-' and b != 'N'}
    if not non_gap:
      consensus.append('N')
      continue
    
    top_base, top_count = max(non_gap.items(), key=lambda x: x[1])
    non_gap_total = sum(non_gap.values())
    
    if top_count / non_gap_total >= threshold:
      consensus.append(top_base)
    else:
      consensus.append('N')
  
  return "".join(consensus)

if __name__ == "__main__":
  aln = sys.argv[1]
  name = sys.argv[2] if len(sys.argv) > 2 else "consensus"
  cons = majority_consensus(aln)
  print(f">{name}")
  for i in range(0, len(cons), 60):
    print(cons[i:i+60])