from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re

def split_fasta(input_file, n_length=100):
    records = list(SeqIO.parse(input_file, "fasta"))
    
    first_seq = str(records[0].seq)
    pattern = f"N{{{n_length},}}"
    
    splits = []
    for m in re.finditer(pattern, first_seq, re.IGNORECASE):
        splits.append((m.start(), m.end()))
    
    if not splits:
        print("No separator regions found.")
        sys.exit(1)
    
    boundaries = []
    prev = 0
    for start, end in splits:
        boundaries.append((prev, start))
        prev = end
    boundaries.append((prev, len(first_seq)))
    
    print(f"Found {len(boundaries)} regions at positions: {boundaries}")
    
    # write each region
    for i, (reg_start, reg_end) in enumerate(boundaries):
        out_records = []
        for rec in records:
            out_records.append(rec[reg_start:reg_end])
        
        out_file = f"region_{i+1}.fasta"
        SeqIO.write(out_records, out_file, "fasta")
        print(f"Wrote {len(out_records)} sequences to {out_file}")
    
    # write concatenated output with Ns stripped
    stripped_records = []
    for rec in records:
        joined_seq = "".join(str(rec.seq[s:e]) for s, e in boundaries)
        stripped_records.append(SeqRecord(seq=type(rec.seq)(joined_seq), id=rec.id, description=rec.description))
    
    SeqIO.write(stripped_records, "all_regions_stripped.fasta", "fasta")
    print(f"Wrote {len(stripped_records)} stripped sequences to all_regions_stripped.fasta")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python split_fasta.py <input.fasta> [n_length]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    n_length = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    split_fasta(input_file, n_length)