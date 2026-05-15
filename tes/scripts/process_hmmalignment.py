import argparse
import Bio.SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="Input hmmalignment file.")
parser.add_argument("output", type=str, help="Output file")
parser.add_argument("-p", "--polya_file", type=str, required=True, help="Polya file")
parser.add_argument("-a", "--alignment_file", type=str, required=False, help="Alignment file")
parser.add_argument("-k", "--keep-non-envelope", action="store_true", help="Keep non-envelope regions")
args = parser.parse_args()

records = list(Bio.SeqIO.parse(args.input, "fasta"))

input_family = records[0].id.split('|')[1]
input_family = input_family.replace('#SINE/Alu', '').strip()

non_polya_hmm_len = -1
with open(args.polya_file, "r") as f:
  for line in f:
    family, accession, nonpolya_len = line.strip().split('\t')
    if family == input_family:
      non_polya_hmm_len = int(nonpolya_len)
      break
if non_polya_hmm_len == -1:
  raise ValueError(f"Family {input_family} not found in {args.polya_file}")
  exit(1)


hmm_position_count = 0
non_polya_global_len = 0
for i, char in enumerate(records[0].seq):
  if char == '.':
    continue
  elif char == '-':
    hmm_position_count += 1
  elif char.isupper():
    hmm_position_count += 1
  elif char.islower():
    continue
  else:
    raise ValueError(f"Invalid character: {char}")

  if hmm_position_count == non_polya_hmm_len:
    non_polya_global_len = i + 1
    break

# 0 based, last position of the non-polyA sequence
non_polya_global_terminal_pos = non_polya_global_len - 1

alignment_records = []
if args.alignment_file:
  with open(args.alignment_file, "r") as f:
    for line in f:
      fields = line.strip().split('\t')
      seq_id, hmmfrom, hmmto, alifrom, alito, envfrom, envto = fields[2], int(fields[4]), int(fields[5]), int(fields[6]), int(fields[7]), int(fields[8]), int(fields[9])
      alignment_records.append((seq_id, hmmfrom, hmmto, alifrom, alito, envfrom, envto))
if alignment_records:
  if len(alignment_records) != len(records):
    raise ValueError(f"Number of alignment records ({len(alignment_records)}) does not match number of records ({len(records)})")
    exit(1)

out_fh = open(args.output, "w")
for i, record in enumerate(records):
  keep_seq = None
  if args.keep_non_envelope or not alignment_records:
    keep_seq = record.seq[:non_polya_global_terminal_pos+1].replace('-', '').replace('.', '').upper()
  else:
    alignment_record = alignment_records[i]
    assert record.id == alignment_record[0]
    envfrom, envto = alignment_record[5], alignment_record[6] # 1 based, inclusive
    envfrom_global, envto_global = -1, -1 # 0 based, inclusive
    if envfrom > envto: continue
    base_count = 0
    for j, char in enumerate(record.seq):
      if   char == '.': continue
      elif char == '-': continue
      else:
        base_count += 1
        if base_count == envfrom:
          envfrom_global = j
        if base_count == envto:
          envto_global = j
          break

    # all 0 based, inclusive
    keep_from_global = envfrom_global if envfrom_global != -1 else 0
    keep_to_global = min(non_polya_global_terminal_pos, envto_global) if envto_global != -1 else non_polya_global_terminal_pos

    if keep_from_global >= keep_to_global: continue
    keep_seq = record.seq[keep_from_global:keep_to_global+1].replace('-', '').replace('.', '').upper()

  out_fh.write(f'>{record.id}\n{keep_seq}\n')

out_fh.close()





