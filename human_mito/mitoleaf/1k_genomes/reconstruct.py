import sys, csv
from collections import defaultdict

def load_fasta(path):
  seq = []
  with open(path) as f:
    for line in f:
      if not line.startswith('>'):
        seq.append(line.strip())
  return ''.join(seq).upper()

def load_variants(tsv):
  by_sample = defaultdict(list)
  with open(tsv, newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
      by_sample[row['sample_id']].append(
        (int(row['pos']), row['type'], row['ref'], row['alt']))
  return by_sample

def validate_ref_bases(ref, by_sample):
  L = len(ref)
  mism = []
  checked = 0
  for sample, variants in by_sample.items():
    for pos, typ, r, a in variants:
      if typ == 'INS':
        continue
      if not (1 <= pos <= L):
        mism.append((sample, pos, typ, r, a, 'OUT_OF_RANGE')); continue
      checked += 1
      obs = ref[pos - 1]
      if obs != r:
        mism.append((sample, pos, typ, r, a, f'rCRS_has_{obs}'))
  return checked, mism

def reconstruct(ref, variants):
  cells = [''] + list(ref)
  L = len(ref)
  for pos, typ, r, a in variants:
    if typ == 'INS' or not (1 <= pos <= L):
      continue
    cells[pos] = '' if typ == 'DEL' else a
  # insertions append to their anchor so downstream rCRS coords never shift
  for pos, typ, r, a in variants:
    if typ == 'INS' and 1 <= pos <= L:
      cells[pos] = cells[pos] + a
  return ''.join(cells)

def wrap(s, w=70):
  return '\n'.join(s[i:i+w] for i in range(0, len(s), w))

def main():
  if len(sys.argv) < 4:
    sys.exit("usage: reconstruct.py rCRS.fasta variants.tsv out.fasta [--nonstrict]")
  ref_path, tsv_path, out_path = sys.argv[1:4]
  strict = '--nonstrict' not in sys.argv
  ref = load_fasta(ref_path)
  by_sample = load_variants(tsv_path)

  checked, mism = validate_ref_bases(ref, by_sample)
  print(f"[sanity] reference length: {len(ref)}")
  print(f"[sanity] substitution/deletion sites checked: {checked}")
  print(f"[sanity] ref-base mismatches: {len(mism)}")
  if mism:
    with open('reconstruct_mismatches.tsv', 'w', newline='') as f:
      w = csv.writer(f, delimiter='\t')
      w.writerow(['sample', 'pos', 'type', 'claimed_ref', 'alt', 'note'])
      w.writerows(mism)
    for row in mism[:5]:
      print(f"    {row[0]} {row[3]}{row[1]}{row[4]} -> {row[5]}")
    print("    (full list in reconstruct_mismatches.tsv)")
    if strict:
      sys.exit("[sanity] FAILED: claimed ref bases disagree with rCRS. "
               "Fix coordinates/reference, or rerun with --nonstrict to skip bad edits.")

  n = 0
  with open(out_path, 'w') as out:
    for sample, variants in by_sample.items():
      seq = reconstruct(ref, variants)
      n_ins = sum(1 for v in variants if v[1] == 'INS')
      n_del = sum(1 for v in variants if v[1] == 'DEL')
      assert len(seq) == len(ref) + n_ins - n_del, sample
      out.write(f">{sample}\n{wrap(seq)}\n")
      n += 1
  print(f"[done] wrote {n} consensus sequences to {out_path}")

if __name__ == '__main__':
  main()