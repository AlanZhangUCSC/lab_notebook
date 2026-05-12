#!/usr/bin/env python3
import sys
import argparse

def parse_assemblies(path):
  table = {}
  with open(path) as f:
    for line in f:
      line = line.rstrip("\n")
      if not line:
        continue
      fields = line.split("\t")
      if len(fields) < 3:
        continue
      asm = fields[0]
      common = fields[1]
      sci = fields[2].split()
      genus = sci[0] if len(sci) >= 1 else "NA"
      species = sci[1] if len(sci) >= 2 else "NA"
      table[asm] = (common, genus, species)
  return table

def parse_nhmm(path):
  table = {}
  with open(path) as f:
    for line in f:
      line = line.rstrip("\n")
      if not line:
        continue
      fields = line.split("\t")
      if len(fields) < 2:
        continue
      key = fields[0]
      paren = key.rfind("(")
      if paren != -1:
        key = key[:paren]
      fam_field = fields[1].split("#", 1)[0]
      table[key] = fam_field
  return table

def header_to_nhmm_key(header):
  parts = header.split("__")
  if len(parts) < 3:
    return None
  asm = parts[0]
  fam = parts[1]
  loc = "__".join(parts[2:])
  us = loc.rfind("_")
  if us == -1:
    return None
  loc_converted = loc[:us] + ":" + loc[us+1:]
  return f"{asm}::{fam}::{loc_converted}"

def emit(out, header, seq_len, asm_table, nhmm_table):
  if header is None:
    return
  parts = header.split("__")
  asm = parts[0] if parts else "NA"
  header_fam = parts[1] if len(parts) >= 2 else "NA"

  if asm in asm_table:
    common, genus, species = asm_table[asm]
  else:
    sys.stderr.write(f"WARNING: assembly '{asm}' not found in assemblies file (header: {header})\n")
    common, genus, species = "NA", "NA", "NA"

  key = header_to_nhmm_key(header)
  if key is None:
    sys.stderr.write(f"WARNING: could not construct nhmm key from header '{header}'\n")
    nhmm_fam = "NA"
  elif key in nhmm_table:
    nhmm_fam = nhmm_table[key]
  else:
    sys.stderr.write(f"WARNING: key '{key}' not found in nhmmscan file (header: {header})\n")
    nhmm_fam = "NA"

  out.write(f"{header}\t{asm}\t{common}\t{genus}\t{species}\t{seq_len}\t{header_fam}\t{nhmm_fam}\n")

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("--fasta", required=True)
  ap.add_argument("--assemblies", required=True)
  ap.add_argument("--nhmm", required=True)
  ap.add_argument("--out", default="-")
  ap.add_argument("--no-header", action="store_true")
  args = ap.parse_args()

  asm_table = parse_assemblies(args.assemblies)
  nhmm_table = parse_nhmm(args.nhmm)

  out = sys.stdout if args.out == "-" else open(args.out, "w")
  try:
    if not args.no_header:
      out.write("seq_id\tassembly\tcommon_name\tgenus\tspecies\tlength\tlabeled_family\thmm_family\n")

    cur_header = None
    cur_len = 0
    with open(args.fasta) as f:
      for line in f:
        if not line:
          continue
        if line[0] == ">":
          emit(out, cur_header, cur_len, asm_table, nhmm_table)
          cur_header = line[1:].rstrip("\n").rstrip("\r").split()[0]
          cur_len = 0
        else:
          cur_len += len(line.rstrip("\n").rstrip("\r"))
      emit(out, cur_header, cur_len, asm_table, nhmm_table)
  finally:
    if out is not sys.stdout:
      out.close()

if __name__ == "__main__":
  main()