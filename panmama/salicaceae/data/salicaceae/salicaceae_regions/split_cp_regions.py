#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

REGIONS = ("LSC", "IRb", "SSC", "IRa")

def parse_coord(field):
  label, spans = field.split(":", 1)
  intervals = []
  for span in spans.split(","):
    start, end = span.split("-", 1)
    intervals.append((int(start), int(end)))
  return label, intervals

def read_single_fasta(path):
  header = None
  seq_parts = []
  with open(path, "r") as fh:
    for line in fh:
      if line.startswith(">"):
        if header is not None:
          raise ValueError(f"{path} contains multiple records; expected one.")
        header = line.rstrip("\n")
      else:
        seq_parts.append(line.strip())
  if header is None:
    raise ValueError(f"{path} has no FASTA header.")
  return header, "".join(seq_parts)

def write_fasta(path, header, seq, width=80):
  with open(path, "w") as fh:
    fh.write(header + "\n")
    for i in range(0, len(seq), width):
      fh.write(seq[i:i + width] + "\n")

def validate_span(coords, seq_len, lineno):
  all_intervals = []
  for label, ivs in coords.items():
    for start, end in ivs:
      all_intervals.append((start, end, label))
  all_intervals.sort(key=lambda x: x[0])

  if all_intervals[0][0] != 1:
    print(f"[warn] line {lineno}: intervals do not start at 1 "
          f"(first interval begins at {all_intervals[0][0]}); skipping record.",
          file=sys.stderr)
    return False

  if all_intervals[-1][1] != seq_len:
    print(f"[warn] line {lineno}: intervals do not end at sequence length "
          f"(last interval ends at {all_intervals[-1][1]}, seq length {seq_len}); "
          f"skipping record.", file=sys.stderr)
    return False

  for i in range(len(all_intervals) - 1):
    cur_end = all_intervals[i][1]
    nxt_start = all_intervals[i + 1][0]
    if nxt_start != cur_end + 1:
      if nxt_start <= cur_end:
        print(f"[warn] line {lineno}: intervals overlap between "
              f"{cur_end} and {nxt_start}; skipping record.", file=sys.stderr)
      else:
        print(f"[warn] line {lineno}: gap between intervals at "
              f"{cur_end + 1}-{nxt_start - 1}; skipping record.", file=sys.stderr)
      return False

  total = sum(end - start + 1 for start, end, _ in all_intervals)
  if total != seq_len:
    print(f"[warn] line {lineno}: interval lengths sum to {total} "
          f"but sequence length is {seq_len}; skipping record.", file=sys.stderr)
    return False

  return True

def process(coord_file, input_dir, output_root, line_width=80):
  input_dir = Path(input_dir)
  output_root = Path(output_root)
  out_dirs = {r: output_root / r for r in REGIONS}
  for d in out_dirs.values():
    d.mkdir(parents=True, exist_ok=True)

  with open(coord_file, "r") as fh:
    for lineno, line in enumerate(fh, 1):
      line = line.strip()
      if not line:
        continue
      fields = line.split()
      if len(fields) < 5:
        print(f"[warn] line {lineno}: expected 5 columns, got {len(fields)}; skipping.",
              file=sys.stderr)
        continue

      fname = fields[0]
      fasta_path = input_dir / fname
      if not fasta_path.is_file():
        print(f"[warn] line {lineno}: {fasta_path} not found; skipping.", file=sys.stderr)
        continue

      prefix = fname
      for ext in (".fasta", ".fa", ".fna"):
        if prefix.endswith(ext):
          prefix = prefix[: -len(ext)]
          break

      try:
        header, seq = read_single_fasta(fasta_path)
      except ValueError as e:
        print(f"[warn] line {lineno}: {e}; skipping.", file=sys.stderr)
        continue

      coords = {}
      ok = True
      for f in fields[1:5]:
        try:
          label, intervals = parse_coord(f)
        except ValueError:
          print(f"[warn] line {lineno}: cannot parse '{f}'; skipping record.",
                file=sys.stderr)
          ok = False
          break
        if label not in out_dirs:
          print(f"[warn] line {lineno}: unknown region '{label}'; skipping record.",
                file=sys.stderr)
          ok = False
          break
        bad = False
        for start, end in intervals:
          if start < 1 or end > len(seq) or start > end:
            print(f"[warn] line {lineno}: {label} coords {start}-{end} out of range "
                  f"(seq length {len(seq)}); skipping record.", file=sys.stderr)
            bad = True
            break
        if bad:
          ok = False
          break
        coords[label] = intervals
      if not ok:
        continue

      if set(coords.keys()) != set(REGIONS):
        missing = set(REGIONS) - set(coords.keys())
        print(f"[warn] line {lineno}: missing regions {missing}; skipping record.",
              file=sys.stderr)
        continue

      if not validate_span(coords, len(seq), lineno):
        continue

      for label, intervals in coords.items():
        sub = "".join(seq[start - 1:end] for start, end in intervals)
        coord_str = ",".join(f"{s}-{e}" for s, e in intervals)
        new_header = f">{prefix}.{label} {header[1:]} [{coord_str}]"
        out_path = out_dirs[label] / f"{prefix}.{label}.fasta"
        write_fasta(out_path, new_header, sub, width=line_width)

def main():
  ap = argparse.ArgumentParser(description="Split chloroplast FASTAs into LSC/IRb/SSC/IRa regions.")
  ap.add_argument("coord_file", help="TSV/whitespace-delimited file with filename + 4 region coords.")
  ap.add_argument("input_dir", help="Directory containing the input FASTA files.")
  ap.add_argument("output_dir", help="Root output directory; LSC/IRb/SSC/IRa subdirs created inside.")
  ap.add_argument("--width", type=int, default=80, help="FASTA line wrap width (default 80).")
  args = ap.parse_args()
  process(args.coord_file, args.input_dir, args.output_dir, line_width=args.width)

if __name__ == "__main__":
  main()