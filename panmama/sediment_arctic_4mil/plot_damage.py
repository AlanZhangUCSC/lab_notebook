#!/usr/bin/env python3
"""
Lightweight BAM -> subs + damage/mismatch/length plots, plus aDNA-likeness scoring.

Subcommands:
  compute  BAM -> subs file + combined plot (damage + mismatch freq + read length).
  score    subs file -> per-end damage/background/signal metrics.
"""

import argparse
import math
import os
import sys
from collections import defaultdict

import pysam

try:
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  HAVE_MPL = True
except ImportError:
  HAVE_MPL = False


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def mismatch_table(seq, cigar_tuples, md, is_reverse, max_pos):
  out = []
  if seq is None or md is None or not cigar_tuples:
    return out

  read_bases = []
  ref_bases = []

  md_i = 0
  md_len = len(md)
  seq_i = 0

  def next_md_match():
    nonlocal md_i
    n = 0
    while md_i < md_len and md[md_i].isdigit():
      n = n * 10 + (ord(md[md_i]) - 48)
      md_i += 1
    return n

  pending_match = next_md_match()

  for op, length in cigar_tuples:
    if op == 0 or op == 7 or op == 8:
      remaining = length
      while remaining > 0:
        if pending_match > 0:
          take = pending_match if pending_match < remaining else remaining
          for k in range(take):
            b = seq[seq_i + k]
            read_bases.append(b)
            ref_bases.append(b)
          seq_i += take
          pending_match -= take
          remaining -= take
        else:
          if md_i < md_len and md[md_i] == "^":
            md_i += 1
            while md_i < md_len and md[md_i].isalpha():
              md_i += 1
            pending_match = next_md_match()
            continue
          if md_i >= md_len:
            break
          ref_base = md[md_i]
          md_i += 1
          read_bases.append(seq[seq_i])
          ref_bases.append(ref_base)
          seq_i += 1
          remaining -= 1
          pending_match = next_md_match()
    elif op == 1:
      seq_i += length
    elif op == 2:
      if md_i < md_len and md[md_i] == "^":
        md_i += 1
        while md_i < md_len and md[md_i].isalpha():
          md_i += 1
        pending_match = next_md_match()
    elif op == 4:
      seq_i += length
    elif op == 5:
      pass

  if is_reverse:
    read_bases = [COMPLEMENT.get(b, "N") for b in reversed(read_bases)]
    ref_bases = [COMPLEMENT.get(b, "N") for b in reversed(ref_bases)]

  L = len(read_bases)
  if L == 0:
    return out

  for i in range(min(max_pos, L)):
    pos = i + 1
    rb = ref_bases[i]
    qb = read_bases[i]
    if rb == "N" or qb == "N":
      continue
    out.append((rb, qb, pos))

  for i in range(min(max_pos, L)):
    j = L - 1 - i
    if j < max_pos:
      break
    pos = -(i + 1)
    rb = ref_bases[j]
    qb = read_bases[j]
    if rb == "N" or qb == "N":
      continue
    out.append((rb, qb, pos))

  return out


def gather_all(bam_path, max_pos=15, min_mapq=0,
               dedup_multimappers=False, best_alignment=False):
  """
  Single pass over the BAM producing substitution counts, NM histogram, and
  read-length histogram. Combining them avoids a second IO pass versus
  running damage analysis and plotbaminfo separately.

  Memory: counts <= 480 keys; nm_hist and length_hist are O(distinct values)
  which is small. Everything else is transient per-alignment.

  Multi-mapper policies match the documented semantics in earlier versions:
    default            : every accepted alignment contributes 1.
    dedup_multimappers : each alignment scaled by 1/N for its read.
    best_alignment     : one alignment per read (NM asc, MAPQ desc), counted 1.

  Grouped modes require consecutive alignments per query_name.
  """
  if dedup_multimappers and best_alignment:
    raise ValueError("dedup_multimappers and best_alignment are mutually exclusive")

  counts = defaultdict(float) if dedup_multimappers else defaultdict(int)
  nm_hist = defaultdict(float) if dedup_multimappers else defaultdict(int)
  length_hist = defaultdict(int)

  bam = pysam.AlignmentFile(bam_path, "rb", require_index=False)
  n_reads = 0
  n_skipped = 0

  def emit(aln, weight):
    if not aln.has_tag("MD") or aln.query_sequence is None:
      return False
    md = aln.get_tag("MD")
    for k in mismatch_table(aln.query_sequence, aln.cigartuples, md,
                            aln.is_reverse, max_pos):
      counts[k] += weight
    try:
      nm = aln.get_tag("NM")
      nm_hist[int(nm)] += weight
    except KeyError:
      pass
    return True

  if not dedup_multimappers and not best_alignment:
    for aln in bam:
      if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
      if aln.mapping_quality < min_mapq:
        continue
      if not aln.has_tag("MD") or aln.query_sequence is None:
        n_skipped += 1
        continue
      if not emit(aln, 1):
        n_skipped += 1
        continue
      length_hist[aln.query_length or 0] += 1
      n_reads += 1
    bam.close()
    return counts, nm_hist, length_hist, n_reads, n_skipped

  if dedup_multimappers:
    current_name = None
    buffered = []

    def flush():
      nonlocal buffered
      if not buffered:
        return
      w = 1.0 / len(buffered)
      for aln in buffered:
        emit(aln, w)
      length_hist[buffered[0].query_length or 0] += 1
      buffered = []

    for aln in bam:
      if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
      if aln.mapping_quality < min_mapq:
        continue
      if not aln.has_tag("MD") or aln.query_sequence is None:
        n_skipped += 1
        continue
      name = aln.query_name
      if name != current_name:
        flush()
        current_name = name
        n_reads += 1
      buffered.append(aln)
    flush()
    bam.close()
    return counts, nm_hist, length_hist, n_reads, n_skipped

  current_name = None
  best_rank = None
  best_aln = None

  def commit():
    if best_aln is None:
      return
    emit(best_aln, 1)
    length_hist[best_aln.query_length or 0] += 1

  for aln in bam:
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
      continue
    if aln.mapping_quality < min_mapq:
      continue
    if not aln.has_tag("MD") or aln.query_sequence is None:
      n_skipped += 1
      continue
    try:
      nm = aln.get_tag("NM")
    except KeyError:
      n_skipped += 1
      continue

    name = aln.query_name
    if name != current_name:
      commit()
      current_name = name
      best_rank = None
      best_aln = None
      n_reads += 1

    rank = (nm, -aln.mapping_quality)
    if best_rank is None or rank < best_rank:
      best_rank = rank
      best_aln = aln

  commit()
  bam.close()
  return counts, nm_hist, length_hist, n_reads, n_skipped


def write_subs_file(counts, out_path, max_pos=15, tax_id=0, tax_name="all"):
  totals = defaultdict(float)
  for (f, t, p), c in counts.items():
    totals[(f, p)] += c

  bases = ("A", "C", "G", "T")
  tokens = []
  for p in list(range(1, max_pos + 1)) + list(range(-1, -max_pos - 1, -1)):
    for f in bases:
      denom = totals.get((f, p), 0)
      if denom == 0:
        continue
      for t in bases:
        c = counts.get((f, t, p), 0)
        if c == 0:
          continue
        tokens.append(f"{f}{t}{p}:{c / denom:.6f}")

  with open(out_path, "w") as fh:
    fh.write(f"{tax_id}\t{tax_name}\t{' '.join(tokens)}\n")


def curves_from_counts(counts, max_pos=15):
  def zeros():
    return [0.0] * max_pos

  ct5 = zeros(); ga5 = zeros(); oth5 = zeros()
  ct3 = zeros(); ga3 = zeros(); oth3 = zeros()

  c5 = [0.0] * max_pos; g5 = [0.0] * max_pos; nd5 = [0.0] * max_pos
  c3 = [0.0] * max_pos; g3 = [0.0] * max_pos; nd3 = [0.0] * max_pos
  ct5r = [0.0] * max_pos; ga5r = [0.0] * max_pos; oth5r = [0.0] * max_pos
  ct3r = [0.0] * max_pos; ga3r = [0.0] * max_pos; oth3r = [0.0] * max_pos

  for (f, t, p), c in counts.items():
    if 1 <= p <= max_pos:
      i = p - 1
      if f == "C":
        c5[i] += c
        if t == "T": ct5r[i] += c
      if f == "G":
        g5[i] += c
        if t == "A": ga5r[i] += c
      if f == t:
        nd5[i] += c
      elif not ((f == "C" and t == "T") or (f == "G" and t == "A")):
        oth5r[i] += c
        nd5[i] += c
    elif -max_pos <= p <= -1:
      i = -p - 1
      if f == "C":
        c3[i] += c
        if t == "T": ct3r[i] += c
      if f == "G":
        g3[i] += c
        if t == "A": ga3r[i] += c
      if f == t:
        nd3[i] += c
      elif not ((f == "C" and t == "T") or (f == "G" and t == "A")):
        oth3r[i] += c
        nd3[i] += c

  for i in range(max_pos):
    ct5[i] = ct5r[i] / c5[i] if c5[i] else 0.0
    ga5[i] = ga5r[i] / g5[i] if g5[i] else 0.0
    oth5[i] = oth5r[i] / nd5[i] if nd5[i] else 0.0
    ct3[i] = ct3r[i] / c3[i] if c3[i] else 0.0
    ga3[i] = ga3r[i] / g3[i] if g3[i] else 0.0
    oth3[i] = oth3r[i] / nd3[i] if nd3[i] else 0.0

  return {"ct5": ct5, "ga5": ga5, "oth5": oth5,
          "ct3": ct3, "ga3": ga3, "oth3": oth3}


def parse_subs_file(path, max_pos=15):
  """
  Parse a bamdam-format subs file into per-line curves. Note that the subs
  file stores per-(from,to,pos) proportions rather than raw counts, so the
  Other rate reconstructed here is the sum of non-C>T / non-G>A mismatch
  rates at that position rather than the joint non-deam denominator used
  when counts are available. For the terminus-vs-interior score this
  distinction doesn't matter as long as Other is computed the same way
  on both ends, which it is.
  """
  results = []
  with open(path) as fh:
    for line in fh:
      line = line.rstrip("\n")
      if not line:
        continue
      parts = line.split("\t")
      if len(parts) < 3:
        continue
      tax_id, tax_name, data = parts[0], parts[1], parts[2]

      ct5 = [0.0] * max_pos; ga5 = [0.0] * max_pos; oth5 = [0.0] * max_pos
      ct3 = [0.0] * max_pos; ga3 = [0.0] * max_pos; oth3 = [0.0] * max_pos

      for tok in data.split():
        try:
          key, val = tok.split(":")
          val = float(val)
        except ValueError:
          continue
        if len(key) < 3:
          continue
        f = key[0]; t = key[1]
        try:
          p = int(key[2:])
        except ValueError:
          continue
        if f == t:
          continue
        if 1 <= p <= max_pos:
          i = p - 1
          if f == "C" and t == "T":
            ct5[i] = val
          elif f == "G" and t == "A":
            ga5[i] = val
          else:
            oth5[i] += val
        elif -max_pos <= p <= -1:
          i = -p - 1
          if f == "C" and t == "T":
            ct3[i] = val
          elif f == "G" and t == "A":
            ga3[i] = val
          else:
            oth3[i] += val

      results.append((tax_id, tax_name, {
        "ct5": ct5, "ga5": ga5, "oth5": oth5,
        "ct3": ct3, "ga3": ga3, "oth3": oth3}))
  return results


def ancient_score(curves, max_pos=15, interior_start=5):
  """
  Terminus-vs-interior damage metric.

  signal_5 = CT[1] - mean(Other_5'[interior_start..max_pos])
  signal_3 = GA[-1] - mean(Other_3'[interior_start..max_pos])

  In a genuine aDNA library the termini are enriched for cytosine
  deamination while the interior Other rate reflects baseline sequencing
  error and true divergence. Subtracting the interior background corrects
  for references of different quality. Double-stranded libraries are
  expected to show signal_5 approximately equal to signal_3; asymmetric
  signatures are suggestive of single-stranded protocols, contamination,
  or artefactual spikes at one end.

  Verdict thresholds are tuned for non-UDG libraries and are intentionally
  conservative. UDG or USER-treated libraries will read lower (terminus
  C>T of order 0.02-0.1 rather than 0.1-0.4); users processing such
  libraries should rely on the numeric signal rather than the verdict.
  """
  ct5 = curves["ct5"]; ga3 = curves["ga3"]
  oth5 = curves["oth5"]; oth3 = curves["oth3"]

  interior_5 = oth5[interior_start - 1:]
  interior_3 = oth3[interior_start - 1:]
  bg5 = sum(interior_5) / len(interior_5) if interior_5 else 0.0
  bg3 = sum(interior_3) / len(interior_3) if interior_3 else 0.0

  term5 = ct5[0]
  term3 = ga3[0]
  sig5 = term5 - bg5
  sig3 = term3 - bg3

  if sig5 > 0 and sig3 > 0:
    ratio = max(sig5, sig3) / min(sig5, sig3)
  else:
    ratio = float("inf")

  if sig5 > 0.05 and sig3 > 0.05 and ratio < 2.0:
    verdict = "HIGH"
  elif sig5 > 0.03 or sig3 > 0.03:
    verdict = "MEDIUM"
  else:
    verdict = "LOW"

  return {
    "term_ct5": term5, "term_ga3": term3,
    "bg_5": bg5, "bg_3": bg3,
    "signal_5": sig5, "signal_3": sig3,
    "symmetry_ratio": ratio, "verdict": verdict,
  }


def make_plot(counts, nm_hist, length_hist, out_path, title,
              max_pos=15, ymax=0.0):
  if not HAVE_MPL:
    print("matplotlib not available; skipping plot.", file=sys.stderr)
    return

  curves = curves_from_counts(counts, max_pos)
  xs_5 = list(range(1, max_pos + 1))
  xs_3 = list(range(-max_pos, 0))

  ct_5 = curves["ct5"]
  ga_5 = curves["ga5"]
  oth_5 = curves["oth5"]
  ct_3 = list(reversed(curves["ct3"]))
  ga_3 = list(reversed(curves["ga3"]))
  oth_3 = list(reversed(curves["oth3"]))

  if ymax <= 0:
    flat = ct_5 + ga_5 + oth_5 + ct_3 + ga_3 + oth_3
    max_y = min(1.0, max(max(flat) * 1.2, 1e-3)) if flat else 0.1
  else:
    max_y = float(ymax)

  colors = {"Other": "#009E73", "CT": "#F8766D", "GA": "#56B4E9"}
  nm_color = "#009E73"
  len_color = "#CC79A7"

  fig = plt.figure(figsize=(12, 8))

  ax1 = plt.subplot(2, 2, 1)
  ax1.plot(xs_5, ct_5, color=colors["CT"], linewidth=2, label="C to T")
  ax1.plot(xs_5, ga_5, color=colors["GA"], linewidth=2, label="G to A")
  ax1.plot(xs_5, oth_5, color=colors["Other"], linewidth=2, label="Other")
  ax1.set_xlabel("Position (5')")
  ax1.set_ylabel("Frequency")
  ax1.set_ylim(0, max_y)
  ax1.set_xticks(xs_5)
  ax1.set_title("Damage (5')")

  ax2 = plt.subplot(2, 2, 2)
  ax2.plot(xs_3, ct_3, color=colors["CT"], linewidth=2, label="C to T")
  ax2.plot(xs_3, ga_3, color=colors["GA"], linewidth=2, label="G to A")
  ax2.plot(xs_3, oth_3, color=colors["Other"], linewidth=2, label="Other")
  ax2.set_xlabel("Position (3')")
  ax2.set_ylim(0, max_y)
  ax2.yaxis.tick_right()
  ax2.set_xticks(xs_3)
  ax2.set_title("Damage (3')")
  ax2.legend(loc="upper right")

  ax3 = plt.subplot(2, 2, 3)
  if nm_hist:
    xs = sorted(nm_hist.keys())
    ys = [nm_hist[x] for x in xs]
    ax3.plot(xs, ys, color=nm_color, linewidth=2)
  ax3.set_xlabel("Number of mismatches (NM)")
  ax3.set_ylabel("Reads")
  ax3.set_title("Mismatch frequency")

  ax4 = plt.subplot(2, 2, 4)
  if length_hist:
    xs = sorted(length_hist.keys())
    ys = [length_hist[x] for x in xs]
    ax4.plot(xs, ys, color=len_color, linewidth=2)
  ax4.set_xlabel("Read length")
  ax4.set_ylabel("Reads")
  ax4.set_title("Read length distribution")

  plt.suptitle(title)
  plt.tight_layout(rect=[0, 0, 1, 0.96])

  ext = os.path.splitext(out_path)[1].lower()
  if ext == ".pdf":
    plt.savefig(out_path, format="pdf")
  else:
    plt.savefig(out_path, format="png")
  plt.close(fig)


def print_score_report(scored, out_tsv=None):
  header = ["tax_id", "tax_name", "CT_pos1", "GA_pos-1",
            "bg_5prime", "bg_3prime", "signal_5prime", "signal_3prime",
            "symmetry_ratio", "verdict"]
  rows = []
  for tax_id, tax_name, s in scored:
    rows.append([
      str(tax_id), tax_name,
      f"{s['term_ct5']:.4f}", f"{s['term_ga3']:.4f}",
      f"{s['bg_5']:.4f}", f"{s['bg_3']:.4f}",
      f"{s['signal_5']:.4f}", f"{s['signal_3']:.4f}",
      (f"{s['symmetry_ratio']:.2f}"
       if math.isfinite(s['symmetry_ratio']) else "inf"),
      s["verdict"],
    ])

  lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
  text = "\n".join(lines) + "\n"
  if out_tsv:
    with open(out_tsv, "w") as fh:
      fh.write(text)
    print(f"Wrote score table: {out_tsv}", file=sys.stderr)
  print(text)


def cmd_compute(args):
  if not os.path.exists(args.bam):
    print(f"Error: BAM file not found: {args.bam}", file=sys.stderr)
    sys.exit(1)

  print(f"Scanning {args.bam} ...", file=sys.stderr)
  counts, nm_hist, length_hist, n_reads, n_skipped = gather_all(
    args.bam, args.max_pos, args.min_mapq,
    args.dedup_multimappers, args.best_alignment)
  print(f"  reads used:    {n_reads}", file=sys.stderr)
  print(f"  reads skipped: {n_skipped}", file=sys.stderr)

  if n_reads == 0:
    print("No usable alignments. Ensure BAM has MD tags "
          "(samtools calmd -b in.bam ref.fa > out.bam).", file=sys.stderr)
    sys.exit(2)

  write_subs_file(counts, args.out_subs, args.max_pos, args.tax_id, args.tax_name)
  print(f"Wrote subs file: {args.out_subs}", file=sys.stderr)

  title = args.title if args.title else os.path.basename(args.bam)
  make_plot(counts, nm_hist, length_hist, args.out_plot, title,
            args.max_pos, args.ymax)
  print(f"Wrote plot:     {args.out_plot}", file=sys.stderr)

  if args.score:
    curves = curves_from_counts(counts, args.max_pos)
    s = ancient_score(curves, args.max_pos, args.interior_start)
    print_score_report([(args.tax_id, args.tax_name, s)])


def cmd_score(args):
  entries = parse_subs_file(args.in_subs, args.max_pos)
  if not entries:
    print("No entries in subs file.", file=sys.stderr)
    sys.exit(1)
  scored = []
  for tax_id, tax_name, curves in entries:
    s = ancient_score(curves, args.max_pos, args.interior_start)
    scored.append((tax_id, tax_name, s))
  print_score_report(scored, args.out_tsv)


def build_parser():
  ap = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  sub = ap.add_subparsers(dest="cmd", required=True)

  c = sub.add_parser("compute", help="BAM -> subs + combined plot")
  c.add_argument("--bam", required=True)
  c.add_argument("--out-subs", required=True)
  c.add_argument("--out-plot", required=True)
  c.add_argument("--max-pos", type=int, default=15)
  c.add_argument("--min-mapq", type=int, default=0)
  c.add_argument("--tax-id", type=int, default=0)
  c.add_argument("--tax-name", type=str, default="all")
  c.add_argument("--ymax", type=float, default=0.0)
  c.add_argument("--title", default=None)
  c.add_argument("--score", action="store_true",
                 help="Also print aDNA-likeness score.")
  c.add_argument("--interior-start", type=int, default=5,
                 help="First position (from terminus) treated as interior "
                      "for background estimation (default 5).")
  g = c.add_mutually_exclusive_group()
  g.add_argument("--dedup-multimappers", action="store_true")
  g.add_argument("--best-alignment", action="store_true")
  c.set_defaults(func=cmd_compute)

  s = sub.add_parser("score", help="subs file -> aDNA-likeness metrics")
  s.add_argument("--in-subs", required=True)
  s.add_argument("--out-tsv", default=None)
  s.add_argument("--max-pos", type=int, default=15)
  s.add_argument("--interior-start", type=int, default=5)
  s.add_argument("--tax-id", type=int, default=0)
  s.add_argument("--tax-name", type=str, default="all")
  s.set_defaults(func=cmd_score)

  return ap


def main():
  args = build_parser().parse_args()
  args.func(args)


if __name__ == "__main__":
  main()