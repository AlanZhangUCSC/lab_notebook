from __future__ import annotations
import os
import subprocess
import threading
from collections import namedtuple
import pysam
from bisect import bisect_right
from collections import namedtuple
from itertools import islice
from operator import attrgetter
import tempfile
from dataclasses import dataclass, field
import numpy as np



RefFragment = namedtuple(
  "RefFragment", "query_name ref_name start end is_reverse ref_seq"
)

MergedFragment = namedtuple(
  "MergedFragment", "ref_name start end ref_seq n_reads query_names"
)

_BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
_ALN_OPTS = ["-t", "1", "-l", "1024", "-n", "0.02"]


def _drain(stream, sink):
  sink.append(stream.read())
  stream.close()


def _ensure_index(ref_fasta_file):
  if not os.path.exists(ref_fasta_file):
    raise FileNotFoundError(ref_fasta_file)

  paths = [ref_fasta_file + s for s in _BWA_INDEX_SUFFIXES]
  if all(os.path.exists(p) for p in paths):
    ref_mtime = os.path.getmtime(ref_fasta_file)
    if all(os.path.getmtime(p) >= ref_mtime for p in paths):
      return

  proc = subprocess.run(["bwa", "index", ref_fasta_file],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  if proc.returncode != 0:
    raise RuntimeError(
      f"bwa index failed (exit {proc.returncode}): "
      f"{proc.stderr.decode(errors='replace')[-2000:]}"
    )


def get_covered_ref_fragments(ref_fasta_file, fastq_file):
  _ensure_index(ref_fasta_file)

  procs, errs, drains = [], {}, []

  def spawn(cmd, stdin=None):
    p = subprocess.Popen(cmd, stdin=stdin, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    procs.append((cmd[0] + " " + cmd[1], p))
    errs[p] = []
    t = threading.Thread(target=_drain, args=(p.stderr, errs[p]), daemon=True)
    t.start()
    drains.append(t)
    return p

  aln = spawn(["bwa", "aln"] + _ALN_OPTS + [ref_fasta_file, fastq_file])
  # bwa samse takes no "-" for the .sai argument; /dev/stdin is fd 0, the pipe.
  samse = spawn(["bwa", "samse", ref_fasta_file, "/dev/stdin", fastq_file],
                stdin=aln.stdout)
  aln.stdout.close()
  view = spawn(["samtools", "view", "-u", "-F", "4"], stdin=samse.stdout)
  samse.stdout.close()

  fragments = []
  try:
    with pysam.AlignmentFile(view.stdout, "rb") as bam:
      for rec in bam:
        fragments.append(RefFragment(
          rec.query_name, rec.reference_name,
          rec.reference_start + 1, rec.reference_end,
          rec.is_reverse, rec.get_reference_sequence().upper(),
        ))
  finally:
    for _, p in procs:
      p.wait()
    for t in drains:
      t.join()

  for name, p in procs:
    if p.returncode != 0:
      raise RuntimeError(
        f"{name} failed (exit {p.returncode}): "
        f"{b''.join(errs[p]).decode(errors='replace')[-2000:]}"
      )

  fragments.sort(key=lambda frag: frag.start)
  return fragments

class FragmentConflictError(ValueError):
  def __init__(self, ref_name, coord, held, held_base, incoming, incoming_base):
    self.ref_name = ref_name
    self.coord = coord
    self.held = held
    self.held_base = held_base
    self.incoming = incoming
    self.incoming_base = incoming_base
    super().__init__(
      f"conflicting reference bases on {ref_name} at 1-based position {coord}: "
      f"read {held.query_name!r} ([{held.start}, {held.end}], "
      f"{'-' if held.is_reverse else '+'}) gives '{held_base}' at ref_seq index "
      f"{coord - held.start}; read {incoming.query_name!r} "
      f"([{incoming.start}, {incoming.end}], "
      f"{'-' if incoming.is_reverse else '+'}) gives '{incoming_base}' at "
      f"ref_seq index {coord - incoming.start}"
    )


def _validate(frag):
  if frag.start < 1:
    raise ValueError(
      f"{frag.ref_name}: read {frag.query_name!r} has start {frag.start} "
      f"below 1 (coordinates are 1-based inclusive)"
    )
  if frag.end < frag.start:
    raise ValueError(
      f"{frag.ref_name}: read {frag.query_name!r} has end {frag.end} "
      f"preceding start {frag.start}"
    )
  span = frag.end - frag.start + 1
  if len(frag.ref_seq) != span:
    raise ValueError(
      f"{frag.ref_name}: read {frag.query_name!r} spans [{frag.start}, "
      f"{frag.end}] ({span} bases) but carries a ref_seq of length "
      f"{len(frag.ref_seq)}; check that the MD tag is present and consistent"
    )


def _emit(ref_name, start, end, cluster, members):
  return MergedFragment(
    ref_name, start, end, cluster.decode('ascii'),
    len(members), tuple(f.query_name for f in members)
  )


def _merge_one_ref(ref_name, group):
  out = []
  first = group[0]
  cur_start, cur_end = first.start, first.end
  cluster = bytearray(first.ref_seq.encode('ascii'))
  part_offsets = [0]
  part_frags = [first]
  members = [first]

  for frag in islice(group, 1, None):
    if frag.start > cur_end + 1:
      out.append(_emit(ref_name, cur_start, cur_end, cluster, members))
      cur_start, cur_end = frag.start, frag.end
      cluster = bytearray(frag.ref_seq.encode('ascii'))
      part_offsets = [0]
      part_frags = [frag]
      members = [frag]
      continue

    incoming = frag.ref_seq.encode('ascii')
    overlap = min(frag.end, cur_end) - frag.start + 1
    off = frag.start - cur_start

    if overlap and cluster[off:off + overlap] != incoming[:overlap]:
      i = next(k for k in range(overlap) if cluster[off + k] != incoming[k])
      held = part_frags[bisect_right(part_offsets, off + i) - 1]
      raise FragmentConflictError(
        ref_name, frag.start + i, held, chr(cluster[off + i]),
        frag, chr(incoming[i])
      )

    members.append(frag)
    if frag.end > cur_end:
      part_offsets.append(len(cluster))
      part_frags.append(frag)
      cluster += incoming[overlap:]
      cur_end = frag.end

  out.append(_emit(ref_name, cur_start, cur_end, cluster, members))
  return out


def merge_ref_fragments(fragments):
  by_ref = {}
  for frag in fragments:
    _validate(frag)
    by_ref.setdefault(frag.ref_name, []).append(frag)

  merged = []
  for ref_name, group in by_ref.items():
    if any(a.start > b.start for a, b in zip(group, islice(group, 1, None))):
      group.sort(key=attrgetter('start'))
    merged.extend(_merge_one_ref(ref_name, group))
  return merged


_ALN_OPTS = ["-t", "1", "-l", "1024", "-n", "0.02"]
_VIEW_CMD = ["samtools", "view", "-u", "-F", "4"]

_ENC = bytes(
  {65: 1, 97: 1, 67: 2, 99: 2, 71: 3, 103: 3, 84: 4, 116: 4}.get(i, 5)
  for i in range(256)
)
_AMBIGUOUS = 5

_BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

_CONSUMES_BOTH = frozenset((0, 7, 8))
_CONSUMES_QUERY = frozenset((1, 4))
_CONSUMES_REF = frozenset((2, 3))

# per reference position; the update is a monotone max, so a MATCHED position
# can never be pulled back down to MISMATCHED by a later read
_UNCOVERED = 0
_AMBIG = 1
_MISMATCHED = 2
_MATCHED = 3


@dataclass
class DivergenceResult:
  substitutions: int = 0
  covered_positions: int = 0
  ambiguous_positions: int = 0
  inserted_bases: int = 0
  deleted_bases: int = 0
  reads_used: int = 0
  conflicting_observations: int = 0

  @property
  def divergence(self) -> float:
    return self.substitutions / self.covered_positions if self.covered_positions else 0.0


def align_and_measure_divergence(reference_path: str, query_path: str) -> DivergenceResult:
  _ensure_indexes(reference_path)

  res = DivergenceResult()
  with tempfile.TemporaryDirectory() as tmp:
    sai = os.path.join(tmp, "query.sai")
    with open(sai, "wb") as fh:
      subprocess.run(
        ["bwa", "aln", *_ALN_OPTS, reference_path, query_path],
        stdout=fh, check=True,
      )

    samse = subprocess.Popen(
      ["bwa", "samse", reference_path, sai, query_path],
      stdout=subprocess.PIPE,
    )
    view = subprocess.Popen(_VIEW_CMD, stdin=samse.stdout, stdout=subprocess.PIPE)
    samse.stdout.close()

    try:
      with pysam.FastaFile(reference_path) as fasta, \
           pysam.AlignmentFile(view.stdout, "rb") as bam:
        _accumulate(bam, fasta, res)
    finally:
      for proc in (view, samse):
        if proc.poll() is None:
          proc.terminate()
        proc.wait()

    if samse.returncode != 0:
      raise subprocess.CalledProcessError(samse.returncode, "bwa samse")
    if view.returncode != 0:
      raise subprocess.CalledProcessError(view.returncode, "samtools view")

  return res


def _ensure_indexes(reference_path: str) -> None:
  if not os.path.exists(reference_path):
    raise FileNotFoundError(reference_path)
  ref_dir = os.path.dirname(os.path.abspath(reference_path)) or "."

  # a partially written index from an interrupted run is unusable, so rebuild
  # the whole set unless every companion file is present
  needs_bwa = any(not os.path.exists(reference_path + s) for s in _BWA_INDEX_SUFFIXES)
  needs_fai = not os.path.exists(reference_path + ".fai")
  if not (needs_bwa or needs_fai):
    return
  if not os.access(ref_dir, os.W_OK):
    raise PermissionError(
      f"{ref_dir} is not writable; cannot build the missing index for {reference_path}"
    )

  if needs_bwa:
    subprocess.run(["bwa", "index", reference_path], check=True)
  if needs_fai:
    pysam.faidx(reference_path)


def _accumulate(bam, fasta, res: DivergenceResult) -> None:
  contigs: dict[int, tuple[np.ndarray, np.ndarray]] = {}

  for read in bam:
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
      continue
    seq = read.query_sequence
    cig = read.cigartuples
    if not seq or not cig:
      continue
    res.reads_used += 1

    tid = read.reference_id
    entry = contigs.get(tid)
    if entry is None:
      raw = fasta.fetch(bam.get_reference_name(tid)).encode("ascii").translate(_ENC)
      entry = (np.frombuffer(raw, dtype=np.uint8), np.zeros(len(raw), dtype=np.uint8))
      contigs[tid] = entry

    query = np.frombuffer(seq.encode("ascii").translate(_ENC), dtype=np.uint8)
    qpos = 0
    rpos = read.reference_start

    for op, length in cig:
      if op in _CONSUMES_BOTH:
        _merge_block(entry, query, qpos, rpos, length, read, res)
        qpos += length
        rpos += length
      elif op in _CONSUMES_QUERY:
        if op == 1:
          res.inserted_bases += length
        qpos += length
      elif op in _CONSUMES_REF:
        if op == 2:
          res.deleted_bases += length
        rpos += length

  _tally(contigs, res)


def _merge_block(entry, query, qpos, rpos, length, read, res: DivergenceResult) -> None:
  ref, state = entry
  qb = query[qpos:qpos + length]
  rb = ref[rpos:rpos + length]
  if qb.size != length or rb.size != length:
    raise ValueError(f"alignment of {read.query_name} runs past the end of the contig")

  usable = (qb != _AMBIGUOUS) & (rb != _AMBIGUOUS)
  value = np.where(usable, np.where(qb == rb, _MATCHED, _MISMATCHED), _AMBIG).astype(np.uint8)

  st = state[rpos:rpos + length]
  res.conflicting_observations += int(np.count_nonzero(
    ((st == _MISMATCHED) & (value == _MATCHED)) | ((st == _MATCHED) & (value == _MISMATCHED))
  ))
  # basic slicing yields a view, so out=st writes through to the backing array
  np.maximum(st, value, out=st)


def _tally(contigs, res: DivergenceResult) -> None:
  for _, state in contigs.values():
    counts = np.bincount(state, minlength=4)
    res.ambiguous_positions += int(counts[_AMBIG])
    res.substitutions += int(counts[_MISMATCHED])
    res.covered_positions += int(counts[_MISMATCHED] + counts[_MATCHED])