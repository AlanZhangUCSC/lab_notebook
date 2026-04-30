#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass
class RegionSet:
    lsc: Tuple[int, int]
    irb: Tuple[int, int]
    ssc: Tuple[int, int]
    ira: Tuple[int, int]


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    current = None
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line)
    return {k: "".join(v).upper() for k, v in records.items()}


def write_fasta(path: Path, records: Dict[str, str], wrap: int = 80) -> None:
    with path.open("w") as out:
        for sid, seq in records.items():
            out.write(f">{sid}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i : i + wrap] + "\n")


def informative_bases(seq: str) -> int:
    return sum(1 for c in seq if c in "ACGT")


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")
    return seq.translate(table)[::-1]


def parse_cpstools_regions(path: Path) -> RegionSet:
    text = path.read_text()
    pat = re.compile(r"(LSC|IRb|SSC|IRa)\s*[:=]\s*(\d+)\s*-\s*(\d+)", re.I)
    found: Dict[str, Tuple[int, int]] = {}
    for label, start, end in pat.findall(text):
        found[label.upper()] = (int(start), int(end))
    missing = {"LSC", "IRB", "SSC", "IRA"} - set(found)
    if missing:
        raise ValueError(f"Missing cpstools regions {sorted(missing)} in {path}")
    return RegionSet(found["LSC"], found["IRB"], found["SSC"], found["IRA"])


def rotate_to_lsc(seq: str, regions: RegionSet) -> Tuple[str, RegionSet]:
    shift0 = regions.lsc[0] - 1
    rotated = seq[shift0:] + seq[:shift0]
    length = len(seq)

    def shift(r: Tuple[int, int]) -> Tuple[int, int]:
        s, e = r
        ns = ((s - 1 - shift0) % length) + 1
        ne = ((e - 1 - shift0) % length) + 1
        if ne < ns:
            raise ValueError(f"Wrapped interval after LSC rotation: {r}")
        return ns, ne

    return rotated, RegionSet(
        lsc=shift(regions.lsc),
        irb=shift(regions.irb),
        ssc=shift(regions.ssc),
        ira=shift(regions.ira),
    )


def parse_regions_tsv(path: Path) -> Dict[str, RegionSet]:
    regions: Dict[str, RegionSet] = {}
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("sample\t"):
                continue
            left, *parts = line.split("\t")
            sample = left.replace("salicaceae.part_", "").replace(".fasta", "")
            region_map: Dict[str, Tuple[int, int]] = {}
            for part in parts:
                if ":" not in part:
                    continue
                name, coords = part.split(":")
                start, end = coords.split("-")
                region_map[name.upper()] = (int(start), int(end))
            regions[sample] = RegionSet(
                lsc=region_map["LSC"],
                irb=region_map["IRB"],
                ssc=region_map["SSC"],
                ira=region_map["IRA"],
            )
    return regions


def slice_region(seq: str, coords: Tuple[int, int]) -> str:
    start, end = coords
    return seq[start - 1 : end]


def normalize_mafft_headers(path: Path) -> None:
    records = read_fasta(path)
    normalized: Dict[str, str] = {}
    for sid, seq in records.items():
        nid = sid[3:] if sid.startswith("_R_") else sid
        old = normalized.get(nid)
        if old is None or informative_bases(seq) > informative_bases(old):
            normalized[nid] = seq
    write_fasta(path, normalized)


def seq_column_map(aln_seq: str) -> List[int | None]:
    coords: List[int | None] = []
    pos = 0
    for ch in aln_seq:
        if ch in "ACGTNacgtn":
            pos += 1
            coords.append(pos)
        else:
            coords.append(None)
    return coords


def columns_to_keep_by_trimmed(aln_records: Dict[str, str], trimmed_records: Dict[str, str]) -> List[int]:
    sample_ids = sorted(set(aln_records) & set(trimmed_records))
    if not sample_ids:
        raise ValueError("No shared sample IDs between aligned and trimmed FASTA")
    aln_len = len(next(iter(aln_records.values())))
    trim_len = len(next(iter(trimmed_records.values())))
    aln_cols = [tuple(aln_records[sid][i] for sid in sample_ids) for i in range(aln_len)]
    trim_cols = [tuple(trimmed_records[sid][i] for sid in sample_ids) for i in range(trim_len)]
    keep_cols: List[int] = []
    j = 0
    for i, col in enumerate(aln_cols):
        if j < trim_len and col == trim_cols[j]:
            keep_cols.append(i)
            j += 1
    if j != trim_len:
        raise ValueError("Unable to map trimmed columns back to original alignment")
    return keep_cols


def contiguous_ranges(indices: Iterable[int]) -> List[Tuple[int, int]]:
    idx = sorted(set(indices))
    if not idx:
        return []
    out: List[Tuple[int, int]] = []
    start = prev = idx[0]
    for x in idx[1:]:
        if x == prev + 1:
            prev = x
            continue
        out.append((start, prev))
        start = prev = x
    out.append((start, prev))
    return out


def make_cactus_chunk_id(sample: str, partition: str, start: int, end: int, chunk_num: int) -> str:
    # Cactus only accepts header first words made of alnum plus "_-:."
    return f"{sample}:{partition}:{start}-{end}:chunk{chunk_num}"


def run(cmd: Sequence[str], stdout_path: Path | None = None) -> None:
    stdout_handle = stdout_path.open("w") if stdout_path else None
    try:
        subprocess.run(list(cmd), check=True, text=True, stdout=stdout_handle)
    finally:
        if stdout_handle:
            stdout_handle.close()


def parse_maf_blocks(maf_path: Path) -> List[List[Tuple[str, int, int, str, int, str]]]:
    blocks: List[List[Tuple[str, int, int, str, int, str]]] = []
    current: List[Tuple[str, int, int, str, int, str]] = []
    for raw in maf_path.open():
        line = raw.strip()
        if not line:
            if current:
                blocks.append(current[:])
                current = []
            continue
        if not line.startswith("s "):
            continue
        cols = line.split()
        sample = cols[1].split(".", 1)[0]
        current.append((sample, int(cols[2]) + 1, int(cols[3]), cols[4], int(cols[5]), cols[6]))
    if current:
        blocks.append(current[:])
    return blocks


def parse_xmfa_blocks(xmfa_path: Path) -> List[List[Tuple[str, int, int, str, int, str]]]:
    blocks: List[List[Tuple[str, int, int, str, int, str]]] = []
    seq_index_to_sample: Dict[str, str] = {}
    current_block: List[Tuple[str, int, int, str, int, str]] = []
    current_header: Tuple[str, int, int, str] | None = None
    current_seq: List[str] = []

    def flush_record() -> None:
        nonlocal current_header, current_seq, current_block
        if current_header is None:
            return
        sample, start, end, strand = current_header
        text = "".join(current_seq).upper()
        size = len(text.replace("-", ""))
        current_block.append((sample, min(start, end), size, strand, 0, text))
        current_header = None
        current_seq = []

    header_re = re.compile(r"^>\s*(\d+):(\d+)-(\d+)\s+([+-])(?:\s+(.+))?$")
    seqfile_re = re.compile(r"^#Sequence(\d+)File\s+(.+)$")
    for raw in xmfa_path.open():
        line = raw.strip()
        if not line:
            continue
        m = seqfile_re.match(line)
        if m:
            idx, path_text = m.groups()
            seq_index_to_sample[idx] = Path(path_text).stem
            continue
        if line.startswith("#"):
            continue
        if line == "=":
            flush_record()
            if current_block:
                blocks.append(current_block[:])
                current_block = []
            continue
        if line.startswith(">"):
            flush_record()
            m = header_re.match(line)
            if not m:
                raise ValueError(f"Unrecognized XMFA header: {line}")
            idx, start, end, strand, path_text = m.groups()
            sample = seq_index_to_sample.get(idx)
            if sample is None:
                if not path_text:
                    raise ValueError(f"Cannot resolve sample name for XMFA header: {line}")
                sample = Path(path_text).stem
            current_header = (sample, int(start), int(end), strand)
            current_seq = []
            continue
        current_seq.append(line)
    flush_record()
    if current_block:
        blocks.append(current_block[:])
    return blocks


def orient_by_outgroup(block_rows: List[Tuple[str, int, int, str, int, str]], outgroup: str) -> List[Tuple[str, int, int, str, int, str]]:
    for row in block_rows:
        if row[0] == outgroup and row[3] == "-":
            oriented = []
            for sample, start, size, strand, src_size, text in block_rows:
                oriented.append((sample, start, size, "+" if strand == "-" else "-", src_size, revcomp(text)))
            return oriented
    return block_rows


def cmd_orient_outgroup(args: argparse.Namespace) -> None:
    seq = next(iter(read_fasta(args.input_fasta).values()))
    regions = parse_cpstools_regions(args.cpstools_txt)
    rotated_seq, rotated_regions = rotate_to_lsc(seq, regions)
    write_fasta(args.output_fasta, {args.sample_id: rotated_seq})
    with args.output_regions.open("w") as out:
        out.write("sample\tLSC\tIRb\tSSC\tIRa\n")
        out.write(
            f"{args.sample_id}\t"
            f"LSC:{rotated_regions.lsc[0]}-{rotated_regions.lsc[1]}\t"
            f"IRb:{rotated_regions.irb[0]}-{rotated_regions.irb[1]}\t"
            f"SSC:{rotated_regions.ssc[0]}-{rotated_regions.ssc[1]}\t"
            f"IRa:{rotated_regions.ira[0]}-{rotated_regions.ira[1]}\n"
        )


def cmd_split_partitions(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    sample_regions = parse_regions_tsv(args.sample_regions_tsv)
    outgroup_regions = parse_regions_tsv(args.outgroup_regions_tsv)[args.outgroup_id]
    outgroup_seq = next(iter(read_fasta(args.outgroup_fasta).values()))
    lsc_records: Dict[str, str] = {}
    irb_records: Dict[str, str] = {}
    ssc_records: Dict[str, str] = {}
    for sample, region_set in sorted(sample_regions.items()):
        fasta_path = args.sample_fasta_dir / f"{sample}.fasta"
        seq = next(iter(read_fasta(fasta_path).values()))
        lsc_records[sample] = slice_region(seq, region_set.lsc)
        irb_records[sample] = slice_region(seq, region_set.irb)
        ssc_records[sample] = slice_region(seq, region_set.ssc)
    lsc_records[args.outgroup_id] = slice_region(outgroup_seq, outgroup_regions.lsc)
    irb_records[args.outgroup_id] = slice_region(outgroup_seq, outgroup_regions.irb)
    ssc_records[args.outgroup_id] = slice_region(outgroup_seq, outgroup_regions.ssc)
    write_fasta(args.outdir / "LSC.fasta", lsc_records)
    write_fasta(args.outdir / "IRb.fasta", irb_records)
    write_fasta(args.outdir / "SSC.fasta", ssc_records)


def cmd_normalize_mafft(args: argparse.Namespace) -> None:
    normalize_mafft_headers(args.input_fasta)


def cmd_concat_backbone(args: argparse.Namespace) -> None:
    lsc = read_fasta(args.lsc_trimmed)
    irb = read_fasta(args.irb_trimmed)
    ssc = read_fasta(args.ssc_trimmed)
    all_ids = sorted(set(lsc) | set(irb) | set(ssc))
    l_len = len(next(iter(lsc.values())))
    i_len = len(next(iter(irb.values())))
    s_len = len(next(iter(ssc.values())))
    spacer = "N" * args.spacer
    records = {
        sid: lsc.get(sid, "N" * l_len) + spacer + irb.get(sid, "N" * i_len) + spacer + ssc.get(sid, "N" * s_len)
        for sid in all_ids
    }
    write_fasta(args.output_fasta, records)
    with args.output_partitions.open("w") as out:
        out.write(f"DNA,LSC=1-{l_len}\n")
        out.write(f"DNA,IRb={l_len + args.spacer + 1}-{l_len + args.spacer + i_len}\n")
        start3 = l_len + args.spacer + i_len + args.spacer + 1
        out.write(f"DNA,SSC={start3}-{start3 + s_len - 1}\n")


def cmd_recover_chunks(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    partitions = [
        ("LSC", args.lsc_aligned, args.lsc_trimmed),
        ("IRb", args.irb_aligned, args.irb_trimmed),
        ("SSC", args.ssc_aligned, args.ssc_trimmed),
    ]
    per_sample_chunks: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    all_samples = set()
    meta_path = args.outdir / "recovered_chunks.tsv"
    with meta_path.open("w") as meta:
        meta.write("sample\tsource_partition\tsource_start\tsource_end\tstrand\tungapped_length\tchunk_id\n")
        for label, aln_path, trim_path in partitions:
            aln_records = read_fasta(aln_path)
            trim_records = read_fasta(trim_path)
            all_samples.update(aln_records)
            keep_cols = set(columns_to_keep_by_trimmed(aln_records, trim_records))
            aln_len = len(next(iter(aln_records.values())))
            removed_cols = [i for i in range(aln_len) if i not in keep_cols]
            for sample, aln_seq in aln_records.items():
                coord_map = seq_column_map(aln_seq)
                chunk_counter = 0
                for a, b in contiguous_ranges(removed_cols):
                    fragment = aln_seq[a : b + 1]
                    ungapped = "".join(ch for ch in fragment if ch in "ACGTN")
                    ungapped_no_n = "".join(ch for ch in ungapped if ch in "ACGT")
                    if len(ungapped_no_n) <= args.min_len:
                        continue
                    valid_coords = [coord_map[i] for i in range(a, b + 1) if coord_map[i] is not None]
                    if not valid_coords:
                        continue
                    chunk_counter += 1
                    chunk_id = make_cactus_chunk_id(
                        sample,
                        label,
                        min(valid_coords),
                        max(valid_coords),
                        chunk_counter,
                    )
                    per_sample_chunks[sample].append((chunk_id, ungapped))
                    meta.write(
                        f"{sample}\t{label}\t{min(valid_coords)}\t{max(valid_coords)}\t+\t{len(ungapped_no_n)}\t{chunk_id}\n"
                    )
    for sample in sorted(all_samples):
        sample_records = {chunk_id: seq for chunk_id, seq in per_sample_chunks.get(sample, [])}
        if not sample_records:
            sample_records = {make_cactus_chunk_id(sample, "empty", 0, 0, 0): "N"}
        write_fasta(args.outdir / f"{sample}.fa", sample_records)


def cmd_write_seqfile(args: argparse.Namespace) -> None:
    samples = sorted(p.stem for p in args.recovered_dir.glob("*.fa"))
    tree_text = args.guide_tree.read_text().strip()
    with args.output_seqfile.open("w") as out:
        out.write(tree_text + "\n")
        for sid in samples:
            out.write(f"{sid}\t{args.recovered_dir / (sid + '.fa')}\n")


def cmd_extract_blocks(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    mafft = shutil.which(args.mafft_bin)
    if not mafft:
        raise FileNotFoundError("mafft not found in PATH for block extraction")
    manifest = args.outdir / "recovered_blocks.tsv"
    block_idx = 0
    if args.xmfa:
        parsed_blocks = parse_xmfa_blocks(args.xmfa)
    elif args.maf:
        parsed_blocks = parse_maf_blocks(args.maf)
    else:
        raise ValueError("One of --xmfa or --maf is required")
    with manifest.open("w") as out:
        out.write("block_id\tsort_key\toutgroup_start\tsample_count\talignment\n")
        for block in parsed_blocks:
            rows = orient_by_outgroup(block, args.outgroup)
            present = {r[0] for r in rows}
            if len(present) <= args.min_samples:
                continue
            block_records: Dict[str, str] = {}
            sample_meta: Dict[str, int] = {}
            for sample, start, size, strand, src_size, text in rows:
                seq = text.replace("-", "").upper()
                if len(re.sub(r"[^ACGT]", "", seq)) == 0:
                    continue
                old = block_records.get(sample)
                if old is None or len(seq) > len(old):
                    block_records[sample] = seq
                    sample_meta[sample] = start
            if len(block_records) < 2:
                continue
            block_idx += 1
            out_start = sample_meta.get(args.outgroup, 10**12)
            if out_start == 10**12:
                fallback = min(block_records)
                out_start = sample_meta[fallback]
                sort_key = f"fallback|{fallback}|{out_start:012d}"
            else:
                sort_key = f"outgroup|{out_start:012d}"
            prefix = args.outdir / f"block_{block_idx:05d}"
            raw = prefix.with_suffix(".input.fasta")
            aln = prefix.with_suffix(".aln.fasta")
            write_fasta(raw, block_records)
            run(
                [
                    mafft,
                    "--thread",
                    str(args.threads),
                    "--6merpair",
                    "--retree",
                    "1",
                    "--maxiterate",
                    "0",
                    "--adjustdirection",
                    str(raw),
                ],
                stdout_path=aln,
            )
            normalize_mafft_headers(aln)
            chosen = aln
            if not chosen.exists() or chosen.stat().st_size == 0:
                continue
            aln_records = read_fasta(chosen)
            if not aln_records:
                continue
            aln_len = len(next(iter(aln_records.values())))
            if aln_len <= args.min_block_len:
                continue
            out.write(f"block_{block_idx:05d}\t{sort_key}\t{out_start}\t{len(block_records)}\t{chosen}\n")


def cmd_concat_final(args: argparse.Namespace) -> None:
    backbone = read_fasta(args.backbone_fasta)
    all_ids = sorted(backbone)
    spacer = "N" * args.spacer
    concat = {sid: backbone[sid] for sid in all_ids}
    parts: List[Tuple[str, int, int]] = [("backbone", 1, len(next(iter(backbone.values()))))]
    cursor = parts[0][2]
    rows = []
    with args.blocks_manifest.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = sorted(reader, key=lambda r: (r["sort_key"], r["block_id"]))
    for row in rows:
        aln_path = Path(row["alignment"])
        if not aln_path.exists() or aln_path.stat().st_size == 0:
            stem = aln_path.stem.replace(".trim", "")
            alt_aln = aln_path.with_name(f"{stem}.aln.fasta")
            if alt_aln.exists() and alt_aln.stat().st_size > 0:
                aln_path = alt_aln
            else:
                continue
        seqs = read_fasta(aln_path)
        if not seqs:
            continue
        blen = len(next(iter(seqs.values())))
        for sid in all_ids:
            concat[sid] += spacer + seqs.get(sid, "-" * blen)
        start = cursor + args.spacer + 1
        end = start + blen - 1
        parts.append((row["block_id"], start, end))
        cursor = end
    write_fasta(args.output_fasta, concat)
    with args.output_partitions.open("w") as out:
        for name, start, end in parts:
            out.write(f"DNA,{name}={start}-{end}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Helper subcommands for stepwise Salicaceae cp pipeline.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("orient-outgroup")
    p.add_argument("--input-fasta", type=Path, required=True)
    p.add_argument("--cpstools-txt", type=Path, required=True)
    p.add_argument("--output-fasta", type=Path, required=True)
    p.add_argument("--output-regions", type=Path, required=True)
    p.add_argument("--sample-id", default="NC_034285.1")
    p.set_defaults(func=cmd_orient_outgroup)

    p = sub.add_parser("split-partitions")
    p.add_argument("--sample-regions-tsv", type=Path, required=True)
    p.add_argument("--sample-fasta-dir", type=Path, required=True)
    p.add_argument("--outgroup-fasta", type=Path, required=True)
    p.add_argument("--outgroup-regions-tsv", type=Path, required=True)
    p.add_argument("--outgroup-id", default="NC_034285.1")
    p.add_argument("--outdir", type=Path, required=True)
    p.set_defaults(func=cmd_split_partitions)

    p = sub.add_parser("normalize-mafft")
    p.add_argument("--input-fasta", type=Path, required=True)
    p.set_defaults(func=cmd_normalize_mafft)

    p = sub.add_parser("concat-backbone")
    p.add_argument("--lsc-trimmed", type=Path, required=True)
    p.add_argument("--irb-trimmed", type=Path, required=True)
    p.add_argument("--ssc-trimmed", type=Path, required=True)
    p.add_argument("--output-fasta", type=Path, required=True)
    p.add_argument("--output-partitions", type=Path, required=True)
    p.add_argument("--spacer", type=int, default=100)
    p.set_defaults(func=cmd_concat_backbone)

    p = sub.add_parser("recover-chunks")
    p.add_argument("--lsc-aligned", type=Path, required=True)
    p.add_argument("--lsc-trimmed", type=Path, required=True)
    p.add_argument("--irb-aligned", type=Path, required=True)
    p.add_argument("--irb-trimmed", type=Path, required=True)
    p.add_argument("--ssc-aligned", type=Path, required=True)
    p.add_argument("--ssc-trimmed", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--min-len", type=int, default=100)
    p.set_defaults(func=cmd_recover_chunks)

    p = sub.add_parser("write-seqfile")
    p.add_argument("--recovered-dir", type=Path, required=True)
    p.add_argument("--guide-tree", type=Path, required=True)
    p.add_argument("--output-seqfile", type=Path, required=True)
    p.set_defaults(func=cmd_write_seqfile)

    p = sub.add_parser("extract-blocks")
    p.add_argument("--maf", type=Path)
    p.add_argument("--xmfa", type=Path)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--outgroup", default="NC_034285.1")
    p.add_argument("--mafft-bin", default="mafft")
    p.add_argument("--threads", type=int, default=12)
    p.add_argument("--min-samples", type=int, default=2, help="Keep blocks with sample count > this threshold.")
    p.add_argument("--min-block-len", type=int, default=100, help="Keep blocks with aligned length > this threshold.")
    p.set_defaults(func=cmd_extract_blocks)

    p = sub.add_parser("concat-final")
    p.add_argument("--backbone-fasta", type=Path, required=True)
    p.add_argument("--blocks-manifest", type=Path, required=True)
    p.add_argument("--output-fasta", type=Path, required=True)
    p.add_argument("--output-partitions", type=Path, required=True)
    p.add_argument("--spacer", type=int, default=100)
    p.set_defaults(func=cmd_concat_final)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
