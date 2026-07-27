#!/usr/bin/env bash
set -euo pipefail

PANMAP=/scratch1/alan/panmap/build/bin/panmap
GARGAMMEL=/home/alan/tools/gargammel/gargammel.pl
LEEHOME=/home/alan/tools/leeHom/src/leeHom
IDS_FILE=/scratch1/alan/lab_notebook/human_mito/mitoleaf/mito_fastas/merged.ids.txt
FASTA_FILE=/scratch1/alan/lab_notebook/human_mito/mitoleaf/mito_fastas/merged.with_chimp.masked.fasta
PANMAN_FILE=/scratch1/alan/lab_notebook/human_mito/mitoleaf/mito.panman

GARGDIR=$(dirname "$GARGAMMEL")
ADPTSIM=$GARGDIR/src/adptSim
ART=$GARGDIR/art_src_MountRainier/art_illumina
ADAPTF=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG
ADAPTR=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT

NS=(100 200 300 400 500 1000)

usage() {
  echo "usage: $0 OUT_DIR P1 [P2 ... PN]" >&2
  echo "  P1..PN are haplotype proportions in (0,1] and must sum to 1.0" >&2
  exit 1
}

[ $# -ge 2 ] || usage
OUT_DIR=$1
shift
PROPS=("$@")
NHAP=${#PROPS[@]}

# ===== ARGUMENT VALIDATION =====
for p in "${PROPS[@]}"; do
  if ! awk -v x="$p" 'BEGIN{ if (x+0 != x || x <= 0 || x > 1) exit 1 }' 2>/dev/null; then
    echo "error: proportion '$p' is not a number in (0,1]" >&2
    exit 1
  fi
done

SUM=$(printf '%s\n' "${PROPS[@]}" | awk '{s+=$1} END{printf "%.12f", s}')
if ! awk -v s="$SUM" 'BEGIN{ if (s < 1-1e-6 || s > 1+1e-6) exit 1 }'; then
  echo "error: proportions sum to $SUM, expected 1.0 (tolerance 1e-6)" >&2
  exit 1
fi

AVAIL=$(grep -c . "$IDS_FILE")
if [ "$AVAIL" -lt "$NHAP" ]; then
  echo "error: requested $NHAP haplotypes but $IDS_FILE only has $AVAIL" >&2
  exit 1
fi

# frac <= 0.4 rounds down, otherwise up
round_reads() {
  awk -v p="$1" -v n="$2" 'BEGIN{ v=p*n; f=v-int(v); printf "%d", (f<=0.4 ? int(v) : int(v)+1) }'
}

# ===== PRE-FLIGHT: no haplotype may round to zero reads at any depth =====
BAD=""
for n in "${NS[@]}"; do
  for i in "${!PROPS[@]}"; do
    r=$(round_reads "${PROPS[$i]}" "$n")
    if [ "$r" -lt 1 ]; then
      BAD="${BAD}  n=${n}\tproportion[$((i+1))]=${PROPS[$i]}\t-> ${r} reads\n"
    fi
  done
done
if [ -n "$BAD" ]; then
  echo "error: the following (depth, proportion) pairs yield zero reads:" >&2
  printf "$BAD" >&2
  echo "fix the proportions or drop the affected depths from NS before rerunning." >&2
  exit 1
fi

# ===== SETUP =====
mkdir -p "$OUT_DIR"
TMP_DIR=$OUT_DIR/tmp
rm -rf "$OUT_DIR/ref" "$TMP_DIR"
mkdir -p "$OUT_DIR/ref" "$TMP_DIR" "$OUT_DIR/og_fragments"

mapfile -t HAPS < <(shuf -n "$NHAP" "$IDS_FILE")

for hap in "${HAPS[@]}"; do
  mkdir -p "$OUT_DIR/ref/${hap}_ref/endo" "$OUT_DIR/ref/${hap}_ref/cont" "$OUT_DIR/ref/${hap}_ref/bact"
  seqkit grep -p "$hap" "$FASTA_FILE" > "$OUT_DIR/ref/${hap}_ref/endo/$hap.fa"
done

RUN_SALT=$(date +%s%N)
MANIFEST=$OUT_DIR/manifest.tsv
printf 'n_total\thaplotype\tproportion\tn_reads\tseed\tmixed_fq\tmixed_nodmg_fq\n' > "$MANIFEST"

derive_seed() {
  local digits
  digits=$(printf '%s:%s:%s' "$RUN_SALT" "$1" "$2" | sha256sum | tr -dc '0-9')
  printf '%s' "${digits:0:9}"
}

# ===== MAIN LOOP =====
for n in "${NS[@]}"; do
  MIXED=$OUT_DIR/mixed_${n}_merged.fq.gz
  MIXED_NODMG=$OUT_DIR/mixed_${n}_nodmg_merged.fq.gz

  DMG_PARTS=()
  NODMG_PARTS=()

  for i in "${!HAPS[@]}"; do
    hap=${HAPS[$i]}
    prop=${PROPS[$i]}
    reads=$(round_reads "$prop" "$n")
    SEED=$(derive_seed "$n" "$hap")
    PREFIX=$TMP_DIR/${hap}_${n}

    # ===== DAMAGED ARM =====
    $GARGAMMEL -n "$reads" --minsize 30 --loc 3.7344 --scale 0.35 \
      -damage 0.024,0.36,0.0097,0.68 -rs "$SEED" \
      --o "$PREFIX" "$OUT_DIR/ref/${hap}_ref" > /dev/null 2>&1

    $LEEHOME --ancientdna -fq1 "${PREFIX}_s1.fq.gz" -fq2 "${PREFIX}_s2.fq.gz" \
      -fqo "${PREFIX}_merged" > /dev/null 2>&1

    cp "${PREFIX}.e.fa.gz" "$OUT_DIR/og_fragments/${hap}_${n}.og_fragment.fa.gz"

    # ===== UNDAMAGED ARM =====
    $ADPTSIM -f "$ADAPTF" -s "$ADAPTR" -l 75 -artp "${PREFIX}_nodmg_a.fa" \
      "${PREFIX}.e.fa.gz" > /dev/null 2>&1
    $ART -ss HS25 -amp -na -p -i "${PREFIX}_nodmg_a.fa" -l 75 -c 1 -qs 0 -qs2 0 \
      -rs "$SEED" -o "${PREFIX}_nodmg_s" > /dev/null 2>&1
    gzip -f "${PREFIX}_nodmg_s1.fq" "${PREFIX}_nodmg_s2.fq"

    $LEEHOME --ancientdna -fq1 "${PREFIX}_nodmg_s1.fq.gz" -fq2 "${PREFIX}_nodmg_s2.fq.gz" \
      -fqo "${PREFIX}_nodmg_merged" > /dev/null 2>&1

    DMG_PARTS+=("${PREFIX}_merged.fq.gz")
    NODMG_PARTS+=("${PREFIX}_nodmg_merged.fq.gz")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$n" "$hap" "$prop" "$reads" "$SEED" "$MIXED" "$MIXED_NODMG" >> "$MANIFEST"
  done

  # ===== POOL =====
  zcat "${DMG_PARTS[@]}" | gzip -c > "$MIXED"
  zcat "${NODMG_PARTS[@]}" | gzip -c > "$MIXED_NODMG"

  # ===== PANMAP ON THE MIXTURE =====
  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 -o "$OUT_DIR/mixed_${n}_panmap" \
    "$PANMAN_FILE" "$MIXED" > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 -o "$OUT_DIR/mixed_${n}_maskends_panmap" \
    "$PANMAN_FILE" "$MIXED" --mask-read-ends 5 > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 -o "$OUT_DIR/mixed_${n}_nodmg_panmap" \
    "$PANMAN_FILE" "$MIXED_NODMG" > /dev/null 2>&1
done

rm -rf "$TMP_DIR"