PANMAP=/scratch1/alan/panmap/build/bin/panmap
GARGAMMEL=/home/alan/tools/gargammel/gargammel.pl
LEEHOME=/home/alan/tools/leeHom/src/leeHom
ALIGN_ANCIENT=/home/alan/tools/misc/align_ancient
IDS_FILE=/scratch1/alan/lab_notebook/human_mito/mitoleaf/mito_fastas/merged.ids.txt
FASTA_FILE=/scratch1/alan/lab_notebook/human_mito/mitoleaf/mito_fastas/merged.with_chimp.masked.fasta
PANMAN_FILE=/scratch1/alan/lab_notebook/human_mito/mitoleaf/mito.panman


GARGDIR=$(dirname "$GARGAMMEL")
ADPTSIM=$GARGDIR/src/adptSim
ART=$GARGDIR/art_src_MountRainier/art_illumina
SEED=13
ADAPTF=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG
ADAPTR=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT

# Argument parsing: support both positional and flags, error on missing OUT_DIR
usage() {
  echo "Usage: $0 OUT_DIR [-v DAMAGE_V] [-l DAMAGE_L] [-d DAMAGE_D] [-s DAMAGE_S]"
  echo "  OUT_DIR           Output directory (required)"
  echo "  -v, --v VALUE     Damage parameter v (default: 0.024)"
  echo "  -l, --l VALUE     Damage parameter l (default: 0.36)"
  echo "  -d, --d VALUE     Damage parameter d (default: 0.0097)"
  echo "  -s, --s VALUE     Damage parameter s (default: 0.68)"
  exit 1
}

# Default values
V=0.024
L=0.36
D=0.0097
S=0.68

OUT_DIR=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--v)
      if [[ -z "$2" ]]; then usage; fi
      V="$2"
      shift 2
      ;;
    -l|--l)
      if [[ -z "$2" ]]; then usage; fi
      L="$2"
      shift 2
      ;;
    -d|--d)
      if [[ -z "$2" ]]; then usage; fi
      D="$2"
      shift 2
      ;;
    -s|--s)
      if [[ -z "$2" ]]; then usage; fi
      S="$2"
      shift 2
      ;;
    -*)
      echo "Unknown flag: $1"
      usage
      ;;
    *)
      if [[ -z "$OUT_DIR" ]]; then
        OUT_DIR="$1"
      else
        echo "Unknown argument or duplicate OUT_DIR: $1"
        usage
      fi
      shift
      ;;
  esac
done

if [[ -z "$OUT_DIR" ]]; then
  echo "Error: OUT_DIR must be specified."
  usage
fi

mkdir -p $OUT_DIR


# grab a random haplotype 
haplotype=$(shuf -n 1 $IDS_FILE)

# simulate reads across 20, 40, 60, 80, 100, 200, 300, 400, 500, 1000 fragments.
mkdir -p $OUT_DIR/ref
if [ -n "$(ls -A $OUT_DIR/ref 2>/dev/null)" ]; then
  echo "Removing existing reference files..."
  rm -rf $OUT_DIR/ref/*
fi
mkdir $OUT_DIR/ref/endo $OUT_DIR/ref/cont $OUT_DIR/ref/bact
mkdir $OUT_DIR/bam_files
mkdir $OUT_DIR/og_fragments
seqkit grep -p $haplotype $FASTA_FILE > $OUT_DIR/ref/endo/$haplotype.fa



# simulate 1000 fragments once, then subsample down to each n
MAX_N=1000
SEED=$(date +%s%N | sha256sum | tr -dc '0-9' | head -c 9)
BASE=$OUT_DIR/${haplotype}_full

# ===== DAMAGED ARM (generated once) =====
$GARGAMMEL -n $MAX_N --minsize 30 --loc 3.913 --scale 0.35 -damage $V,$L,$D,$S -rs $SEED --o ${BASE} $OUT_DIR/ref 2> /dev/null 2>&1
$LEEHOME --ancientdna -fq1 ${BASE}_s1.fq.gz -fq2 ${BASE}_s2.fq.gz -fqo ${BASE}_merged > /dev/null 2>&1

# ===== UNDAMAGED ARM (generated once) =====
$ADPTSIM -f $ADAPTF -s $ADAPTR -l 75 -artp ${BASE}_nodmg_a.fa ${BASE}.e.fa.gz > /dev/null 2>&1
$ART -ss HS25 -amp -na -p -i ${BASE}_nodmg_a.fa -l 75 -c 1 -qs 0 -qs2 0 -rs $SEED -o ${BASE}_nodmg_s > /dev/null 2>&1
gzip -f ${BASE}_nodmg_s1.fq ${BASE}_nodmg_s2.fq
$LEEHOME --ancientdna -fq1 ${BASE}_nodmg_s1.fq.gz -fq2 ${BASE}_nodmg_s2.fq.gz -fqo ${BASE}_nodmg_merged > /dev/null 2>&1

for n in 20 40 60 80 100 200 300 400 500 1000; do
  PREFIX=$OUT_DIR/${haplotype}_$n
  SUBSEED=$((10#${SEED} % 1000000 + n))

  seqkit sample -s $SUBSEED -n $n ${BASE}_merged.fq.gz -o ${PREFIX}_merged.fq.gz 2> /dev/null

  # propagate the sampled read set to the undamaged arm and to the source fragments
  seqkit seq -n -i ${PREFIX}_merged.fq.gz > ${PREFIX}.read_ids.txt
  seqkit grep -f ${PREFIX}.read_ids.txt ${BASE}_nodmg_merged.fq.gz -o ${PREFIX}_nodmg_merged.fq.gz 2> /dev/null
  sed -E 's/-[0-9]+$//' ${PREFIX}.read_ids.txt | sort -u > ${PREFIX}.frag_ids.txt
  seqkit grep -f ${PREFIX}.frag_ids.txt ${BASE}.e.fa.gz -o $OUT_DIR/og_fragments/${haplotype}_${n}.og_fragment.fa.gz 2> /dev/null

  # ===== DAMAGED ARM ANALYSIS =====
  $ALIGN_ANCIENT $OUT_DIR/ref/endo/$haplotype.fa ${PREFIX}_merged.fq.gz $OUT_DIR/bam_files/${haplotype}_${n}_to_ref > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 --em-convergence-threshold 0.000001 -o ${PREFIX}_panmap $PANMAN_FILE ${PREFIX}_merged.fq.gz > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 --em-convergence-threshold 0.000001 -o ${PREFIX}_maskends_panmap $PANMAN_FILE ${PREFIX}_merged.fq.gz --mask-read-ends 5 > /dev/null 2>&1

  # ===== UNDAMAGED ARM ANALYSIS =====
  $ALIGN_ANCIENT $OUT_DIR/ref/endo/$haplotype.fa ${PREFIX}_nodmg_merged.fq.gz $OUT_DIR/bam_files/${haplotype}_${n}_nodmg_to_ref > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 --em-convergence-threshold 0.000001 -o ${PREFIX}_nodmg_panmap $PANMAN_FILE ${PREFIX}_nodmg_merged.fq.gz > /dev/null 2>&1
done

# clean up unneeded files
rm $OUT_DIR/*fa.gz $OUT_DIR/*fail.fq.gz $OUT_DIR/*r1.fq.gz $OUT_DIR/*r2.fq.gz $OUT_DIR/*s1.fq.gz $OUT_DIR/*s2.fq.gz $OUT_DIR/*_nodmg_a.fa
