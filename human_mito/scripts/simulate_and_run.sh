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

OUT_DIR=$1

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



# for n in 500; do
for n in 20 40 60 80 100 200 300 400 500 1000; do
  SEED=$(date +%s%N | sha256sum | tr -dc '0-9' | head -c 9)
  PREFIX=$OUT_DIR/${haplotype}_$n

  # ===== DAMAGED ARM =====
  # Using gargammel to simulate reads (seeded so the undamaged arm can match).
  $GARGAMMEL -n $n --minsize 30 --loc 3.7344 --scale 0.35 -damage 0.024,0.36,0.0097,0.68 -rs $SEED --o ${PREFIX} $OUT_DIR/ref 2> /dev/null 2>&1

  # Using leeHom to merge reads.
  $LEEHOME --ancientdna -fq1 ${PREFIX}_s1.fq.gz -fq2 ${PREFIX}_s2.fq.gz -fqo ${PREFIX}_merged > /dev/null 2>&1
  
  cp ${PREFIX}.e.fa.gz $OUT_DIR/og_fragments/${haplotype}_${n}.og_fragment.fa.gz

  # align the reads to the reference and generate the damage signal for sanity check.
  $ALIGN_ANCIENT $OUT_DIR/ref/endo/$haplotype.fa ${PREFIX}_merged.fq.gz $OUT_DIR/bam_files/${haplotype}_${n}_to_ref > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 -o ${PREFIX}_panmap $PANMAN_FILE ${PREFIX}_merged.fq.gz > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 -o ${PREFIX}_maskends_panmap $PANMAN_FILE ${PREFIX}_merged.fq.gz --mask-read-ends 5 > /dev/null 2>&1

  # ===== UNDAMAGED ARM =====
  # Re-sequence the raw (undamaged) fragments gargammel wrote to .e.fa.gz, same adapters/platform/seed.
  $ADPTSIM -f $ADAPTF -s $ADAPTR -l 75 -artp ${PREFIX}_nodmg_a.fa ${PREFIX}.e.fa.gz > /dev/null 2>&1
  $ART -ss HS25 -amp -na -p -i ${PREFIX}_nodmg_a.fa -l 75 -c 1 -qs 0 -qs2 0 -rs $SEED -o ${PREFIX}_nodmg_s > /dev/null 2>&1
  gzip -f ${PREFIX}_nodmg_s1.fq ${PREFIX}_nodmg_s2.fq

  $LEEHOME --ancientdna -fq1 ${PREFIX}_nodmg_s1.fq.gz -fq2 ${PREFIX}_nodmg_s2.fq.gz -fqo ${PREFIX}_nodmg_merged > /dev/null 2>&1

  $ALIGN_ANCIENT $OUT_DIR/ref/endo/$haplotype.fa ${PREFIX}_nodmg_merged.fq.gz $OUT_DIR/bam_files/${haplotype}_${n}_nodmg_to_ref > /dev/null 2>&1

  $PANMAP --meta -t 4 -k 15 -s 8 -l 1 -o ${PREFIX}_nodmg_panmap $PANMAN_FILE ${PREFIX}_nodmg_merged.fq.gz > /dev/null 2>&1
done

# clean up unneeded files
rm $OUT_DIR/*fa.gz $OUT_DIR/*fail.fq.gz $OUT_DIR/*r1.fq.gz $OUT_DIR/*r2.fq.gz $OUT_DIR/*s1.fq.gz $OUT_DIR/*s2.fq.gz $OUT_DIR/*_nodmg_a.fa
