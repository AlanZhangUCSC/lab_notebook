#!/bin/bash

set -x

source $(conda info --base)/etc/profile.d/conda.sh
conda deactivate

usage() {
  cat << EOF
Usage: $0 --seqtype SEQTYPE --numhap NUMHAP --numsnp NUMSNP --permut PERMUT --numreads NUMREADS --rep REP --panmap PANMAP --panman PANMAN --pmi PMI --out-prefix DIR [OPTIONS]

Required arguments:
  --seqtype SEQTYPE                 Sequence type
  --numhap NUMHAP                   Number of haplotypes
  --numsnp NUMSNP                   Number of SNPs
  --permut PERMUT                   Permutation parameter
  --numreads NUMREADS               Number of reads
  --rep REP                         Replication parameter
  --panmap PANMAP                   Pan map parameter
  --panman PANMAN                   Pan man parameter
  --pmi PMI                         PMI parameter
  --random-seed STRING              Random seed string
  --out-prefix DIR                  Output prefix
  --cpus                            Number of CPUs to use (default: 1)
  --test-commands                   Print commands to stdout

Optional arguments:
  --swampy PATH                     Path to SWAMPY executable (enables SWAMPY mode)
  --reference-primer-bed-file FILE  Reference primer BED file (required with --swampy)
  --reference-fasta-file FILE       Reference FASTA file (required with --swampy)
  --jvarkit JVARKIT                 JVarKit parameter (required with --swampy)
  -h, --help                        Display this help message

EOF
  exit 1
}

SEQTYPE=""
NUMHAP=""
NUMSNP=""
PERMUT=""
NUMREADS=""
REP=""
PANMAP=""
PANMAN=""
PMI=""
RANDOM_SEED=""
SWAMPY=""
REFERENCE_PRIMER_BED_FILE=""
REFERENCE_FASTA_FILE=""
JVARKIT=""
OUT_PREFIX=""
TEST_COMMANDS=0
CPUS=1

while [[ $# -gt 0 ]]; do
  case $1 in
    --test-commands)
      TEST_COMMANDS=1
      shift
      ;;
    --seqtype)
      SEQTYPE="$2"
      shift 2
      ;;
    --cpus)
      CPUS="$2"
      shift 2
      ;;
    --numhap)
      NUMHAP="$2"
      shift 2
      ;;
    --numsnp)
      NUMSNP="$2"
      shift 2
      ;;
    --permut)
      PERMUT="$2"
      shift 2
      ;;
    --numreads)
      NUMREADS="$2"
      shift 2
      ;;
    --rep)
      REP="$2"
      shift 2
      ;;
    --panmap)
      PANMAP="$2"
      shift 2
      ;;
    --panman)
      PANMAN="$2"
      shift 2
      ;;
    --pmi)
      PMI="$2"
      shift 2
      ;;
    --random-seed)
      RANDOM_SEED="$2"
      shift 2
      ;;
    --swampy)
      SWAMPY="$2"
      shift 2
      ;;
    --reference-primer-bed-file)
      REFERENCE_PRIMER_BED_FILE="$2"
      shift 2
      ;;
    --reference-fasta-file)
      REFERENCE_FASTA_FILE="$2"
      shift 2
      ;;
    --jvarkit)
      JVARKIT="$2"
      shift 2
      ;;
    --out-prefix)
      OUT_PREFIX="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Error: Unknown option $1"
      usage
      ;;
  esac
done

if [[ -z "$SEQTYPE" || -z "$NUMHAP" || -z "$NUMSNP" || -z "$PERMUT" || -z "$NUMREADS" || -z "$REP" || -z "$PANMAP" || -z "$PANMAN" || -z "$PMI" || -z "$RANDOM_SEED" || -z "$OUT_PREFIX" ]]; then
  echo "Error: All of --seqtype, --numhap, --numsnp, --permut, --numreads, --rep, --panmap, --panman, --pmi, and --out-prefix are required."
  usage
fi

if [[ -n "$SWAMPY" ]]; then
  if [[ -z "$REFERENCE_PRIMER_BED_FILE" || -z "$REFERENCE_FASTA_FILE" || -z "$JVARKIT" ]]; then
    echo "Error: --reference-primer-bed-file, --reference-fasta-file, and --jvarkit are required when --swampy is enabled."
    exit 1
  fi
fi

logfile=${OUT_PREFIX}.genreads.log
if [[ -e "$logfile" ]]; then
  rm "$logfile"
fi

SCRIPT_PATH="$(realpath "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

{
echo "Configuration:"
echo "  SEQTYPE: $SEQTYPE"
echo "  NUMHAP: $NUMHAP"
echo "  NUMSNP: $NUMSNP"
echo "  PERMUT: $PERMUT"
echo "  NUMREADS: $NUMREADS"
echo "  REP: $REP"
echo "  PANMAP: $PANMAP"
echo "  PANMAN: $PANMAN"
echo "  PMI: $PMI"
echo "  RANDOM_SEED: $RANDOM_SEED"
echo "  SWAMPY: $SWAMPY"
if [[ -n "$SWAMPY" ]]; then
  echo "  REFERENCE_PRIMER_BED_FILE: $REFERENCE_PRIMER_BED_FILE"
  echo "  REFERENCE_FASTA_FILE: $REFERENCE_FASTA_FILE"
  echo "  JVARKIT: $JVARKIT"
fi
echo "  OUT_PREFIX: $OUT_PREFIX"
} >> "$logfile"
echo "" >> "$logfile"

# Hash the RANDOM_SEED to a 32-bit unsigned int
HASHED_SEED=$(echo -n "$RANDOM_SEED" | md5sum | awk '{print $1}' | xxd -r -p | od -An -tu4 | awk '{print $1}' | awk '{print $1 % 4294967296}')

echo "Hashed seed: $HASHED_SEED" >> "$logfile"




cluster_sizes_1hap=(100)
cluster_sizes_5haps=(50 20 15 10 5)
cluster_sizes_10haps=(15 15 15 15 10 10 5 5 5 5)
cluster_sizes_20haps=($(printf "5 %.0s" {1..20}))
cluster_sizes_50haps=($(printf "2 %.0s" {1..50}))
cluster_sizes_100haps=($(printf "1 %.0s" {1..100}))
cluster_sizes=()
if [[ "$NUMHAP" -eq 1 ]]; then
  cluster_sizes=("${cluster_sizes_1hap[@]}")
elif [[ "$NUMHAP" -eq 5 ]]; then
  cluster_sizes=("${cluster_sizes_5haps[@]}")
elif [[ "$NUMHAP" -eq 10 ]]; then
  cluster_sizes=("${cluster_sizes_10haps[@]}")
elif [[ "$NUMHAP" -eq 20 ]]; then
  cluster_sizes=("${cluster_sizes_20haps[@]}")
elif [[ "$NUMHAP" -eq 50 ]]; then
  cluster_sizes=("${cluster_sizes_50haps[@]}")
elif [[ "$NUMHAP" -eq 100 ]]; then
  cluster_sizes=("${cluster_sizes_100haps[@]}")
else
  echo "Error: Unsupported NUMHAP value: $NUMHAP" >> "$logfile"
  exit 1
fi
cluster_sizes_string="${cluster_sizes[@]}"


abundance_file=${OUT_PREFIX}.abundance.txt


# abundance_file=${OUT_PREFIX}.abundance.txt

# OUT_PREFIX=${OUT_PREFIX}_perfect

# random_nodes_file=${OUT_PREFIX}.randomNodeIDs.txt
# touch $random_nodes_file


OUT_DIR=$(dirname "$OUT_PREFIX")
BASE_DIR=$(basename "$OUT_DIR")
PREFIX_BASE=$(basename "$OUT_PREFIX")

OUT_PREFIX="/data/tmp/${PREFIX_BASE}_perfect"




random_nodes_file=${OUT_PREFIX}.randomNodeIDs.txt
touch $random_nodes_file

cat $abundance_file | cut -f 1 > $random_nodes_file


dump_sequences_parameter_string=$(cat $random_nodes_file | sed -e 's/^/\"/g' -e 's/$/\"/g' | tr '\n' ' ' )
simulate_snps_parameter_string=""
num_mutated=$(awk -v p="$PERMUT" 'BEGIN {print int(100 * p)}')
num_unmutated=$((100 - num_mutated))

for i in $(seq 1 $num_mutated); do
  simulate_snps_parameter_string+=" $NUMSNP"
done

for i in $(seq 1 $num_unmutated); do
  simulate_snps_parameter_string+=" 0"
done

simulate_snps_parameter_string=$(python3 $SCRIPT_DIR/shuff.py $simulate_snps_parameter_string --seed $HASHED_SEED)

cmd=(docker run --rm
  -v "$(realpath $PANMAP):/panmap"
  -v "$(realpath $(dirname $PANMAN)):/panmans"
  -v "$(realpath $(dirname $OUT_PREFIX)):/output"
  -w /panmap
  --user "$(id -u):$(id -g)"
  panmap-dev
  bash -c "/panmap/build/bin/panmap /panmans/$(basename $PANMAN) \
           --dump-sequences $dump_sequences_parameter_string \
           --simulate-snps $simulate_snps_parameter_string  \
           --random-seed $HASHED_SEED \
           --prefix /output/$(basename $OUT_PREFIX)")

echo "${cmd[@]}" >> "$logfile"
if [[ $TEST_COMMANDS -eq 0 ]]; then
  "${cmd[@]}"
fi


genomes_files=(${OUT_PREFIX}.*snps.fa)

mutation_info_file=${OUT_PREFIX}.mutation_info.txt
genomes_mutation_info=$(for fasta_file in "${genomes_files[@]}"; do
  head -n 1 "$fasta_file" | sed 's/^>//g'
done)
echo "$genomes_mutation_info" > $mutation_info_file

awk 'NR==FNR {order[$1]=NR; next} {print order[$1], $0}' \
  $abundance_file \
  $mutation_info_file | \
  sort -n | \
  cut -d' ' -f2- > ${mutation_info_file}.tmp && mv ${mutation_info_file}.tmp $mutation_info_file



if [[ "$SEQTYPE" == "shotgun" ]]; then
  for fasta_file in "${genomes_files[@]}"; do
    echo $fasta_file
    sed -i 's/\//_/g' $fasta_file
    sed -i 's/|/_/g' $fasta_file
  done
  abundance_file_clean=${abundance_file}.clean
  sed 's/|/_/g' $abundance_file | sed 's/\//_/g' > $abundance_file_clean
  conda activate /private/home/bzhan146/miniconda3/envs/snakemake
  cmd=(iss generate --genomes "${genomes_files[@]}" \
    --abundance_file $abundance_file_clean \
    --n_reads $NUMREADS \
    --mode perfect \
    --seed $HASHED_SEED \
    --cpus $CPUS \
    -o ${OUT_PREFIX})
  echo "${cmd[@]}" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    "${cmd[@]}"
  fi
  conda deactivate
  rm $abundance_file_clean

elif [[ "$SEQTYPE" == "amplicon" ]]; then
  conda activate /private/home/bzhan146/miniconda3/envs/swampy
  merged_fasta=${OUT_PREFIX}.merged.fasta
  touch $merged_fasta
  swampy_temp_folder=${OUT_PREFIX}_swampy.tmp
  cat "${genomes_files[@]}" | awk '/^>/ {print $1; next} {print}' > $merged_fasta
  sed -i 's/|/_/g' $merged_fasta
  sed -i 's/\//_/g' $merged_fasta
  swampy_abundance_file=${OUT_PREFIX}.swampy.abundance.txt
  touch $swampy_abundance_file
  sed 's/|/_/g' $abundance_file | sed 's/\//_/g' > $swampy_abundance_file
  OUT_DIR=$(dirname "$OUT_PREFIX")
  PREFIX_BASE=$(basename "$OUT_PREFIX")
  cmd="python3 $SWAMPY \
          --genomes_file $merged_fasta \
          --genome_abundances $swampy_abundance_file \
          --temp_folder $swampy_temp_folder \
          --primer_set n2 \
          --output_folder $OUT_DIR \
          --output_filename_prefix $PREFIX_BASE \
          --n_reads $((NUMREADS / 2)) \
          --read_length 149 \
          --amplicon_pseudocounts 1000 \
          --autoremove \
          --no_pcr_errors"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    $cmd
  fi
  conda deactivate

  conda activate /private/home/bzhan146/miniconda3/envs/snakemake
  sorted_bam=${OUT_PREFIX}.sorted.bam
  trimmed_bam=${OUT_PREFIX}.trimmed.bam
  trimmed_sorted_bam=${OUT_PREFIX}.trimmed.sorted.bam

  cmd="bwa mem -t 4 ${REFERENCE_FASTA_FILE} ${OUT_PREFIX}_R1.fastq ${OUT_PREFIX}_R2.fastq | samtools sort | samtools view -h -F 4 -o ${sorted_bam}"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    eval "$cmd"
  fi
  cmd="samtools index --threads 4 ${sorted_bam}"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    $cmd
  fi
  cmd="ivar trim -e -q 1 -i ${sorted_bam} -b ${REFERENCE_PRIMER_BED_FILE} -p ${trimmed_bam} -x 3"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    $cmd
  fi
  cmd="samtools sort $trimmed_bam -o $trimmed_sorted_bam"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    $cmd
  fi
  conda deactivate

  conda activate /private/home/bzhan146/miniconda3/envs/jvarkit
  trimmed_sc_removed_bam=${OUT_PREFIX}.trimmed.sc.removed.bam
  cmd="java -jar $JVARKIT biostar84452 --samoutputformat BAM ${trimmed_sorted_bam} > ${trimmed_sc_removed_bam}"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    eval "$cmd"
  fi
  conda deactivate

  trimmed_fastq=${OUT_PREFIX}.trimmed.fastq
  cmd="samtools fastq --threads 4 ${trimmed_sc_removed_bam} -o ${trimmed_fastq}"
  echo "$cmd" >> "$logfile"
  if [[ $TEST_COMMANDS -eq 0 ]]; then
    $cmd
  fi

  rm ${OUT_PREFIX}_amplicon_abundances_summary.tsv
  rm ${OUT_PREFIX}.log
  rm ${OUT_PREFIX}_PCR_errors.vcf
  rm ${OUT_PREFIX}_R1.fastq
  rm ${OUT_PREFIX}_R2.fastq
  rm -r ${OUT_PREFIX}_swampy.tmp
  rm $merged_fasta
  rm $sorted_bam
  rm ${sorted_bam}.bai
  rm $swampy_abundance_file
  rm $trimmed_bam
  rm $trimmed_sorted_bam
  rm $trimmed_sc_removed_bam


fi


for fasta_file in "${genomes_files[@]}"; do
  rm $fasta_file
done
rm $random_nodes_file



# output files:
# shotgun: ${OUT_PREFIX}_R1.fastq, ${OUT_PREFIX}_R2.fastq ${OUT_PREFIX}.mutation_info.txt ${OUT_PREFIX}.abundance.txt
# amplicon: ${OUT_PREFIX}.trimmed.fastq ${OUT_PREFIX}.mutation_info.txt ${OUT_PREFIX}.abundance.txt
