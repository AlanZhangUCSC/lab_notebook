#!/bin/bash

#SBATCH --job-name=gen-data
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-59%10

set -x

OUT_DIR=$1
PANMAN_PATH=$2
PMAI_PATH=$3
OUT_DIR=$(realpath "$OUT_DIR")
PANMAN_PATH=$(realpath "$PANMAN_PATH")
PMAI_PATH=$(realpath "$PMAI_PATH")

mkdir -p $OUT_DIR

mapfile -t combinations < <(python3 /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gencomb.py \
  --snps 0 \
  --haplotypes 1 10 100 \
  --percent-mutated 1.0 \
  --seq-types shotgun amplicon \
  --num-reads 100000 1500000 \
  --num-rep 5 | tail -n +2)

read seqType numhap numsnps percentmutated numreads rep <<< "${combinations[$SLURM_ARRAY_TASK_ID]}"

docker load -i /private/groups/corbettlab/alan/panmap/panmap-dev.tar

prefix="${seqType}_${numhap}_${numsnps}_${percentmutated}_${numreads}_${rep}"
bash /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/genreads.sh \
  --seqtype $seqType \
  --numhap $numhap \
  --numsnp $numsnps \
  --permut $percentmutated \
  --numreads $numreads \
  --rep $rep \
  --cpus 8 \
  --panmap /private/groups/corbettlab/alan/panmap/ \
  --panman $PANMAN_PATH \
  --pmi $PMAI_PATH  \
  --random-seed $prefix \
  --swampy /private/home/bzhan146/tools/SWAMPy/src/simulate_metagenome.py \
  --reference-primer-bed-file /private/home/bzhan146/tools/SWAMPy/primer_sets/nimagenV2.bed \
  --reference-fasta-file /private/home/bzhan146/tools/SWAMPy/ref/MN908947.3.fasta \
  --jvarkit /private/home/bzhan146/tools/jvarkit/dist/jvarkit.jar \
  --out-prefix ${OUT_DIR}/${prefix}

