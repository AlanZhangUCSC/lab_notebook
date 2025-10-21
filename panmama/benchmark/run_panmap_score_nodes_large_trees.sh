#!/bin/bash

#SBATCH --job-name=gen-data
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-29%3

# task array 0-29

set -x

OUT_DIR=$1
PANMAP_PATH=$2
PANMAN_PATH=$3
PMAI_PATH=$4
DATA_DIR=$5
OUT_DIR=$(realpath "$OUT_DIR")
PANMAN_PATH=$(realpath "$PANMAN_PATH")
PMAI_PATH=$(realpath "$PMAI_PATH")
DATA_DIR=$(realpath "$DATA_DIR")

mkdir -p $OUT_DIR

mapfile -t combinations < <(python3 /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gencomb.py \
  --snps 0 \
  --haplotypes 1 10 100 \
  --percent-mutated 1.0 \
  --seq-types amplicon \
  --num-reads 100000 1500000 \
  --num-rep 5 | tail -n +2)

read seqType numhap numsnps percentmutated numreads rep <<< "${combinations[$SLURM_ARRAY_TASK_ID]}"

docker load -i /private/groups/corbettlab/alan/panmap/panmap-dev.tar

prefix="${seqType}_${numhap}_${numsnps}_${percentmutated}_${numreads}_${rep}"
output_prefix=$(basename "$PANMAN_PATH" .panman)_${prefix}

if [[ "$seqType" == "amplicon" ]]; then
  readpath="${prefix}.trimmed.fastq"
  docker run --rm \
    -v "$(realpath $PANMAP_PATH):/panmap" \
    -v "$(realpath $(dirname $PANMAN_PATH)):/panmans" \
    -v "$(realpath $(dirname $PMAI_PATH)):/pmais" \
    -v "$(realpath $DATA_DIR):/data" \
    -v "$(realpath $OUT_DIR):/output" \
    -w /panmap \
    --user "$(id -u):$(id -g)" \
    panmap-dev \
    bash -c "/panmap/build/bin/panmap \
                /panmans/$(basename $PANMAN_PATH) \
                /data/$readpath \
                -m /pmais/$(basename $PMAI_PATH) \
                --prefix /output/$output_prefix \
                --cpus 32"
elif [[ "$seqType" == "shotgun" ]]; then
  readpath1="${prefix}_R1.fastq"
  readpath2="${prefix}_R2.fastq"
  docker run --rm \
    -v "$(realpath $PANMAP_PATH):/panmap" \
    -v "$(realpath $(dirname $PANMAN_PATH)):/panmans" \
    -v "$(realpath $(dirname $PMAI_PATH)):/pmais" \
    -v "$(realpath $DATA_DIR):/data" \
    -v "$(realpath $OUT_DIR):/output" \
    -w /panmap \
    --user "$(id -u):$(id -g)" \
    panmap-dev \
    bash -c "/panmap/build/bin/panmap \
                /panmans/$(basename $PANMAN_PATH) \
                /data/$readpath1 \
                /data/$readpath2 \
                -m /pmais/$(basename $PMAI_PATH) \
                --prefix /output/$output_prefix \
                --cpus 32"
fi