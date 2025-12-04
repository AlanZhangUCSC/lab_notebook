#!/bin/bash

#SBATCH --job-name=ancient-mito
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-651%50

set -x

mapfile -t fastq_paths < /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/final_analysis/fastq_paths.txt

fastq_path=${fastq_paths[$SLURM_ARRAY_TASK_ID]}

basename_fastq=$(basename ${fastq_path} .fq)

unique_identifier=$(uuidgen)
fastq_filtered="/data/tmp/${unique_identifier}.${basename_fastq}.filtered.fq"

# remove low complexity reads
~/tools/BBTools/bbduk.sh in=${fastq_path} out=${fastq_filtered} entropy=0.7 

# run panmap
PANMAP_PATH=/private/groups/corbettlab/alan/panmap
PANMAN_PATH=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/input_data/v_mtdna.panman
PMAI_PATH=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/input_data/v_mtdna.pmai
OUT_DIR=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/ancient_mito_out

docker load -i /private/groups/corbettlab/alan/panmap/panmap-dev.tar
docker run --rm \
  -v "$(realpath $PANMAP_PATH):/panmap" \
  -v "$(realpath $(dirname $PANMAN_PATH)):/panmans" \
  -v "$(realpath $(dirname $PMAI_PATH)):/pmais" \
  -v "$(realpath $fastq_filtered):/${basename_fastq}.fq" \
  -v "$(realpath $OUT_DIR):/output" \
  -w /panmap \
  --user "$(id -u):$(id -g)" \
  panmap-dev \
  bash -c "/panmap/build/bin/panmap \
              /panmans/$(basename $PANMAN_PATH) \
              /${basename_fastq}.fq \
              -m /pmais/$(basename $PMAI_PATH) \
              --read-scores \
              --prefix /output/${basename_fastq}.read_scores \
              --cpus 8"

docker run --rm \
  -v "$(realpath $PANMAP_PATH):/panmap" \
  -v "$(realpath $(dirname $PANMAN_PATH)):/panmans" \
  -v "$(realpath $(dirname $PMAI_PATH)):/pmais" \
  -v "$(realpath $fastq_filtered):/${basename_fastq}.fq" \
  -v "$(realpath $OUT_DIR):/output" \
  -w /panmap \
  --user "$(id -u):$(id -g)" \
  panmap-dev \
  bash -c "/panmap/build/bin/panmap \
              /panmans/$(basename $PANMAN_PATH) \
              /${basename_fastq}.fq \
              -m /pmais/$(basename $PMAI_PATH) \
              --overlap-coefficients 1000 \
              --prefix /output/${basename_fastq}.overlap_coefficients \
              --cpus 8"


rm ${fastq_filtered}
