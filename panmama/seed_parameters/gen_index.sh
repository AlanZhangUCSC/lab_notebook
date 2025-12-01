#!/bin/bash

#SBATCH --job-name=gen-index
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/seed_parameters/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/seed_parameters/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-29

ls=(1 3 5)
ks=(10 15 20 25 30 35 40 45 50)
ss=(5 7 11 13 17 19 23 29 31 37 41 43 47)
ts=(1 3 5 7)

PANMAP_PATH=$1
PANMAN_PATH=$2
PMAI_PREFIX=$3
OUTPUT_DIR=$4
docker run --rm \
  -v "$(realpath $PANMAP_PATH):/panmap" \
  -v "$(realpath $(dirname $PANMAN_PATH)):/panmans" \
  -v "$(realpath $OUTPUT_DIR:/output_dir" \
  -w /panmap \
  --user "$(id -u):$(id -g)" \
  panmap-dev \
  bash -c "/panmap/build/bin/panmap \
              /panmans/$(basename $PANMAN_PATH) \
              --index-mgsr /output_dir/${PMAI_PREFIX}.pmai"