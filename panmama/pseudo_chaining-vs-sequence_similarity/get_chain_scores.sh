#!/bin/bash

#SBATCH --job-name=rsv_panmama_with_mutations
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=25gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/pseudo_chaining-vs-sequence_similarity/logs/%x.%j.%A.%a.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/pseudo_chaining-vs-sequence_similarity/logs/%x.%j.%A.%a.err
#SBATCH --partition=short
#SBATCH --time=01:00:00

tree_name=rsv4000
panman_path=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/rsv/k19_s8_null/panmap_outputs/panmap2_index/${tree_name}.panman

panmap \
  ${panman_path} \
  /private/groups/corbettlab/alan/lab_notebook/panmama/pseudo_chaining-vs-sequence_similarity/simulated_reads.fastq \
  -k 19 -s 8 \
  --place-per-read \
  --prefix simulated_reads \
  --check-frequency 20 \
  --remove-iteration 20 \
  --redo-read-threshold 0 \
  --em-filter-round 2 \
  --remove-threshold 0.01 \
  --rounds-remove 5 \
  --preem-filter-method mbc \
  --save-kminmer-binary-coverage \
  --call-subconsensus \
  --cpus 1
