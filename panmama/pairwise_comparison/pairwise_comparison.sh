#!/bin/bash

#SBATCH --job-name=pairwise-comparison
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/pairwise_comparison/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/pairwise_comparison/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-30

set -x
pairs_file=$1
data_dir=$2
output_dir=$3

mkdir -p "$output_dir"

mapfile -t pairs < <(cat $pairs_file)

total_pairs=${#pairs[@]}

echo "Total pairs: $total_pairs"

pairs_per_task=$((total_pairs / SLURM_ARRAY_TASK_COUNT))
remainder=$((total_pairs % SLURM_ARRAY_TASK_COUNT))

if [ $SLURM_ARRAY_TASK_ID -lt $remainder ]; then
  start_idx=$((SLURM_ARRAY_TASK_ID * (pairs_per_task + 1)))
  end_idx=$((start_idx + pairs_per_task + 1))
else
  start_idx=$((SLURM_ARRAY_TASK_ID * pairs_per_task + remainder))
  end_idx=$((start_idx + pairs_per_task))
fi


echo "Processing pairs from index $start_idx to $end_idx"

output_file="${output_dir}/comparisons_${SLURM_ARRAY_TASK_ID}.txt"

if [ -f "$output_file" ]; then
  rm "$output_file"
fi

for ((i=start_idx; i<end_idx; i++)); do
  IFS=$'\t' read -r seq1 seq2 <<< "${pairs[$i]}"
  seq1fasta_prefix=$(echo $seq1 | tr '\/\-\.\|' '____')
  seq2fasta_prefix=$(echo $seq2 | tr '\/\-\.\|' '____')
  fasta1path="${data_dir}/${seq1fasta_prefix}.unaligned.fasta"
  fasta2path="${data_dir}/${seq2fasta_prefix}.unaligned.fasta"
  if [ ! -f "$fasta1path" ] || [ ! -f "$fasta2path" ]; then
    echo "One of the fasta files does not exist: $fasta1path or $fasta2path"
    exit 1
  fi

  distance=$(bash calculate_errors_from_mafft.sh "$fasta1path" "$fasta2path")
  mapfile -t values < <(echo "$distance" | cut -f 1 -d ' ')
  errors="${values[0]}"
  snps="${values[1]}"
  snps_ambiguous="${values[2]}"
  gaps="${values[3]}"
  gap_edge_corrected="${values[4]}"
  echo -e "${seq1}\t${seq2}\t${snps}\t${snps_ambiguous}\t${gaps}\t${gap_edge_corrected}" >> "$output_file"

done