#!/bin/bash

getRank="$1"
READ_SAMPLES_DIR="$2"
TRUE_HAPLOTYPES_DIR="$3"

scores_files=$(find $READ_SAMPLES_DIR -name "*.testScores.txt")
for scores_file in $scores_files; do
  prefix=$(basename $scores_file .testScores.txt | cut -f 1,2 -d '_')
  rep=$(basename $scores_file .testScores.txt | cut -f 4 -d '_')
  output_prefix=$(basename $scores_file .testScores.txt)
  true_haplotypes_file=$TRUE_HAPLOTYPES_DIR/${prefix}_abundance_${rep}.tsv
  python3 $getRank $true_haplotypes_file $scores_file <(sort -k3,3 -gr $scores_file) | column -t > ${READ_SAMPLES_DIR}/${output_prefix}_rank_stats.tsv
done

