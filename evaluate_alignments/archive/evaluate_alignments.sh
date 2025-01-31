#!/bin/bash

pairs_tsv=$1
data_dir=$2
output_file=$3

parents=($(head -n 20 $pairs_tsv | cut -f1 -d ' '))
children=($(head -n 20 $pairs_tsv | cut -f2 -d ' '))
panman_alignment_diffs=($(for i in $parents; do echo 0; done))
mafft_alignment_diffs=($(for i in $parents; do echo 0; done))

for i in $(seq 0 $((${#parents[@]} - 1))); do
    parent=${parents[$i]}
    child=${children[$i]}
    parent_cleaned=$(echo $parent | tr -cd '[:alnum:]_.' | tr '.' '_')
    child_cleaned=$(echo $child | tr -cd '[:alnum:]_.' | tr '.' '_')
    unaligned_parent_path=$data_dir/fasta_unaligned/$parent_cleaned.fasta
    unaligned_child_path=$data_dir/fasta_unaligned/$child_cleaned.fasta
    aligned_parent_path=$data_dir/fasta_aligned/$parent_cleaned.fasta
    aligned_child_path=$data_dir/fasta_aligned/$child_cleaned.fasta
    panman_alignment_diff=$(python3 compare_aligned_sequences.py $aligned_parent_path $aligned_child_path)
    tmp_joined_fasta=$(mktemp)
    tmp_aligned_fasta=$(mktemp)
    cat $unaligned_parent_path $unaligned_child_path > $tmp_joined_fasta
    mafft --auto --quiet $tmp_joined_fasta > $tmp_aligned_fasta
    mafft_alignment_diff=$(python3 compare_aligned_sequences.py $tmp_aligned_fasta)
    rm $tmp_aligned_fasta
    rm $tmp_joined_fasta
    panman_alignment_diffs[$i]=$panman_alignment_diff
    mafft_alignment_diffs[$i]=$mafft_alignment_diff
done

echo -e "Parent\tChild\tPanMAN_Diff\tMAFFT_Diff" > $output_file
for i in $(seq 0 $((${#parents[@]} - 1))); do
    parent=${parents[$i]}
    child=${children[$i]}
    panman_diff=${panman_alignment_diffs[$i]}
    mafft_diff=${mafft_alignment_diffs[$i]}
    echo -e "$parent\t$child\t$panman_diff\t$mafft_diff" >> $output_file
done