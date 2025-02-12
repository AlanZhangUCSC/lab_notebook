pairs_tsv=/private/groups/corbettlab/alan/panmania/panman/panmans/tb400_sequences/selected_pairs.tsv
data_dir=/private/groups/corbettlab/alan/panmania/panman/panmans/tb400_sequences/
output_file=out/tb400_alignment_exepcted_correction_new_panman.tsv

# Create temporary directory for parallel processing
tmp_dir=$(mktemp -d)
mkdir -p "$tmp_dir/results"

# Create a temporary file to store the pairs and paths
pairs_file="$tmp_dir/pairs.txt"
while IFS=$'\t' read -r child parent; do
    parent_cleaned=$(echo -n "$parent" | tr -c '[:alnum:]' '_')
    child_cleaned=$(echo -n "$child" | tr -c '[:alnum:]' '_')
    echo -e "$child $parent $parent_cleaned $child_cleaned"

    unaligned_parent_path="$data_dir/$child_cleaned.unaligned.fasta"
    unaligned_child_path="$data_dir/$parent_cleaned.unaligned.fasta"
    aligned_parent_path="$data_dir/$child_cleaned.aligned.fastags"
    aligned_child_path="$data_dir/$parent_cleaned.aligned.fastags"


    minimap_output=$(minimap2 -x asm20 "$unaligned_parent_path" "$unaligned_child_path")
    readarray -t minimap_values < <(echo -e "$minimap_output" | tr '\t' '\n')
    query_length=${minimap_values[1]}
    query_start=${minimap_values[2]}
    query_end=${minimap_values[3]}
    target_length=${minimap_values[6]}
    target_start=${minimap_values[7]}
    target_end=${minimap_values[8]}
    matching_bases_in_mapping=${minimap_values[9]}
    number_of_bases_and_gaps_in_mapping=${minimap_values[10]}

    echo -e "$query_length\t$query_start\t$query_end\t$target_length\t$target_start\t$target_end\t$matching_bases_in_mapping\t$number_of_bases_and_gaps_in_mapping"
    
    minimap_diff=$((query_length - query_end + query_start + target_length - target_end + target_start + number_of_bases_and_gaps_in_mapping - matching_bases_in_mapping))

    echo -e "$parent\t$child\t$minimap_diff"
    break
done < "$pairs_tsv"

# Cleanup
rm -rf "$tmp_dir"