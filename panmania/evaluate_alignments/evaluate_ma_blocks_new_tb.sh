#!/bin/bash

#SBATCH --job-name=tb_ma_blocks_correction_new_panman
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmania/evaluate_alignments/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmania/evaluate_alignments/logs/%x.%A.%a.%j.err
#SBATCH --partition=medium
#SBATCH --time=04:00:00

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
done < "$pairs_tsv" > "$pairs_file"

# Create the process_pair function
process_pair() {
    IFS=' ' read -r child parent parent_cleaned child_cleaned <<< "$1"
    
    # Check if variables are empty
    if [[ -z "$parent_cleaned" || -z "$child_cleaned" ]]; then
        echo "Error: Empty filename detected for pair $child - $parent" >&2
        return 1
    fi
    
    unaligned_parent_path="$data_dir/$child_cleaned.unaligned.fasta"
    unaligned_child_path="$data_dir/$parent_cleaned.unaligned.fasta"
    aligned_parent_path="$data_dir/$child_cleaned.aligned.fastags"
    aligned_child_path="$data_dir/$parent_cleaned.aligned.fastags"
    
    # Check if files exist
    for file in "$unaligned_parent_path" "$unaligned_child_path" "$aligned_parent_path" "$aligned_child_path"; do
        if [[ ! -f "$file" ]]; then
            echo "Error: File not found: $file" >&2
            return 1
        fi
    done
    
    # Process PanMAN alignment
    output=$(python3 search_ma_blocks_new_panman.py "$aligned_parent_path" "$aligned_child_path")

    readarray -t values < <(echo -e "$output" | cut -d ':' -f 2)
    all_mismatch_count=${values[0]}
    expected_mismatch_count_after_correction=${values[1]}
    mismatch_count_in_ma_blocks=${values[2]}

    # minimap alignment
    minimap_output=$(minimap2 -x asm20 -w 5 "$unaligned_parent_path" "$unaligned_child_path")
    readarray -t minimap_values < <(echo -e "$minimap_output" | tr '\t' '\n')
    query_length=${minimap_values[1]}
    query_start=${minimap_values[2]}
    query_end=${minimap_values[3]}
    target_length=${minimap_values[6]}
    target_start=${minimap_values[7]}
    target_end=${minimap_values[8]}
    matching_bases_in_mapping=${minimap_values[9]}
    number_of_bases_and_gaps_in_mapping=${minimap_values[10]}

    minimap_diff=$((query_length - query_end + query_start + target_length - target_end + target_start + number_of_bases_and_gaps_in_mapping - matching_bases_in_mapping))
    
    # Write results to temporary file
    echo -e "$parent\t$child\t$all_mismatch_count\t$minimap_diff\t$expected_mismatch_count_after_correction\t$mismatch_count_in_ma_blocks"
}

export -f process_pair
export data_dir

# Run the processing in parallel
echo -e "Parent\tChild\tPanMAN_Diff\tMinimap_Diff\tExpected_Mismatch_Count_After_Correction\tMismatch_Count_In_MA_Blocks" > "$output_file"
parallel --jobs 100 process_pair :::: "$pairs_file" >> "$output_file"

# Check if any errors occurred
if [[ $? -ne 0 ]]; then
    echo "Error: Some parallel processes failed. Check the output for details." >&2
    rm -rf "$tmp_dir"
    exit 1
fi

# Cleanup
rm -rf "$tmp_dir"