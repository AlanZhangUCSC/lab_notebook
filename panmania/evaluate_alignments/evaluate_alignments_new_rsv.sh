#!/bin/bash

#SBATCH --job-name=new_rsv_evaluate_alignments
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

pairs_tsv=data/rsv/rsv4000_pairs_new.tsv
data_dir=data/rsv
output_file=out/rsv4000_alignment_differences_new.tsv

# Create temporary directory for parallel processing
tmp_dir=$(mktemp -d)
mkdir -p "$tmp_dir/results"

# Create a temporary file to store the pairs and paths
pairs_file="$tmp_dir/pairs.txt"
while IFS=' ' read -r parent child; do
    parent_cleaned=$(echo -n "$parent" | tr -c '[:alnum:]' '_')
    child_cleaned=$(echo -n "$child" | tr -c '[:alnum:]' '_')
    echo -e "$parent $child $parent_cleaned $child_cleaned"
done < "$pairs_tsv" > "$pairs_file"

# Create the process_pair function
process_pair() {
    IFS=' ' read -r parent child parent_cleaned child_cleaned <<< "$1"
    
    # Check if variables are empty
    if [[ -z "$parent_cleaned" || -z "$child_cleaned" ]]; then
        echo "Error: Empty filename detected for pair $parent - $child" >&2
        return 1
    fi
    
    unaligned_parent_path="$data_dir/rsv4000_fasta_unaligned/$parent_cleaned.fasta"
    unaligned_child_path="$data_dir/rsv4000_fasta_unaligned/$child_cleaned.fasta"
    aligned_parent_path="$data_dir/rsv4000_fasta_aligned/$parent_cleaned.fasta"
    aligned_child_path="$data_dir/rsv4000_fasta_aligned/$child_cleaned.fasta"

    echo -e "Unaligned parent path: $unaligned_parent_path\nUnaligned child path: $unaligned_child_path\nAligned parent path: $aligned_parent_path\nAligned child path: $aligned_child_path" >&2
    
    # Check if files exist
    for file in "$unaligned_parent_path" "$unaligned_child_path" "$aligned_parent_path" "$aligned_child_path"; do
        if [[ ! -f "$file" ]]; then
            echo "Error: File not found: $file" >&2
            return 1
        fi
    done
    
    # Process PanMAN alignment
    panman_alignment_diff=$(python3 compare_aligned_sequences.py "$aligned_parent_path" "$aligned_child_path")
    if [[ $? -ne 0 ]]; then
        echo "Error: PanMAN comparison failed for $parent - $child" >&2
        return 1
    fi
    
    # Process MAFFT alignment
    tmp_joined_fasta=$(mktemp)
    tmp_aligned_fasta=$(mktemp)
    
    cat "$unaligned_parent_path" "$unaligned_child_path" > "$tmp_joined_fasta"
    mafft --auto --quiet "$tmp_joined_fasta" > "$tmp_aligned_fasta"
    mafft_alignment_diff=$(python3 compare_aligned_sequences.py "$tmp_aligned_fasta")
    
    # Cleanup temporary files
    rm "$tmp_joined_fasta" "$tmp_aligned_fasta"
    
    # Write results to temporary file
    echo -e "$parent\t$child\t$panman_alignment_diff\t$mafft_alignment_diff"
}

export -f process_pair
export data_dir

# Run the processing in parallel
echo -e "Parent\tChild\tPanMAN_Diff\tMAFFT_Diff" > "$output_file"
parallel --jobs 100 process_pair :::: "$pairs_file" >> "$output_file"

# Check if any errors occurred
if [[ $? -ne 0 ]]; then
    echo "Error: Some parallel processes failed. Check the output for details." >&2
    rm -rf "$tmp_dir"
    exit 1
fi

# Cleanup
rm -rf "$tmp_dir"