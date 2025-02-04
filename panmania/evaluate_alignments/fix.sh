#!
pairs_tsv=data/hiv/hiv20000_pairs.tsv
data_dir=data/hiv
output_file=out/hiv20000_alignment_differences.tsv

# Create temporary directory for parallel processing
tmp_dir=$(mktemp -d)
mkdir -p "$tmp_dir/results"

# Create a temporary file to store the pairs and paths
pairs_file="$tmp_dir/pairs.txt"
while IFS=' ' read -r parent child; do
    parent_cleaned=$(echo -n "$parent" | tr -c '[:alnum:]' '_')
    child_cleaned=$(echo -n "$child" |  tr -c '[:alnum:]' '_')
    echo -e "$parent $child $parent_cleaned $child_cleaned"
done < "$pairs_tsv"
