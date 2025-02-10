#!/bin/bash

#SBATCH --job-name=sars_evaluate_consensus_heuristic_7
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/sars/with_mutations_1000x_top2_score/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/sars/with_mutations_1000x_top2_score/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-149%50


num_replicates=10
mixed_haplotypes=(1 2 3 5 10)
num_snps=(1 5 10)

# Create a list of all parameter combinations
declare -a param_combinations
index=0
for ((replicate=1; replicate<=num_replicates; replicate++)); do
    for haplotype in "${mixed_haplotypes[@]}"; do
        for snp in "${num_snps[@]}"; do
            param_combinations[$index]="($replicate $haplotype $snp)"
            ((index++))
        done
    done
done

test_combination=${param_combinations[$SLURM_ARRAY_TASK_ID]}
read rep hap snp <<< $(echo $test_combination | tr -d '()')
tree_name_short=sars
tree_name=""
if [ "$tree_name_short" == "hiv" ]; then
    reads_per_coverage=60
    tree_name=hiv20000
elif [ "$tree_name_short" == "rsv" ]; then
    reads_per_coverage=100
    tree_name=rsv4000
elif [ "$tree_name_short" == "sars" ]; then
    reads_per_coverage=200
    tree_name=sars20000
else
    echo "Error: Unknown tree name $tree_name"
    exit 1
fi

echo "Evaluating ${tree_name}: replicate $rep with $hap haplotypes and $snp snps"
prefix=${tree_name}_rep${rep}_${hap}hap_${snp}snp
data_dir=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/${tree_name_short}/with_mutations_1000x_top2_score/rep${rep}
ref_dir=/private/groups/corbettlab/alan/panmap/dev/panmama_eval/${tree_name}/fasta/info/fasta_unaligned
panmap_output_dir=${data_dir}/panmap_outputs
abundance_file=${panmap_output_dir}/${prefix}_abundance.tsv
assigned_reference_file=${panmap_output_dir}/${prefix}_assigned_references.tsv
assigned_reads_dir=${panmap_output_dir}/${prefix}_reads
consensus_output_dir=${panmap_output_dir}/${prefix}_consensus_heuristic_merged
command_text=${consensus_output_dir}/heuristic7_command_text.txt

if [ -f "$command_text" ]; then
  rm "$command_text"
fi
mkdir -p $consensus_output_dir

merged_abundance_file=${panmap_output_dir}/${prefix}_merged_abundance.tsv
if [ ! -f "$merged_abundance_file" ]; then
  python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/merge_abundance.py $abundance_file /private/groups/corbettlab/alan/panmap/dev/panmama_eval/${tree_name}/fasta/info/fasta_unaligned > $merged_abundance_file
fi


estimated_haplotypes=($(cut -f 1 $merged_abundance_file | awk -F',' '{print $1}'))
estimated_abundance=($(cut -f 2 $merged_abundance_file))



for ((i=0; i<${#estimated_haplotypes[@]}; i++)); do
  current_estimated_haplotype=${estimated_haplotypes[$i]}
  current_estimated_haplotype_sanitize_for_ref=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_' | tr '-' '_' | tr '.' '_')
  ref_fasta=${ref_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta
  cp $ref_fasta $consensus_output_dir
  ref_fasta=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta
done

aligned_fasta=${consensus_output_dir}/aligned.fasta
cat "${consensus_output_dir}"/*.fasta | mafft --quiet - > $aligned_fasta

# calling consensus for each estimated haplotype
pids=()
for ((i=0; i<${#estimated_haplotypes[@]}; i++)); do
  current_estimated_haplotype=${estimated_haplotypes[$i]}
  current_estimated_abundance=${estimated_abundance[$i]}
  echo "processing estimated haplotype ${current_estimated_haplotype} with abundance ${current_estimated_abundance}"
  
  current_estimated_haplotype_sanitize_for_reads=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_')
  read1=${assigned_reads_dir}/${current_estimated_haplotype_sanitize_for_reads}_R1.fastq.gz
  read2=${assigned_reads_dir}/${current_estimated_haplotype_sanitize_for_reads}_R2.fastq.gz
  
  current_estimated_haplotype_sanitize_for_ref=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_' | tr '-' '_' | tr '.' '_')
  ref_fasta=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta

  consensus_output_prefix=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_reads}
  # if [ -f "${consensus_output_prefix}.varonly.vcf.gz" ]; then
  #   continue
  # fi
  (
    # Generate minimap2 index (optional, but can speed up multiple alignments)
    minimap2 -d "${ref_fasta}.mmi" "$ref_fasta"

    # Align reads using minimap2 and pipe to samtools for sorting
    minimap2 -ax sr "${ref_fasta}" "$read1" "$read2" | \
        samtools sort -o "${consensus_output_prefix}.sorted.bam"
    
    samtools index "${consensus_output_prefix}.sorted.bam"
    
    bcftools mpileup -f "$ref_fasta" "${consensus_output_prefix}.sorted.bam" \
        --max-depth 10000 --max-idepth 10000 -B \
        -o "${consensus_output_prefix}.mpileup"
    
    cat "${consensus_output_prefix}.mpileup" | \
        bcftools call -m \
        --ploidy 1 \
        --keep-alts \
        -Oz -o "${consensus_output_prefix}.all.vcf.gz"
    
    cat "${consensus_output_prefix}.mpileup" | \
        bcftools call -vm \
        --ploidy 1 \
        --keep-alts \
        -Oz -o "${consensus_output_prefix}.varonly.vcf.gz"

    bcftools index "${consensus_output_prefix}.all.vcf.gz"
    bcftools index "${consensus_output_prefix}.varonly.vcf.gz"

    python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/subconsensus/gen_pcf.py \
      --vcf "${consensus_output_prefix}.all.vcf.gz" \
      --bam "${consensus_output_prefix}.sorted.bam" \
      --abundances "${merged_abundance_file}" \
      --assignments "${assigned_reference_file}" \
      --fasta "$aligned_fasta" \
      --min_coverage 10 \
      --min_allele_depth 5 \
      --error_rate 0.05 \
      > "${consensus_output_prefix}.pcf"

    rm "${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta.amb"
    rm "${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta.ann"
    rm "${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta.bwt"
    rm "${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta.pac"
    rm "${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta.sa"
    rm "${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta.fai"
    rm "${consensus_output_prefix}.varonly.vcf.gz.csi"
    rm "${consensus_output_prefix}.all.vcf.gz.csi"
  ) &
  pids+=($!)

  # Limit number of parallel processes
  if [ ${#pids[@]} -ge 4 ]; then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
  fi
done

# Wait for remaining processes to complete
for pid in "${pids[@]}"; do
  wait $pid
done

all_vcf_gzs=""
all_sorted_bams=""
all_pcfs=""
all_output=""
for ((i=0; i<${#estimated_haplotypes[@]}; i++)); do
    current_estimated_haplotype=${estimated_haplotypes[$i]}
    current_estimated_haplotype_sanitize_for_reads=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_')
    current_estimated_haplotype_sanitize_for_ref=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_' | tr '-' '_' | tr '.' '_')
    consensus_output_prefix=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_reads}
    ref_fasta=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta


    all_vcf_gzs="${all_vcf_gzs} ${consensus_output_prefix}.all.vcf.gz"
    all_sorted_bams="${all_sorted_bams} ${consensus_output_prefix}.sorted.bam"
    all_pcfs="${all_pcfs} ${consensus_output_prefix}.pcf"
    all_output="${all_output} ${consensus_output_prefix}.heuristic_kl_gap_merged_heuristic7.varonly.vcf"

done

cmd="python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/subconsensus/heuristic_v7.py \
  --vcfs ${all_vcf_gzs} \
  --bams ${all_sorted_bams} \
  --pcfs ${all_pcfs} \
  --abundances ${merged_abundance_file} \
  --assignments ${assigned_reference_file} \
  --fasta $aligned_fasta \
  --confidence_threshold 2 \
  --outputs ${all_output}"

echo "$cmd"
echo "$cmd" >> $command_text
eval "$cmd" 

for ((i=0; i<${#estimated_haplotypes[@]}; i++)); do
  current_estimated_haplotype=${estimated_haplotypes[$i]}
  current_estimated_haplotype_sanitize_for_reads=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_')
  current_estimated_haplotype_sanitize_for_ref=$(echo $current_estimated_haplotype | tr '/' '_' | tr '|' '_' | tr '-' '_' | tr '.' '_')
  consensus_output_prefix=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_reads}
  ref_fasta=${consensus_output_dir}/${current_estimated_haplotype_sanitize_for_ref}.fasta

  rm "${consensus_output_prefix}.heuristic_kl_gap_merged_heuristic7.varonly.vcf.gz"
  rm "${consensus_output_prefix}.heuristic_kl_gap_merged_heuristic7.varonly.vcf.gz.csi"
  bgzip "${consensus_output_prefix}.heuristic_kl_gap_merged_heuristic7.varonly.vcf"

  bcftools index "${consensus_output_prefix}.heuristic_kl_gap_merged_heuristic7.varonly.vcf.gz"

  if [ ! -f "${ref_fasta}.fai" ]; then
      samtools faidx "$ref_fasta"
  fi

  bcftools consensus -f "$ref_fasta" \
      -H 1 \
      "${consensus_output_prefix}.heuristic_kl_gap_merged_heuristic7.varonly.vcf.gz" \
      -o "${consensus_output_prefix}.consensus.fa"
done


# Pairwise comparison of consensus sequences
consensus_files=(${consensus_output_dir}/*.consensus.fa)
unique_groups=()
processed=()

for ((i=0; i<${#consensus_files[@]}; i++)); do
    if [[ " ${processed[@]} " =~ " ${consensus_files[i]} " ]]; then
        continue
    fi
    
    current_group=("${consensus_files[i]}")
    
    for ((j=i+1; j<${#consensus_files[@]}; j++)); do
        if [[ " ${processed[@]} " =~ " ${consensus_files[j]} " ]]; then
            continue
        fi
        
        # Create unique temporary filename using PID and timestamp
        tmp_aligned="tmp_aligned_$$_$(date +%s%N).fasta"
        
        # Compare sequences using MAFFT
        cat "${consensus_files[i]}" "${consensus_files[j]}" | mafft --quiet - > "$tmp_aligned"
        
        # Get comparison results
        comparison_output=$(python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/calculate_errors_from_mafft.py "$tmp_aligned" 2>/dev/null)
        rm "$tmp_aligned"
        
        # Parse the last 4 lines for stats
        num_snps=$(echo "$comparison_output" | tail -n 3 | head -n 1 | awk '{print $1}')
        gaps_edge_corrected=$(echo "$comparison_output" | tail -n 1 | awk '{print $1}')
        
        # If sequences are identical (sum of SNPs and edge-corrected gaps is 0)
        if [ $((num_snps + gaps_edge_corrected)) -eq 0 ]; then
            current_group+=("${consensus_files[j]}")
            processed+=("${consensus_files[j]}")
        fi
    done
    
    if [ ${#current_group[@]} -gt 1 ]; then
        unique_groups+=("${current_group[*]}")
    fi
    processed+=("${consensus_files[i]}")
done

# Print results (for information)
echo -e "\nGroups of identical sequences:"
for group in "${unique_groups[@]}"; do
    echo -e "\nGroup:"
    for member in $group; do
        echo "  - $(basename "$member")"
    done
done

# Create final list of unique sequences
echo -e "\nFinal list of unique consensus sequences:"
unique_consensus_list=()

# Add one representative from each group (using the first sequence in each group)
for group in "${unique_groups[@]}"; do
    representative=$(echo "$group" | awk '{print $1}')
    unique_consensus_list+=("$representative")
done

# Add sequences that aren't in any group
for file in "${consensus_files[@]}"; do
    # Check if the file is in any group
    is_grouped=false
    for group in "${unique_groups[@]}"; do
        if [[ "$group" == *"$file"* ]]; then
            is_grouped=true
            break
        fi
    done
    
    # If not in any group, add it to the unique list
    if [ "$is_grouped" = false ]; then
        unique_consensus_list+=("$file")
    fi
done

# Print and save the final list
echo "Unique sequences (including group representatives):"
for seq in "${unique_consensus_list[@]}"; do
    echo "$(basename "$seq")"
done

# Save the final list to a file
output_list="${consensus_output_dir}/unique_consensus_list.txt"
printf "%s\n" "${unique_consensus_list[@]}" > "$output_list"
echo -e "\nFinal list saved to: $output_list"




# compare to true sequences
true_sequences_file=${data_dir}/var_fasta/${prefix}.fa
comparison_output="${consensus_output_dir}/consensus_to_true_comparison_heuristic_kl_gap_merged_heuristic7.tsv"
echo -e "true_haplotype\tconsensus\tclosest_true\tsnp_dist\tgap_edge_corrected_dist" > "$comparison_output"

# Create temporary directory for intermediate files
tmp_dir="${consensus_output_dir}/tmp_$$"
mkdir -p "$tmp_dir"

# Split true sequences into separate files
current_file=""
while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        current_header=${line#>}
        current_header=$(echo "$current_header" | tr '/' '_' | tr '|' '_' | tr '-' '_' | tr '.' '_')
        current_file="${tmp_dir}/true_${current_header}.fa"
        echo "$line" > "$current_file"
    else
        echo "$line" >> "$current_file"
    fi
done < "$true_sequences_file"

# Get list of true sequence files
true_sequence_files=("${tmp_dir}"/true_*.fa)

# Track which consensus sequences are assigned to true haplotypes
declare -A assigned_consensus
for consensus in "${unique_consensus_list[@]}"; do
    assigned_consensus["$(basename "$consensus")"]="false"
done

# First, find closest consensus for each true sequence
for true_file in "${true_sequence_files[@]}"; do
    true_seq=$(basename "$true_file" .fa | sed 's/^true_//')
    
    # Variables to track closest match
    closest_consensus=""
    min_distance=999999
    closest_snps=0
    closest_gaps=0
    closest_snp_positions=""
    closest_gap_positions=""
    closest_left_nongap_start=""
    closest_right_nongap_end=""
    
    # Compare with each unique consensus
    for consensus in "${unique_consensus_list[@]}"; do
        tmp_aligned="${tmp_dir}/aligned_$$_$(date +%s%N).fasta"
        
        cat "$true_file" "$consensus" | mafft --quiet - > "$tmp_aligned"
        
        mafft_distance=$(python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/calculate_errors_from_mafft.py "$tmp_aligned" 2>/dev/null)
        rm "$tmp_aligned"
        
        num_snps=$(echo "$mafft_distance" | tail -n 3 | head -n 1 | awk '{print $1}')
        snp_positions=$(echo "$mafft_distance" | tail -n 3 | head -n 1 | awk '{print $3}')
        gaps_edge_corrected=$(echo "$mafft_distance" | tail -n 1 | awk '{print $1}')
        gap_positions=$(echo "$mafft_distance" | tail -n 2 | head -n 1 | awk '{print $5}')
        left_nongap_start=$(echo "$mafft_distance" | tail -n 2 | head -n 1 | awk '{print $3}')
        right_nongap_end=$(echo "$mafft_distance" | tail -n 2 | head -n 1 | awk '{print $4}')

        total_distance=$((num_snps + gaps_edge_corrected))
        
        if [ "$total_distance" -lt "$min_distance" ]; then
            min_distance=$total_distance
            closest_consensus=$(basename "$consensus")
            closest_snps=$num_snps
            closest_gaps=$gaps_edge_corrected
            closest_snp_positions=$snp_positions
            closest_gap_positions=$gap_positions
            closest_left_nongap_start=$left_nongap_start
            closest_right_nongap_end=$right_nongap_end
        fi
    done
    
    # Write result for true haplotype
    echo -e "${true_seq}\t${closest_consensus}\t.\t${closest_snps}\t${closest_gaps}\t${closest_left_nongap_start}\t${closest_right_nongap_end}\t${closest_snp_positions}\t${closest_gap_positions}" >> "$comparison_output"
    assigned_consensus["$closest_consensus"]="true"
done

# Now process unassigned consensus sequences
for consensus in "${unique_consensus_list[@]}"; do
    consensus_base=$(basename "$consensus")
    if [ "${assigned_consensus[$consensus_base]}" = "false" ]; then
        # Find closest true sequence for this unassigned consensus
        min_distance=999999
        closest_true=""
        closest_snps=0
        closest_gaps=0
        closest_snp_positions=""
        closest_gap_positions=""
        closest_left_nongap_start=""
        closest_right_nongap_end=""
        
        for true_file in "${true_sequence_files[@]}"; do
            tmp_aligned="${tmp_dir}/aligned_$$_$(date +%s%N).fasta"
            
            cat "$consensus" "$true_file" | mafft --quiet - > "$tmp_aligned"
            
            mafft_distance=$(python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/calculate_errors_from_mafft.py "$tmp_aligned" 2>/dev/null)
            rm "$tmp_aligned"
            
            num_snps=$(echo "$mafft_distance" | tail -n 3 | head -n 1 | awk '{print $1}')
            snp_positions=$(echo "$mafft_distance" | tail -n 3 | head -n 1 | awk '{print $3}')
            gaps_edge_corrected=$(echo "$mafft_distance" | tail -n 1 | awk '{print $1}')
            gap_positions=$(echo "$mafft_distance" | tail -n 2 | head -n 1 | awk '{print $5}')
            left_nongap_start=$(echo "$mafft_distance" | tail -n 2 | head -n 1 | awk '{print $3}')
            right_nongap_end=$(echo "$mafft_distance" | tail -n 2 | head -n 1 | awk '{print $4}')
            total_distance=$((num_snps + gaps_edge_corrected))
            
            if [ "$total_distance" -lt "$min_distance" ]; then
                min_distance=$total_distance
                closest_true=$(basename "$true_file" .fa | sed 's/^true_//')
                closest_snps=$num_snps
                closest_gaps=$gaps_edge_corrected
                closest_snp_positions=$snp_positions
                closest_gap_positions=$gap_positions
                closest_left_nongap_start=$left_nongap_start
                closest_right_nongap_end=$right_nongap_end
            fi
        done
        
        # Write result for unassigned consensus
        echo -e ".\t${consensus_base}\t${closest_true}\t${closest_snps}\t${closest_gaps}\t${closest_left_nongap_start}\t${closest_right_nongap_end}\t${closest_snp_positions}\t${closest_gap_positions}" >> "$comparison_output"
    fi
done

# Clean up
rm -rf "$tmp_dir"

echo "Comparison results saved to: $comparison_output"