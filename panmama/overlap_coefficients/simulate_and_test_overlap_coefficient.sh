#! /bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda deactivate

PANMAP=$1
PANMAN=$2
PMI=$3
SWAMPY=$4
REFERENCE_PRIMER_BED_FILE=$5
REFERENCE_FASTA_FILE=$6
JVARKIT=$7
OUT_DIR=$8
NUM_SNP=$9

NUM_SNPS=($NUM_SNP)
NUM_HAPS=(1 5 10 20 50 100)
SEQ_TYPES=(0 1) # 0: shot-gun 1: tiled amplicon  
DEPTH=(10 50 100 500 1000)

abundances_1hap=(1.0)
abundances_5haps=(0.5 0.2 0.15 0.1 0.05)
abundances_10haps=(0.15 0.15 0.15 0.15 0.1 0.1 0.05 0.05 0.05 0.05)
abundances_20haps=($(printf "0.05 %.0s" {1..20}))
abundances_50haps=($(printf "0.02 %.0s" {1..50}))
abundances_100haps=($(printf "0.01 %.0s" {1..100}))

parameter_combinations=()
for haps in "${NUM_HAPS[@]}"; do
  for snps in "${NUM_SNPS[@]}"; do
    for seq in "${SEQ_TYPES[@]}"; do
      for depth in "${DEPTH[@]}"; do
        parameter_combinations+=("$haps $snps $seq $depth")
      done
    done
  done
done


for parameter in "${parameter_combinations[@]}"; do
  read -r num_hap num_snp seq_type depth <<< "$parameter"

  echo "Processing: haps=$num_hap snps=$num_snp seq=$seq_type depth=$depth"

  output_prefix="${OUT_DIR}/${num_hap}_${num_snp}_${seq_type}_${depth}"
  random_seed_string="${num_hap}_${num_snp}_${seq_type}_${depth}"

  $PANMAP $PANMAN \
          --dump-random-nodeIDs $num_hap \
          --random-seed $random_seed_string \
          --prefix $output_prefix
  
  random_nodes_file=${output_prefix}.randomNodeIDs.txt
  dump_sequences_parameter_string=$(cat $random_nodes_file | tr '\n' ' ')
  simulate_snps_parameter_string=""
  for i in $(seq 0 $(($num_hap - 1))); do
    simulate_snps_parameter_string+=" $num_snp"
  done

  $PANMAP $PANMAN \
          --dump-sequences $dump_sequences_parameter_string \
          --simulate-snps $simulate_snps_parameter_string  \
          --random-seed $random_seed_string \
          --prefix $output_prefix

  abundance_file=${output_prefix}.abundance.txt

  if [[ "$num_hap" -eq 1 ]]; then
    paste "$random_nodes_file" <(printf '%s\n' "${abundances_1hap[@]}") > "$abundance_file"
  elif [[ "$num_hap" -eq 5 ]]; then
    paste "$random_nodes_file" <(printf '%s\n' "${abundances_5haps[@]}") > "$abundance_file"
  elif [[ "$num_hap" -eq 10 ]]; then
    paste "$random_nodes_file" <(printf '%s\n' "${abundances_10haps[@]}") > "$abundance_file"
  elif [[ "$num_hap" -eq 20 ]]; then
    paste "$random_nodes_file" <(printf '%s\n' "${abundances_20haps[@]}") > "$abundance_file"
  elif [[ "$num_hap" -eq 50 ]]; then
    paste "$random_nodes_file" <(printf '%s\n' "${abundances_50haps[@]}") > "$abundance_file"
  elif [[ "$num_hap" -eq 100 ]]; then
    paste "$random_nodes_file" <(printf '%s\n' "${abundances_100haps[@]}") > "$abundance_file"
  fi

  genomes_files=(${OUT_DIR}/${random_seed_string}.*.${num_snp}snps.fa)
  for file in "${genomes_files[@]}"; do
    echo "$file"
  done

  seed_hash=$(echo -n $random_seed_string | cksum | cut -d' ' -f1)
  if [[ "$seq_type" -eq 0 ]]; then
    conda activate /home/alan/miniforge3/envs/panmap
    iss generate --genomes "${genomes_files[@]}" \
      --abundance_file $abundance_file \
        --n_reads $((depth * 200)) \
        --model NovaSeq \
        --seed $seed_hash \
        -o ${output_prefix}
    conda deactivate

    $PANMAP $PANMAN \
            ${output_prefix}_R1.fastq \
            ${output_prefix}_R2.fastq \
            -m $PMI \
            --overlap-coefficients \
            --prefix $output_prefix
  elif [[ "$seq_type" -eq 1 ]]; then
    conda activate /home/alan/miniforge3/envs/SWAMPy
    merged_fasta=${output_prefix}.merged.fasta
    swampy_temp_folder=${output_prefix}_swampy.tmp
    cat "${genomes_files[@]}" | awk '/^>/ {print $1; next} {print}' > $merged_fasta
    sed -i 's/|/_/g' $merged_fasta
    sed -i 's/\//_/g' $merged_fasta
    swampy_abundance_file=${output_prefix}.swampy.abundance.txt
    sed 's/|/_/g' $abundance_file | sed 's/\//_/g' > $swampy_abundance_file
    python3 $SWAMPY \
            --genomes_file $merged_fasta \
            --genome_abundances $swampy_abundance_file \
            --temp_folder $swampy_temp_folder \
            --primer_set n2 \
            --output_folder $OUT_DIR \
            --output_filename_prefix $random_seed_string \
            --n_reads $((depth * 200)) \
            --read_length 149 \
            --amplicon_pseudocounts 1000 \
            --autoremove
    conda deactivate

    conda activate /home/alan/miniforge3/envs/panmap
    sorted_bam=${output_prefix}.sorted.bam
    trimmed_bam=${output_prefix}.trimmed.bam
    trimmed_sorted_bam=${output_prefix}.trimmed.sorted.bam
    bwa mem -t 4 ${REFERENCE_FASTA_FILE} ${output_prefix}_R1.fastq ${output_prefix}_R2.fastq | samtools sort | samtools view -h -F 4 -o ${sorted_bam}
    samtools index --threads 4 ${sorted_bam}
    ivar trim -e -q 1 -i ${sorted_bam} -b ${REFERENCE_PRIMER_BED_FILE} -p ${trimmed_bam}
    samtools sort $trimmed_bam -o $trimmed_sorted_bam
    conda deactivate

    conda activate /home/alan/miniforge3/envs/jvarkit-env
    trimmed_sc_removed_bam=${output_prefix}.trimmed.sc.removed.bam
    java -jar $JVARKIT biostar84452 --samoutputformat BAM ${trimmed_sorted_bam}  > ${trimmed_sc_removed_bam}
    conda deactivate

    trimmed_fastq=${output_prefix}.trimmed.fastq
    samtools fastq --threads 4 ${trimmed_sc_removed_bam} -o ${trimmed_fastq}
    $PANMAP $PANMAN \
            ${trimmed_fastq} \
            -m $PMI \
            --overlap-coefficients \
            --prefix $output_prefix
  fi


  node_rank_file=${output_prefix}.nodeRank.txt
  if [[ -f $node_rank_file ]]; then
    rm $node_rank_file
  fi
  for strain in $(cut -f 1 $abundance_file); do
    grep $strain -n ${output_prefix}.overlapCoefficients.txt >> ${node_rank_file}
  done

  sed -i 's/:/\t/g' ${node_rank_file}

  # clean up
  rm ${output_prefix}.randomNodeIDs.txt
  rm ${output_prefix}.*.${num_snp}snps.fa

  if [[ "$seq_type" -eq 0 ]]; then
    rm ${output_prefix}.iss.tmp*vcf
    rm ${output_prefix}_R1.fastq
    rm ${output_prefix}_R2.fastq
  fi

  if [[ "$seq_type" -eq 1 ]]; then
    rm ${output_prefix}.merged.fasta
    rm -rf $swampy_temp_folder
    rm ${output_prefix}_R1.fastq
    rm ${output_prefix}_R2.fastq
    rm $trimmed_fastq
    rm $trimmed_sorted_bam
    rm $trimmed_bam
    rm $sorted_bam
    rm ${sorted_bam}.bai
    rm $trimmed_sc_removed_bam
    rm $swampy_abundance_file
    rm ${output_prefix}_amplicon_abundances_summary.tsv
    rm ${output_prefix}_PCR_errors.vcf
    rm ${output_prefix}.log
  fi

  touch $node_rank_file
done