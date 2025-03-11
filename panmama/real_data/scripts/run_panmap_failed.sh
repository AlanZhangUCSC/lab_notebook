#!/bin/bash

#SBATCH --job-name=rsv_real_data
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/logs/%x.%A.%a.%j.err
#SBATCH --partition=medium
#SBATCH --time=02:00:00
#SBATCH --array=4,5,10,11,14,17,21,24,27,29,41,44,47,62,81,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191%30

rsv_samples=(RSV00059 RSV00060 RSV00061 RSV00062 RSV00063 RSV00064 RSV00065 RSV00066 RSV00067 RSV00068 RSV00069 RSV00070 RSV00071 RSV00072 RSV00073 RSV00074 RSV00075 RSV00076 RSV00077 RSV00078 RSV00079 RSV00080 RSV00081 RSV00082 RSV00083 RSV00084 RSV00085 RSV00086 RSV00087 RSV00088 RSV00089 RSV00090 RSV00091 RSV00092 RSV00093 RSV00094 RSV00095 RSV00096 RSV00097 RSV00098 RSV00099 RSV00100 RSV00101 RSV00102 RSV00103 RSV00104 RSV00105 RSV00106 RSV00107 RSV00108 RSV00109 RSV00111 RSV00112 RSV00113 RSV00114 RSV00115 RSV00116 RSV00117 RSV00118 RSV00119 RSV00120 RSV00121 RSV00122 RSV00123 RSV00124 RSV00125 RSV00126 RSV00127 RSV00128 RSV00129 RSV00130 RSV00131 RSV00132 RSV00133 RSV00134 RSV00135 RSV00136 RSV00137 RSV00138 RSV00139 RSV00140 RSV00141 RSV00142 RSV00143 RSV00144 RSV00145 RSV00146 RSV00147 RSV00148 RSV00149 RSV00150 RSV00151 RSV00152 RSV00153 RSV00154 RSV00159 RSV00164 RSV00165 RSV00166 RSV00167 RSV00168 RSV00169 RSV00170 RSV00171 RSV00172 RSV00173 RSV00174 RSV00175 RSV00176 RSV00177 RSV00178 RSV00179 RSV00180 RSV00181 RSV00182 RSV00183 RSV00184 RSV00185 RSV00186 RSV00187 RSV00188 RSV00189 RSV00190 RSV00191 RSV00192 RSV00193 RSV00194 RSV00195 RSV00196 RSV00197 RSV00198 RSV00199 RSV00200 RSV00201 RSV00202 RSV00203 RSV00204 RSV00205 RSV00206 RSV00207 RSV00208 RSV00209 RSV00210 RSV00211 RSV00212 RSV00213 RSV00214 RSV00215 RSV00216 RSV00217 RSV00218 RSV00219 RSV00220 RSV00221 RSV00222 RSV00223 RSV00224 RSV00225 RSV00226 RSV00227 RSV00228 RSV00229 RSV00230 RSV00231 RSV00232 RSV00233 RSV00234 RSV00235 RSV00236 RSV00237 RSV00238 RSV00239 RSV00240 RSV00241 RSV00242 RSV00243 RSV00244 RSV00245 RSV00246 RSV00247 RSV00248 RSV00249 RSV00250 RSV00251 RSV00252 RSV00253 RSV00254 RSV00255 RSV00256 RSV00257 RSV00258 RSV00259)

data_dir=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data
RSV_A_forward_primer_file=$data_dir/other_files_you_might_need/RSVA_forward.primers.fa
RSV_A_reverse_primer_file=$data_dir/other_files_you_might_need/RSVA_reverse.primers.fa
RSV_B_forward_primer_file=$data_dir/other_files_you_might_need/RSVB_forward.primers.fa
RSV_B_reverse_primer_file=$data_dir/other_files_you_might_need/RSVB_reverse.primers.fa
RSV_A_ref_fasta=$data_dir/other_files_you_might_need/RS20000581_RSVA_reference.fasta
RSV_B_ref_fasta=$data_dir/other_files_you_might_need/RE20000104_RSVB_reference.fasta
RSV_fasta_dir=/private/groups/corbettlab/alan/panmap/dev/panmama_eval/rsv4000/fasta/info/fasta_unaligned

panman_file=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/rsv/k19_s8_null/panmap_outputs/panmap2_index/rsv4000.panman

rsv_dir=$data_dir/${rsv_samples[$SLURM_ARRAY_TASK_ID]}

for type in A B
do
  R1_file=$(ls ${rsv_dir}/${type}/*RSV*R1*fastq*)
  R2_file=$(ls ${rsv_dir}/${type}/*RSV*R2*fastq*)
  R1_trimmed_strict=${rsv_dir}/${type}/trimmed_strict_R1.fastq
  R2_trimmed_strict=${rsv_dir}/${type}/trimmed_strict_R2.fastq
  R1_trimmed_loose=${rsv_dir}/${type}/trimmed_loose_R1.fastq
  R2_trimmed_loose=${rsv_dir}/${type}/trimmed_loose_R2.fastq
  R1_untrimmed_strict=${rsv_dir}/${type}/untrimmed_strict_R1.fastq
  R2_untrimmed_strict=${rsv_dir}/${type}/untrimmed_strict_R2.fastq
  R1_untrimmed_loose=${rsv_dir}/${type}/untrimmed_loose_R1.fastq
  R2_untrimmed_loose=${rsv_dir}/${type}/untrimmed_loose_R2.fastq
  trimmed_strict_type_result=${rsv_dir}/${type}/trimmed_strict_type_result.txt
  trimmed_loose_type_result=${rsv_dir}/${type}/trimmed_loose_type_result.txt
  rm ${trimmed_strict_type_result} ${trimmed_loose_type_result}

  # trim reads
  if [ $type == "A" ]; then
    cutadapt -g file:${RSV_A_forward_primer_file} -G file:${RSV_A_reverse_primer_file} \
      --pair-adapters --pair-filter any \
      --untrimmed-output ${R1_untrimmed_strict} --untrimmed-paired-output ${R2_untrimmed_strict} \
      -o ${R1_trimmed_strict} -p ${R2_trimmed_strict} \
      -j 10 ${R1_file} ${R2_file}

    cutadapt -g file:${RSV_A_forward_primer_file} -G file:${RSV_A_reverse_primer_file} \
      --pair-filter any \
      --untrimmed-output ${R1_untrimmed_loose} --untrimmed-paired-output ${R2_untrimmed_loose} \
      -o ${R1_trimmed_loose} -p ${R2_trimmed_loose} \
      -j 10 ${R1_file} ${R2_file}
  else
    cutadapt -g file:${RSV_B_forward_primer_file} -G file:${RSV_B_reverse_primer_file} \
      --pair-adapters --pair-filter any \
      --untrimmed-output ${R1_untrimmed_strict} --untrimmed-paired-output ${R2_untrimmed_strict} \
      -o ${R1_trimmed_strict} -p ${R2_trimmed_strict} \
      -j 10 ${R1_file} ${R2_file}

    cutadapt -g file:${RSV_B_forward_primer_file} -G file:${RSV_B_reverse_primer_file} \
      --pair-filter any \
      --untrimmed-output ${R1_untrimmed_loose} --untrimmed-paired-output ${R2_untrimmed_loose} \
      -o ${R1_trimmed_loose} -p ${R2_trimmed_loose} \
      -j 10 ${R1_file} ${R2_file}
  fi

  # run panmap
  /private/groups/corbettlab/alan/panmap/build/bin/panmap ${panman_file} ${R1_trimmed_strict} ${R2_trimmed_strict}  --place-per-read --redo-read-threshold 0 --em-filter-round 2 --remove-threshold 0.01 --rounds-remove 5 --preem-filter-method mbc --save-kminmer-binary-coverage --prefix ${rsv_dir}/${type}/trimmed_strict > ${rsv_dir}/${type}/trimmed_strict_panmap.log 2> ${rsv_dir}/${type}/trimmed_strict_panmap.err
  /private/groups/corbettlab/alan/panmap/build/bin/panmap ${panman_file} ${R1_trimmed_loose} ${R2_trimmed_loose}  --place-per-read --redo-read-threshold 0 --em-filter-round 2 --remove-threshold 0.01 --rounds-remove 5 --preem-filter-method mbc --save-kminmer-binary-coverage --prefix ${rsv_dir}/${type}/trimmed_loose > ${rsv_dir}/${type}/trimmed_loose_panmap.log 2> ${rsv_dir}/${type}/trimmed_loose_panmap.err

  declare -A trimmed_strict_type_proportion=( ["A"]=0 ["B"]=0 )
  declare -A trimmed_loose_type_proportion=( ["A"]=0 ["B"]=0 )

  while IFS= read -r line
  do
    haplotype=$(echo -e "$line" | cut -f1 | cut -f1 -d ',')
    proportion=$(echo -e "$line" | cut -f2 | cut -f1 -d ',')

    haplotype_cleaned=$(echo $haplotype | tr '.' '_')
    distance_to_A=$(bash /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/scripts/calculate_errors_from_mafft.sh ${RSV_fasta_dir}/${haplotype_cleaned}.fasta ${RSV_A_ref_fasta} | head -n1 | cut -f1 -d ' ')
    distance_to_B=$(bash /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/scripts/calculate_errors_from_mafft.sh ${RSV_fasta_dir}/${haplotype_cleaned}.fasta ${RSV_B_ref_fasta} | head -n1 | cut -f1 -d ' ')
    
    if [ $distance_to_A -lt $distance_to_B ]; then
      trimmed_strict_type_proportion['A']=$(echo "${trimmed_strict_type_proportion['A']} + $proportion" | bc)
      echo "$line\tA" >> ${trimmed_strict_type_result}
    else
      trimmed_strict_type_proportion['B']=$(echo "${trimmed_strict_type_proportion['B']} + $proportion" | bc)
      echo "$line\tB" >> ${trimmed_strict_type_result}
    fi
  done < trimmed_strict.abundance
  echo "A: ${trimmed_strict_type_proportion['A']}" >> ${trimmed_strict_type_result}
  echo "B: ${trimmed_strict_type_proportion['B']}" >> ${trimmed_strict_type_result}

  while IFS= read -r line
  do
    haplotype=$(echo -e "$line" | cut -f1 | cut -f1 -d ',')
    proportion=$(echo -e "$line" | cut -f2 | cut -f1 -d ',')

    haplotype_cleaned=$(echo $haplotype | tr '.' '_')
    distance_to_A=$(bash /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/scripts/calculate_errors_from_mafft.sh ${RSV_fasta_dir}/${haplotype_cleaned}.fasta ${RSV_A_ref_fasta} | head -n1 | cut -f1 -d ' ')
    distance_to_B=$(bash /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/scripts/calculate_errors_from_mafft.sh ${RSV_fasta_dir}/${haplotype_cleaned}.fasta ${RSV_B_ref_fasta} | head -n1 | cut -f1 -d ' ')
    
    if [ $distance_to_A -lt $distance_to_B ]; then
      trimmed_loose_type_proportion['A']=$(echo "${trimmed_loose_type_proportion['A']} + $proportion" | bc)
      echo "$line\tA" >> ${trimmed_loose_type_result}
    else
      trimmed_loose_type_proportion['B']=$(echo "${trimmed_loose_type_proportion['B']} + $proportion" | bc)
      echo "$line\tB" >> ${trimmed_loose_type_result}
    fi
  done < trimmed_loose.abundance

  echo "A: ${trimmed_loose_type_proportion['A']}" >> ${trimmed_loose_type_result}
  echo "B: ${trimmed_loose_type_proportion['B']}" >> ${trimmed_loose_type_result}

  rm ${R1_trimmed_strict} ${R2_trimmed_strict} ${R1_trimmed_loose} ${R2_trimmed_loose}
  rm ${R1_untrimmed_strict} ${R2_untrimmed_strict} ${R1_untrimmed_loose} ${R2_untrimmed_loose}
done
