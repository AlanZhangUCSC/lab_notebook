#!/bin/bash

#SBATCH --job-name=rsv_real_ivar_trim_separate
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/logs/%x.%A.%a.%j.err
#SBATCH --partition=medium
#SBATCH --time=02:00:00
#SBATCH --array=82,164


set -x

rsv_samples=(RSV00059 RSV00060 RSV00061 RSV00062 RSV00063 RSV00064 RSV00065 RSV00066 RSV00067 RSV00068 RSV00069 RSV00070 RSV00071 RSV00072 RSV00073 RSV00074 RSV00075 RSV00076 RSV00077 RSV00078 RSV00079 RSV00080 RSV00081 RSV00082 RSV00083 RSV00084 RSV00085 RSV00086 RSV00087 RSV00088 RSV00089 RSV00090 RSV00091 RSV00092 RSV00093 RSV00094 RSV00095 RSV00096 RSV00097 RSV00098 RSV00099 RSV00100 RSV00101 RSV00102 RSV00103 RSV00104 RSV00105 RSV00106 RSV00107 RSV00108 RSV00109 RSV00111 RSV00112 RSV00113 RSV00114 RSV00115 RSV00116 RSV00117 RSV00118 RSV00119 RSV00120 RSV00121 RSV00122 RSV00123 RSV00124 RSV00125 RSV00126 RSV00127 RSV00128 RSV00129 RSV00130 RSV00131 RSV00132 RSV00133 RSV00134 RSV00135 RSV00136 RSV00137 RSV00138 RSV00139 RSV00140 RSV00141 RSV00142 RSV00143 RSV00144 RSV00145 RSV00146 RSV00147 RSV00148 RSV00149 RSV00150 RSV00151 RSV00152 RSV00153 RSV00154 RSV00159 RSV00164 RSV00165 RSV00166 RSV00167 RSV00168 RSV00169 RSV00170 RSV00171 RSV00172 RSV00173 RSV00174 RSV00175 RSV00176 RSV00177 RSV00178 RSV00179 RSV00180 RSV00181 RSV00182 RSV00183 RSV00184 RSV00185 RSV00186 RSV00187 RSV00188 RSV00189 RSV00190 RSV00191 RSV00192 RSV00193 RSV00194 RSV00195 RSV00196 RSV00197 RSV00198 RSV00199 RSV00200 RSV00201 RSV00202 RSV00203 RSV00204 RSV00205 RSV00206 RSV00207 RSV00208 RSV00209 RSV00210 RSV00211 RSV00212 RSV00213 RSV00214 RSV00215 RSV00216 RSV00217 RSV00218 RSV00219 RSV00220 RSV00221 RSV00222 RSV00223 RSV00224 RSV00225 RSV00226 RSV00227 RSV00228 RSV00229 RSV00230 RSV00231 RSV00232 RSV00233 RSV00234 RSV00235 RSV00236 RSV00237 RSV00238 RSV00239 RSV00240 RSV00241 RSV00242 RSV00243 RSV00244 RSV00245 RSV00246 RSV00247 RSV00248 RSV00249 RSV00250 RSV00251 RSV00252 RSV00253 RSV00254 RSV00255 RSV00256 RSV00257 RSV00258 RSV00259)
types=(A B)


data_dir=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data

all_rsv_dirs=()
all_rsv_samples=()
all_rsv_types=()
for rsv_sample in ${rsv_samples[@]}
do
  for type in ${types[@]}
  do
    all_rsv_dirs+=(${data_dir}/${rsv_sample}/${type})
    all_rsv_samples+=(${rsv_sample})
    all_rsv_types+=(${type})
  done
done



rsv_dir=${all_rsv_dirs[$SLURM_ARRAY_TASK_ID]}
rsv_sample=${all_rsv_samples[$SLURM_ARRAY_TASK_ID]}
rsv_type=${all_rsv_types[$SLURM_ARRAY_TASK_ID]}


RSV_type_primer_bed_file=$data_dir/other_files_you_might_need/RSV${rsv_type}.primer.bed
RSV_type_ref_fasta=""
RSV_A_ref_fasta=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data/other_files_you_might_need/RS20000581_RSVA_reference.fasta
RSV_B_ref_fasta=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data/other_files_you_might_need/RE20000104_RSVB_reference.fasta
if [ $rsv_type == "A" ]; then
  RSV_type_ref_fasta=${RSV_A_ref_fasta}
elif [ $rsv_type == "B" ]; then
  RSV_type_ref_fasta=${RSV_B_ref_fasta}
fi
RSV_fasta_dir=/private/groups/corbettlab/alan/panmap/dev/panmama_eval/rsv4000/fasta/info/fasta_unaligned

panman_file=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/rsv/k19_s8_null/panmap_outputs/panmap2_index/rsv4000.panman

R1_file=$(ls ${rsv_dir}/*RSV*R1*fastq*)
R2_file=$(ls ${rsv_dir}/*RSV*R2*fastq*)

sorted_bam=${rsv_dir}/pipeline_alignment_${rsv_type}.sorted.bam

trimmed_bam=${rsv_dir}/pipeline_alignment_${rsv_type}.trimmed.bam

trimmed_fastq=${rsv_dir}/pipeline_alignment_${rsv_type}.trimmed.fastq

# # align reads generated from each primer type to the respective reference
# bwa mem -t 8 ${RSV_type_ref_fasta} ${R1_file} ${R2_file} | samtools sort | samtools view -h -F 4 -o ${sorted_bam}

# # index bams
# samtools index --threads 8 ${sorted_bam}

# # trim reads
# ivar trim -e -i ${sorted_bam} -b ${RSV_type_primer_bed_file} -p ${rsv_dir}/pipeline_alignment_${rsv_type}.trimmed 

# # convert trimmed bams back to fastq and then merge
# samtools fastq --threads 8 ${trimmed_bam} -o ${trimmed_fastq}


# run panmap
/private/groups/corbettlab/alan/panmap/build/bin/panmap \
  ${panman_file} ${trimmed_fastq} \
  --place-per-read --redo-read-threshold 0 \
  --em-filter-round 2 --remove-threshold 0.01 \
  --rounds-remove 5 --preem-filter-method mbc \
  --save-kminmer-binary-coverage --prefix ${rsv_dir}/trimmed_${rsv_type} \
  > ${rsv_dir}/trimmed_${rsv_type}_panmap.log 2> ${rsv_dir}/trimmed_${rsv_type}_panmap.err

declare -A trimmed_type_proportion=( ["A"]=0 ["B"]=0 )

trimmed_type_result=${rsv_dir}/trimmed_${rsv_type}_type_result.txt
rm ${trimmed_type_result}
while IFS= read -r line
do
  haplotype=$(echo -e "$line" | cut -f1 | cut -f1 -d ',')
  proportion=$(echo -e "$line" | cut -f2 | cut -f1 -d ',')

  haplotype_cleaned=$(echo $haplotype | tr '.' '_')
  distance_to_A=$(bash /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/scripts/calculate_errors_from_mafft.sh ${RSV_fasta_dir}/${haplotype_cleaned}.fasta ${RSV_A_ref_fasta} | head -n1 | cut -f1 -d ' ')
  distance_to_B=$(bash /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/scripts/calculate_errors_from_mafft.sh ${RSV_fasta_dir}/${haplotype_cleaned}.fasta ${RSV_B_ref_fasta} | head -n1 | cut -f1 -d ' ')
  
  if [ $distance_to_A -lt $distance_to_B ]; then
    trimmed_type_proportion['A']=$(echo "${trimmed_type_proportion['A']} + $proportion" | bc)
    echo -e "$line\tA" >> ${trimmed_type_result}
  else
    trimmed_type_proportion['B']=$(echo "${trimmed_type_proportion['B']} + $proportion" | bc)
    echo -e "$line\tB" >> ${trimmed_type_result}
  fi
done < "${rsv_dir}/trimmed_${rsv_type}.abundance"

echo "A: ${trimmed_type_proportion['A']}" >> ${trimmed_type_result}
echo "B: ${trimmed_type_proportion['B']}" >> ${trimmed_type_result}

# rm ${sorted_bam} \
#   ${sorted_bam}.bai \
#   ${trimmed_bam} \
#   ${trimmed_fastq}



