#!/bin/bash

#SBATCH --job-name=rsv_real_ivar_trim
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --array=0-191%50

set -x

rsv_samples=(RSV00059 RSV00060 RSV00061 RSV00062 RSV00063 RSV00064 RSV00065 RSV00066 RSV00067 RSV00068 RSV00069 RSV00070 RSV00071 RSV00072 RSV00073 RSV00074 RSV00075 RSV00076 RSV00077 RSV00078 RSV00079 RSV00080 RSV00081 RSV00082 RSV00083 RSV00084 RSV00085 RSV00086 RSV00087 RSV00088 RSV00089 RSV00090 RSV00091 RSV00092 RSV00093 RSV00094 RSV00095 RSV00096 RSV00097 RSV00098 RSV00099 RSV00100 RSV00101 RSV00102 RSV00103 RSV00104 RSV00105 RSV00106 RSV00107 RSV00108 RSV00109 RSV00111 RSV00112 RSV00113 RSV00114 RSV00115 RSV00116 RSV00117 RSV00118 RSV00119 RSV00120 RSV00121 RSV00122 RSV00123 RSV00124 RSV00125 RSV00126 RSV00127 RSV00128 RSV00129 RSV00130 RSV00131 RSV00132 RSV00133 RSV00134 RSV00135 RSV00136 RSV00137 RSV00138 RSV00139 RSV00140 RSV00141 RSV00142 RSV00143 RSV00144 RSV00145 RSV00146 RSV00147 RSV00148 RSV00149 RSV00150 RSV00151 RSV00152 RSV00153 RSV00154 RSV00159 RSV00164 RSV00165 RSV00166 RSV00167 RSV00168 RSV00169 RSV00170 RSV00171 RSV00172 RSV00173 RSV00174 RSV00175 RSV00176 RSV00177 RSV00178 RSV00179 RSV00180 RSV00181 RSV00182 RSV00183 RSV00184 RSV00185 RSV00186 RSV00187 RSV00188 RSV00189 RSV00190 RSV00191 RSV00192 RSV00193 RSV00194 RSV00195 RSV00196 RSV00197 RSV00198 RSV00199 RSV00200 RSV00201 RSV00202 RSV00203 RSV00204 RSV00205 RSV00206 RSV00207 RSV00208 RSV00209 RSV00210 RSV00211 RSV00212 RSV00213 RSV00214 RSV00215 RSV00216 RSV00217 RSV00218 RSV00219 RSV00220 RSV00221 RSV00222 RSV00223 RSV00224 RSV00225 RSV00226 RSV00227 RSV00228 RSV00229 RSV00230 RSV00231 RSV00232 RSV00233 RSV00234 RSV00235 RSV00236 RSV00237 RSV00238 RSV00239 RSV00240 RSV00241 RSV00242 RSV00243 RSV00244 RSV00245 RSV00246 RSV00247 RSV00248 RSV00249 RSV00250 RSV00251 RSV00252 RSV00253 RSV00254 RSV00255 RSV00256 RSV00257 RSV00258 RSV00259)

data_dir=/private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data
RSV_A_primer_bed_file=$data_dir/other_files_you_might_need/RSVA.primer.bed
RSV_B_primer_bed_file=$data_dir/other_files_you_might_need/RSVB.primer.bed
RSV_A_ref_fasta=$data_dir/other_files_you_might_need/RS20000581_RSVA_reference.fasta
RSV_B_ref_fasta=$data_dir/other_files_you_might_need/RE20000104_RSVB_reference.fasta
RSV_fasta_dir=/private/groups/corbettlab/alan/panmap/dev/panmama_eval/rsv4000/fasta/info/fasta_unaligned

panman_file=/private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/rsv/k19_s8_null/panmap_outputs/panmap2_index/rsv4000.panman

rsv_dir=$data_dir/${rsv_samples[$SLURM_ARRAY_TASK_ID]}

A_R1_file=$(ls ${rsv_dir}/A/*RSV*R1*fastq*)
A_R2_file=$(ls ${rsv_dir}/A/*RSV*R2*fastq*)
B_R1_file=$(ls ${rsv_dir}/B/*RSV*R1*fastq*)
B_R2_file=$(ls ${rsv_dir}/B/*RSV*R2*fastq*)

sorted_A_bam=${rsv_dir}/A/pipeline_alignment_A.sorted.bam
sorted_B_bam=${rsv_dir}/B/pipeline_alignment_B.sorted.bam

trimmed_A_bam=${rsv_dir}/A/pipeline_alignment_A.trimmed.bam
trimmed_B_bam=${rsv_dir}/B/pipeline_alignment_B.trimmed.bam

trimmed_A_fastq=${rsv_dir}/A/pipeline_alignment_A.trimmed.fastq
trimmed_B_fastq=${rsv_dir}/B/pipeline_alignment_B.trimmed.fastq
trimmed_merged_fastq=${rsv_dir}/pipeline_alignment.trimmed.merged.fastq

trimmed_A_fastq_length=${rsv_dir}/A/pipeline_alignment_A.trimmed.fastq.length
trimmed_B_fastq_length=${rsv_dir}/B/pipeline_alignment_B.trimmed.fastq.length

# align reads generated from each primer type to the respective reference
bwa mem -t 8 ${RSV_A_ref_fasta} ${A_R1_file} ${A_R2_file} | samtools sort | samtools view -F 4 -o ${sorted_A_bam} &
bwa mem -t 8 ${RSV_B_ref_fasta} ${B_R1_file} ${B_R2_file} | samtools sort | samtools view -F 4 -o ${sorted_B_bam} &
wait


# index bams
samtools index --threads 8 ${sorted_A_bam} &
samtools index --threads 8 ${sorted_B_bam} &
wait

# trim reads
ivar trim -e -i ${sorted_A_bam} -b ${RSV_A_primer_bed_file} -p ${rsv_dir}/A/pipeline_alignment_A.trimmed &
ivar trim -e -i ${sorted_B_bam} -b ${RSV_B_primer_bed_file} -p ${rsv_dir}/B/pipeline_alignment_B.trimmed &
wait

# convert trimmed bams back to fastq and then merge
samtools fastq --threads 8 ${trimmed_A_bam} -o ${trimmed_A_fastq} &
samtools fastq --threads 8 ${trimmed_B_bam} -o ${trimmed_B_fastq} &
wait






