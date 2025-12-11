sample="$1"
output_dir="$2"

python3 /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/assign_reads_by_indices.py \
  ${output_dir}/${sample}.mammoth.read_scores.mgsr.assignedReads.out \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/final_analysis/data/reads/${sample}.mammoth.fq \
  ${sample}

python3 /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/gen_read_placement_table.py \
  ${output_dir}/${sample}.mammoth.read_scores.mgsr.abundance.out \
  ${output_dir}/${sample}.mammoth.read_scores.mgsr.assignedReads.out \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/final_analysis/data/reads/${sample}.mammoth.fq \
  ${sample}.read_placement.tsv

cp ${output_dir}/${sample}.mammoth.read_scores.read_scores_info.tsv .
cp ${output_dir}/${sample}.mammoth.read_scores.mgsr.abundance.out .