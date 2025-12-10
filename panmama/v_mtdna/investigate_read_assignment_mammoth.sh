sample="$1"

python3 /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/assign_reads_by_indices.py \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/ancient_mito_out/${sample}.mammoth.read_scores.mgsr.assignedReads.out \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/final_analysis/data/reads/${sample}.mammoth.fq \
  ${sample}

python3 /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/gen_read_placement_table.py \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/ancient_mito_out/${sample}.mammoth.read_scores.mgsr.abundance.out \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/ancient_mito_out/${sample}.mammoth.read_scores.mgsr.assignedReads.out \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/final_analysis/data/reads/${sample}.mammoth.fq \
  ${sample}.read_placement.tsv

cp /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/ancient_mito_out/${sample}.mammoth.read_scores.read_scores_info.tsv .
cp /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/ancient_mito_out/${sample}.mammoth.read_scores.mgsr.abundance.out .