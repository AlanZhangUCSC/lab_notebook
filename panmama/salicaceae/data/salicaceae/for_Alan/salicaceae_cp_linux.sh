work_path=/home/rsd717/rsd717/Myriophyllum/modern_panel/cpDNA
cd $work_path
script_path=/home/rsd717/rsd717/scripts
micromamba activate genome


# Split outgroup plastomes into LSC / IRb / SSC multifasta files
mkdir -p $work_path/outgroup_cpstools
outgroup_cp=$work_path/NC_034285.1.fasta

cpstools IR -i $outgroup_cp > $work_path/outgroup_cpstools/outgroup.cpstools_four_sec.txt
python $script_path/cp_helpers.py orient-outgroup \
  --input-fasta $outgroup_cp \
  --cpstools-txt $work_path/outgroup_cpstools/outgroup.cpstools_four_sec.txt \
  --output-fasta $work_path/outgroup_cpstools/outgroup.LSC_start.fasta \
  --output-regions $work_path/outgroup_cpstools/outgroup.regions.tsv

# Split all plastomes into LSC / IRb / SSC multifasta files
mkdir -p $work_path/partitions

python $script_path/cp_helpers.py split-partitions \
  --sample-regions-tsv $work_path/salicaceae_fastas.final.oriented.tsv \
  --sample-fasta-dir $work_path/final_oriented_fastas \
  --outgroup-fasta $work_path/outgroup_cpstools/outgroup.LSC_start.fasta \
  --outgroup-regions-tsv $work_path/outgroup_cpstools/outgroup.regions.tsv \
  --outdir $work_path/partitions

# Run partition-wise MAFFT with fast reverse-complement detection
mkdir -p $work_path/align_trim
cd $work_path/align_trim

(mafft --thread 24 --auto --adjustdirection $work_path/partitions/LSC.fasta > $work_path/align_trim/LSC.aligned.fasta && \
python $script_path/cp_helpers.py normalize-mafft --input-fasta $work_path/align_trim/LSC.aligned.fasta && \
trimal -automated1 -in $work_path/align_trim/LSC.aligned.fasta -out $work_path/align_trim/LSC_trimmed.fasta ) > $work_path/align_trim/LSC.log 2>&1 &

(mafft --thread 24 --auto --adjustdirection $work_path/partitions/IRb.fasta > $work_path/align_trim/IRb.aligned.fasta && \
python $script_path/cp_helpers.py normalize-mafft --input-fasta $work_path/align_trim/IRb.aligned.fasta && \
trimal -automated1 -in $work_path/align_trim/IRb.aligned.fasta -out $work_path/align_trim/IRb_trimmed.fasta ) > $work_path/align_trim/IR.log 2>&1 &

(mafft --thread 24 --auto --adjustdirection $work_path/partitions/SSC.fasta > $work_path/align_trim/SSC.aligned.fasta && \
python $script_path/cp_helpers.py normalize-mafft --input-fasta $work_path/align_trim/SSC.aligned.fasta && \
trimal -automated1 -in $work_path/align_trim/SSC.aligned.fasta -out $work_path/align_trim/SSC_trimmed.fasta ) > $work_path/align_trim/SSC.log 2>&1 &


# Concatenate the three trimmed backbone partitions and build a guide tree
mkdir -p $work_path/tree
cd $work_path/tree

python $script_path/cp_helpers.py concat-backbone \
  --lsc-trimmed $work_path/align_trim/LSC_trimmed.fasta \
  --irb-trimmed $work_path/align_trim/IRb_trimmed.fasta \
  --ssc-trimmed $work_path/align_trim/SSC_trimmed.fasta \
  --output-fasta $work_path/align_trim/backbone_trimmed_concat.fasta \
  --output-partitions $work_path/align_trim/backbone_trimmed_concat.partitions.tsv

nohup iqtree -s $work_path/align_trim/backbone_trimmed_concat.fasta -m MFP -B 1000 -alrt 1000 -nt AUTO -o NC_034285.1 -pre guide_backbone &
micromamba activate popgen
OrganPath RenameTree -i guide_backbone.treefile -m acc2tax.txt -o rename_guide_backbone.treefile

# Recover >100 bp continuous fragments removed by trimAl
mkdir -p $work_path/recovered_chunks
cd $work_path/recovered_chunks

python $script_path/cp_helpers.py recover-chunks \
  --lsc-aligned $work_path/align_trim/LSC.aligned.fasta \
  --lsc-trimmed $work_path/align_trim/LSC_trimmed.fasta \
  --irb-aligned $work_path/align_trim/IRb.aligned.fasta \
  --irb-trimmed $work_path/align_trim/IRb_trimmed.fasta \
  --ssc-aligned $work_path/align_trim/SSC.aligned.fasta \
  --ssc-trimmed $work_path/align_trim/SSC_trimmed.fasta \
  --outdir $work_path/recovered_chunks \
  --min-len 100

# Run progressive Mauve on the recovered chunks
mkdir -p $work_path/recovered_blocks
progressiveMauve \
  --output=$work_path/recovered_blocks/recovered.xmfa \
  $work_path/recovered_chunks/*.fa

python $script_path/cp_helpers.py extract-blocks \
  --xmfa $work_path/recovered_blocks/recovered.xmfa \
  --outdir $work_path/recovered_blocks \
  --outgroup NC_034285.1 \
  --threads 12

python $script_path/cp_helpers.py concat-final \
  --backbone-fasta $work_path/align_trim/backbone_trimmed_concat.fasta \
  --blocks-manifest $work_path/recovered_blocks/recovered_blocks.tsv \
  --output-fasta $work_path/Salicaceae_cp_final_MSA.fasta \
  --output-partitions $work_path/Salicaceae_cp_final_partitions.tsv


nohup iqtree -s $work_path/Salicaceae_cp_final_MSA.fasta \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -nt AUTO \
  -o NC_034285.1 \
  -pre $work_path/Salicaceae_cp_final &
