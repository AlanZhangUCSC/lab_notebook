#!/bin/bash

#SBATCH --job-name=gen_stats_mito_out
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/logs/%x.%A.%a.%j.err
#SBATCH --partition=medium
#SBATCH --time=04:00:00

set -x


conda run -n snakemake python3 gen_stats_mito_out.py  ancient_mito_out/ mito_assignment_stats.MPS0.5_entropyk5.tsv &

conda run -n snakemake python3 gen_stats_mito_out.py ancient_mito_out_MPS0.5_entropyk4/ mito_assignment_stats.MPS0.5_entropyk4.tsv &

conda run -n snakemake python3 gen_stats_mito_out.py ancient_mito_out_MPS0.6_entropyk4/ mito_assignment_stats.MPS0.6_entropyk4.tsv &

conda run -n snakemake python3 gen_stats_mito_out.py ancient_mito_out_MPS0.6_entropyk5/ mito_assignment_stats.MPS0.6_entropyk5.tsv &

conda run -n snakemake python3 gen_stats_mito_out.py ancient_mito_out_MPS0.7_entropyk5/ mito_assignment_stats.MPS0.7_entropyk5.tsv &

wait


