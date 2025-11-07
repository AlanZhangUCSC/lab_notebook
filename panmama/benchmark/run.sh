#!/bin/bash

#SBATCH --job-name=gen-data
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.err
#SBATCH --partition=short
#SBATCH --time=01:00:00


rm -r data_sars_20K_clustered/*perfect* &
rm -r data_sars_8M_clustered/*perfect*