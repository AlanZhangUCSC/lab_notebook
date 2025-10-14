#!/bin/bash

#SBATCH --job-name=fast-mode-8M
#SBATCH --mail-user=bzhan146@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --mem=2tb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --output=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.log
#SBATCH --error=/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/logs/%x.%A.%a.%j.err
#SBATCH --partition=medium
#SBATCH --time=02:00:00
#SBATCH --array=0-3




mapfile -t combinations < <(python3 /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gencomb.py \
  --snps 0 \
  --haplotypes 100 \
  --percent-mutated 1.0 \
  --seq-types amplicon \
  --num-reads 1000000 2000000 \
  --num-rep 2 | tail -n +2)

# amplicon        100     0       0       1000000 0
# amplicon        100     0       0       1000000 1
# amplicon        100     0       0       2000000 0
# amplicon        100     0       0       2000000 1


read seqType numhap numsnps percentmutated numreads rep <<< "${combinations[$SLURM_ARRAY_TASK_ID]}"

prefix="${seqType}_${numhap}_${numsnps}_${percentmutated}_${numreads}_${rep}"
bash /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/genreads.sh \
  --seqtype $seqType \
  --numhap $numhap \
  --numsnp $numsnps \
  --permut $percentmutated \
  --numreads $numreads \
  --rep $rep \
  --cpus 8 \
  --panmap /private/groups/corbettlab/alan/panmap/ \
  --panman /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  --pmi /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  --random-seed $prefix \
  --swampy /private/home/bzhan146/tools/SWAMPy/src/simulate_metagenome.py \
  --reference-primer-bed-file /private/home/bzhan146/tools/SWAMPy/primer_sets/nimagenV2.bed \
  --reference-fasta-file /private/home/bzhan146/tools/SWAMPy/ref/MN908947.3.fasta \
  --jvarkit /private/home/bzhan146/tools/jvarkit/dist/jvarkit.jar \
  --out-prefix /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/fast_mode_test/${prefix}

docker run --rm \
  -v "/private/groups/corbettlab/alan/panmap/:/panmap" \
  -v "/private/groups/corbettlab/alan/panmap/panmans/:/panmans" \
  -v "/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/fast_mode_test/:/data" \
  -v "/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/fast_mode_test/:/output" \
  -w /panmap \
  --user "$(id -u):$(id -g)" \
  panmap-dev \
  bash -c "/panmap/build/bin/panmap /panmans/sars_8M.panman \
          /data/${prefix}.trimmed.fastq \
          -m /panmans/sars_8M.pmai \
          --low-memory \
          --cpus 64 \
          --prefix /output/${prefix} > /dev/null 2> /output/${prefix}.panmap.log"

rm /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/fast_mode_test/${prefix}.trimmed.fastq
rm /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/fast_mode_test/${prefix}.mutation_info.txt
