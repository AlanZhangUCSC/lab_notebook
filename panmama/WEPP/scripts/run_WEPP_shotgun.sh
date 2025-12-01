#!/bin/bash

set -euo pipefail

usage() {
  cat << EOF
Usage: $0 <tree_path> <read_1> <read_2> <ref_path> <project_name> <output_prefix> <modified_WEPP_dir> <output_dir>

Arguments:
  tree_path          Path to phylogenetic tree file
  read_1             Path to forward reads (R1) FASTQ file
  read_2             Path to reverse reads (R2) FASTQ file
  ref_path           Path to reference genome FASTA file
  project_name       Name for this analysis project
  output_prefix      Prefix for output files
  modified_WEPP_dir  Path to modified WEPP directory
  output_dir        Directory to store results

Example:
  $0 tree.pb reads_R1.fq.gz reads_R2.fq.gz ref.fasta my_project output ./WEPP output_dir
EOF
  exit 1
}

if [ $# -ne 8 ]; then
  echo "Error: Expected 8 arguments, got $#" >&2
  usage
fi

tree_path="$1"
read_1="$2"
read_2="$3"
ref_path="$4"
project_name="$5"
output_prefix="$6"
modified_WEPP_dir="$7"
output_dir="$8"

if [ ! -f "${tree_path}" ]; then
  echo "Error: Tree file not found: ${tree_path}" >&2
  exit 1
fi

if [ ! -f "${read_1}" ]; then
  echo "Error: Read 1 file not found: ${read_1}" >&2
  exit 1
fi

if [ ! -f "${read_2}" ]; then
  echo "Error: Read 2 file not found: ${read_2}" >&2
  exit 1
fi

if [ ! -f "${ref_path}" ]; then
  echo "Error: Reference file not found: ${ref_path}" >&2
  exit 1
fi

if [ ! -d "${modified_WEPP_dir}" ]; then
  echo "Error: Modified WEPP directory not found: ${modified_WEPP_dir}" >&2
  exit 1
fi


required_files=(
  "${modified_WEPP_dir}/src/WEPP/qc_preprocess.py"
  "${modified_WEPP_dir}/config/config.yaml"
  "${modified_WEPP_dir}/workflow/rules/qc.smk"
)

for file in "${required_files[@]}"; do
  if [ ! -f "${file}" ]; then
    echo "Error: Required WEPP file not found: ${file}" >&2
    exit 1
  fi
done

tree_path=$(realpath "${tree_path}")
read_1=$(realpath "${read_1}")
read_2=$(realpath "${read_2}")
ref_path=$(realpath "${ref_path}")
modified_WEPP_dir=$(realpath "${modified_WEPP_dir}")


output_dir=$(realpath "${output_dir}")
results_real_path="${output_dir}/results"
mkdir -p "${results_real_path}"

tree_path_in_docker="/input/$(basename "${tree_path}")"
read_1_in_docker="/input/$(basename "${read_1}")"
read_2_in_docker="/input/$(basename "${read_2}")"
ref_path_in_docker="/input/$(basename "${ref_path}")"
tree_basename=$(basename "${tree_path}")
ref_basename=$(basename "${ref_path}")
uid=$(id -u)
gid=$(id -g)

echo "Running WEPP analysis:"
echo "  Project: ${project_name}"
echo "  Output prefix: ${output_prefix}"
echo "  Tree: ${tree_path}"
echo "  Reads: ${read_1}, ${read_2}"
echo "  Reference: ${ref_path}"
echo ""

docker run --rm \
  -v "${modified_WEPP_dir}/src/WEPP/qc_preprocess.py:/WEPP/src/WEPP/qc_preprocess.py" \
  -v "${modified_WEPP_dir}/config/config.yaml:/WEPP/config/config.yaml" \
  -v "${modified_WEPP_dir}/workflow/rules/qc.smk:/WEPP/workflow/rules/qc.smk" \
  -v "${results_real_path}:/WEPP/results" \
  -v "${tree_path}:${tree_path_in_docker}:ro" \
  -v "${read_1}:${read_1_in_docker}:ro" \
  -v "${read_2}:${read_2_in_docker}:ro" \
  -v "${ref_path}:${ref_path_in_docker}:ro" \
  pranavgangwar/wepp:latest \
  bash -c "
    cd /WEPP
    mkdir -p data/${project_name}

    cp ${tree_path_in_docker} \
       ${ref_path_in_docker} \
       ${read_1_in_docker} \
       ${read_2_in_docker} \
       data/${project_name}/

    snakemake \
      --config DIR=${project_name} \
               FILE_PREFIX=${output_prefix} \
               TREE=${tree_basename} \
               REF=${ref_basename} \
               CLADE_IDX=-1 \
               SHOTGUN=True \
      --cores 32 \
      --use-conda

    chown -R ${uid}:${gid} results/
  "

