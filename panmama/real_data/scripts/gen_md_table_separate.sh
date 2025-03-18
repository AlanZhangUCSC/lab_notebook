#! /bin/bash

echo '| RSV sample | RSV A \| A primers | RSV B \| A primers | RSV A primer num reads | RSV A \| B primers | RSV B \| B primers | RSV B primer num reads | potentially mixed | Marc mixed | '
echo '| --- | --- | --- | --- | --- | --- | --- | --- | --- |'

for rsv_dir in /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data/RSV00*
do
  if [ ! -d "${rsv_dir}" ]; then
    continue
  fi
  rsv_basename=$(basename ${rsv_dir})
  
  marc_mixed=""
  case "${rsv_basename}" in
    RSV00062|RSV00072|RSV00113|RSV00120|RSV00166|RSV00171|RSV00174|RSV00189|RSV00216|RSV00220|RSV00226|RSV00250)
      marc_mixed="true"
      ;;
    *)
      marc_mixed="false"
      ;;
  esac

  outstring="| ${rsv_basename} |"

  trimmed_A_type_result=${rsv_dir}/A/trimmed_A_type_result.txt
  trimmed_B_type_result=${rsv_dir}/B/trimmed_B_type_result.txt
  if [ ! -f "${trimmed_A_type_result}" ] || [ ! -f "${trimmed_B_type_result}" ]; then
    echo "$outstring | | | | | | | | | ${marc_mixed}|"
    continue
  fi
  type_A_from_A_primers=$(tail -n2 ${trimmed_A_type_result} | head -n1 | cut -f 2 -d ' ' | xargs printf "%.3f")
  type_B_from_A_primers=$(tail -n1 ${trimmed_A_type_result} | head -n1 | cut -f 2 -d ' ' | xargs printf "%.3f")
  type_A_from_B_primers=$(tail -n2 ${trimmed_B_type_result} | head -n1 | cut -f 2 -d ' ' | xargs printf "%.3f")
  type_B_from_B_primers=$(tail -n1 ${trimmed_B_type_result} | head -n1 | cut -f 2 -d ' ' | xargs printf "%.3f")


  num_reads_A=$(grep 'Total reads' ${rsv_dir}/A/trimmed_A_panmap.log | cut -f 3 -d ' ')
  num_reads_B=$(grep 'Total reads' ${rsv_dir}/B/trimmed_B_panmap.log | cut -f 3 -d ' ')
  
  if [ -z "${num_reads_A// }" ]; then
    num_reads_A=0
  fi
  
  if [ -z "${num_reads_B// }" ]; then
    num_reads_B=0
  fi

  mixed=""
  if [ "${type_A_from_A_primers}" == "${type_A_from_B_primers}" ] && [ "${type_B_from_A_primers}" == "${type_B_from_B_primers}" ]; then
    mixed="false"
  else
    if [ "${num_reads_A}" -eq 0 ] || [ "${num_reads_B}" -eq 0 ]; then
      mixed="false"
    else
      mixed="true"
    fi
  fi



  echo "${outstring} ${type_A_from_A_primers} | ${type_B_from_A_primers} | ${num_reads_A} | ${type_A_from_B_primers} | ${type_B_from_B_primers} | ${num_reads_B} | ${mixed} | ${marc_mixed}"
done
