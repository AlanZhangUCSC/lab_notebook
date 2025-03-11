#! /bin/bash

echo '| RSV sample | RSV A | RSV B | Num reads |'
echo '| --- | --- | --- | --- |'

for rsv_dir in /private/groups/corbettlab/alan/lab_notebook/panmama/real_data/data/RSV00*
do
  if [ ! -d "${rsv_dir}" ]; then
    continue
  fi
  rsv_basename=$(basename ${rsv_dir})

  outstring="| ${rsv_basename} |"
  trimmed_merged_type_result=${rsv_dir}/trimmed_merged_type_result.txt
  if [ ! -f "${trimmed_merged_type_result}" ]; then
    continue
  fi
  type_A=$(tail -n2 ${trimmed_merged_type_result} | head -n1 | cut -f 2 -d ' ' | xargs printf "%.3f")
  type_B=$(tail -n1 ${trimmed_merged_type_result} | head -n1 | cut -f 2 -d ' ' | xargs printf "%.3f")
  if [ "${type_A}" == "${type_B}" ]; then
    echo "$outstring | | |"
    continue
  fi

  num_reads=$(grep 'Total reads' ${rsv_dir}/trimmed_merged_panmap.log | cut -f 3 -d ' ')



  echo "${outstring} ${type_A} | ${type_B} | ${num_reads} |"
done
