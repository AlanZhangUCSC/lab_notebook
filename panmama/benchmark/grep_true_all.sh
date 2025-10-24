SCORES_DIR=$1

files_to_process=$(ls $SCORES_DIR/*nodeScores.tsv | sort -t'_' -k1,1 -k2,2 -k3,3n -k4,4n -k5,5n -k6,6n)

for file in $files_to_process; do
  filename=$(basename "$file")
  read tree treesize seqType numhap nummut permut numreads rep <<< $(echo "$filename" | tr '_' ' ')
  if [[ "$treesize" == "optimized" ]]; then
    treesize="20K"
  fi

  echo ">${filename}"
  if [[ "$filename" == sars_optimized* ]]; then
    bash /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/grep_true.sh \
    "$file" \
    "/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K"
  elif [[ "$filename" == sars_8M* ]]; then
    bash /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/grep_true.sh \
    "$file" \
    "/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M"
  elif [[ "$filename" == hiv_optimized* ]]; then
    bash /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/grep_true.sh \
    "$file" \
    "/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_hiv"
  fi
done