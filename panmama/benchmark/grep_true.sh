NODE_SCORES_PATH=$1
DATA_DIR=$2

sort -k 4,4 -gr $NODE_SCORES_PATH > ${NODE_SCORES_PATH}.sorted.tmp


node_scores_path_base=$(basename "$NODE_SCORES_PATH" .nodeScores.tsv)


abundance_prefix=$(echo $node_scores_path_base | cut -f 3- -d '_')


cut -f 1 ${DATA_DIR}/${abundance_prefix}.abundance.txt | grep -F -n -f - ${NODE_SCORES_PATH}.sorted.tmp | tr ':' '\t' | column -t


rm ${NODE_SCORES_PATH}.sorted.tmp