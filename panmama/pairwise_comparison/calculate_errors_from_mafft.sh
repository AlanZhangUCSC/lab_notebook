fasta1=$1
fasta2=$2

cat $fasta1 $fasta2 | mafft - 2> /dev/null | \
python3 /private/groups/corbettlab/alan/lab_notebook/panmama/pairwise_comparison/calculate_errors_from_mafft.py 2> /dev/null | \
tail -n5
