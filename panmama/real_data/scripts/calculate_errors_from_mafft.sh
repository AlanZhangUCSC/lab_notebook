fasta1=$1
fasta2=$2
tmp_aligned=$(mktemp)
tmp_distance=$(mktemp)

cat $fasta1 $fasta2 | mafft - > $tmp_aligned 2> /dev/null
python3 /private/groups/corbettlab/alan/panmama-snakemake/workflow/evals/scripts/calculate_errors_from_mafft.py $tmp_aligned > $tmp_distance 2> /dev/null
rm $tmp_aligned

tail -n4 $tmp_distance | cut -f 1-2 -d ' '

rm $tmp_distance
