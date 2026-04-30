#! /bin/bash

set -x

input_fasta=$1
output_prefix=$2
iterations=$3

dipper_cpu -i r -I $input_fasta -O $output_prefix.it0.nwk 

for i in $(seq 1 $iterations); do
  if [ $i -eq 1 ]; then
	  twilight -t $output_prefix.it$((i-1)).nwk -i $input_fasta -o $output_prefix.it${i}.aln -c -v -w --check -r 0.999 --cpu-only -C 32 > $output_prefix.it${i}.twilight.log 2> $output_prefix.it${i}.twilight.err
	else
	  twilight -t $output_prefix.it$((i-1)).nwk -i $output_prefix.it$((i-1)).aln.gz -o $output_prefix.it${i}.aln -c -v -w --check -r 0.999 --cpu-only -C 32 > $output_prefix.it${i}.twilight.log 2> $output_prefix.it${i}.twilight.err
  fi
  dipper_cpu -i m -I $output_prefix.it${i}.aln.gz -O $output_prefix.it${i}.nwk -d 4 --threads 32
done
