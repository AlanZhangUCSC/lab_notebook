#!/bin/bash

echo -e "num_hap\tnum_snp\tseq_type\tdepth\tmin_rank\tmax_rank\tmedian_rank"
for file in out/archive/*_*_1_*.nodeRank.txt; do
  read -r num_hap num_snp seq_type depth <<< $(basename $file | cut -f 1 -d '.' | tr '_' ' ')
  read -r min max median <<< $(awk '{val=$4; if(NR==1){min=max=val} if(val<min) min=val; if(val>max) max=val; a[NR]=val} END{asort(a); n=length(a); if(n%2==1) med=a[(n+1)/2]; else med=(a[n/2]+a[n/2+1])/2; print min, max, med}' "$file")
  echo -e "$num_hap\t$num_snp\t$seq_type\t$depth\t$min\t$max\t$median"
done