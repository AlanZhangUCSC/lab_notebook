#!/bin/bash
FAMILY=$1

assemblies=$(curl -s "https://dfam.org/api/families/${FAMILY}/assemblies" \
  | grep -oP '"id"\s*:\s*"\K[^"]+')

for asm in $assemblies; do
  for nrph in true false; do
    echo "Downloading ${FAMILY} / ${asm} / nrph=${nrph}..."
    wget --content-disposition \
      -O "${FAMILY}_${asm}_nrph-${nrph}.tsv.gz" \
      "https://dfam.org/api/families/${FAMILY}/assemblies/${asm}/annotations?nrph=${nrph}&download=true"
    sleep 1
  done
done