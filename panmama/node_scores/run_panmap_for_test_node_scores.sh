# use panmap branch mgsr-development commit fff14e03f12b609f1877d845299862cb9b71b12d

#!/bin/bash
PANMAP_EXECUTABLE="$1"
PANMAN="$2"
MGSR_INDEX="$3"
READ_SAMPLES_DIR="$4"


if [ -z "$PANMAP_EXECUTABLE" ]; then
  echo "Usage: $0 <path-to-executable> <path-to-mgsr-index> <path-to-read-samples-dir>"
  exit 1
fi

if [ ! -x "$PANMAP_EXECUTABLE" ]; then
  echo "Error: $PANMAP_EXECUTABLE is not executable or doesn't exist"
  exit 1
fi


read_prefixes=$(find $READ_SAMPLES_DIR -name "*R1.fastq" | sed 's/_R1.fastq//g')
for read_prefix in $read_prefixes; do
  "$PANMAP_EXECUTABLE" $PANMAN \
    ${read_prefix}_R1.fastq \
    ${read_prefix}_R2.fastq \
    -m $MGSR_INDEX \
    --prefix $read_prefix \
    --cpus 4  > /dev/null 2>&1
done


