import os
import sys
import subprocess

def merge_to_ranges(nums):
  if not nums:
    return []
  
  ranges = []
  start = nums[0]
  end = nums[0]
  
  for i in range(1, len(nums)):
    if nums[i] == end + 1:
      end = nums[i]
    else:
      ranges.append([start, end])
      start = end = nums[i]
  
  ranges.append([start, end])
  return ranges

read_assignments_path = sys.argv[1]
fastq_path = sys.argv[2]

fastq_filtered = f'{fastq_path}.filtered.fq'
filter_cmd = f'~/tools/BBTools/bbduk.sh in={fastq_path} out={fastq_filtered} entropy=0.7'
subprocess.run(filter_cmd, shell=True, check=True)

fastq_filtered_sorted = f'{fastq_path}.filtered.sorted.fq'
sort_cmd = f"cat {fastq_filtered} | paste - - - - | sort -k1,1 -S 3G | tr '\\t' '\\n' > {fastq_filtered_sorted}"
subprocess.run(sort_cmd, shell=True, check=True)

fastq_path = fastq_filtered_sorted

with open(read_assignments_path, 'r') as fh:
  for line in fh.readlines():
    fields = line.strip().split()
    node = fields[0]
    indices = list(map(int, fields[2].split(',')))
    indices.sort()

    if (len(indices) == 0): continue

    merged_indices = merge_to_ranges(indices)

    fastq_out_path = f'{node}.fq'
    fastq_fh = open(fastq_out_path, 'w')
    for start, end in merged_indices:
      cmd = f'sed -n "{start*4+1},{end*4+4}p" {fastq_path}'
      reads = subprocess.run(cmd, shell=True, check=True, capture_output=True)
      fastq_fh.write(reads.stdout.decode('utf-8'))
    fastq_fh.close()

os.remove(fastq_filtered)
os.remove(fastq_filtered_sorted)