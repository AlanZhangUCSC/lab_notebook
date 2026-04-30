import os
import sys
import signal
import argparse
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

def cleanup(out_dir):
  print(f'Cleaning up {out_dir}')
  for f in os.listdir(out_dir):
    if f.endswith('.mapped.sorted.bam'):
      os.remove(os.path.join(out_dir, f))

def align_family(refs, query_file, ref_dir, out_dir):
  out_bams = []
  for ref in refs:
    fasta_file = os.path.join(ref_dir, ref + '.fasta')
    if not os.path.exists(fasta_file):
      continue
    out_prefix = os.path.join(out_dir, ref)
    out_bam = out_prefix + '.mapped.sorted.bam'
    subprocess.run(['/home/alan/tools/misc/align_ancient', fasta_file, query_file, out_prefix], check=True)
    if not os.path.exists(out_bam):
      raise FileNotFoundError(f'Error: {out_bam} not found')
      exit(1)
    out_bams.append(out_bam)

  if len(out_bams) == 0:
    return
  family = os.path.basename(query_file).split('.')[0]
  merged_bam = os.path.join(out_dir, family + '.merged.bam')
  subprocess.run(['samtools', 'merge', '-f', '-o', merged_bam, *out_bams], check=True)
  for f in out_bams:
    os.remove(f)

  merged_bam_lines = subprocess.run(['samtools', 'view', merged_bam], capture_output=True, text=True).stdout.split('\n')
  if len([l for l in merged_bam_lines if l]) == 0:
    os.remove(merged_bam)
    return

parser = argparse.ArgumentParser(description='Run alignments for multiple families')
parser.add_argument('ref_dir', type=str, help='Directory containing reference fasta files')
parser.add_argument('query_dir', type=str, help='Directory containing query fastq files')
parser.add_argument('meta_data', type=str, help='Path to metadata file')
parser.add_argument('--out_dir', type=str, default='align_family_out', help='Directory to store output BAM files (default: align_family_out)')
parser.add_argument('-t', '--threads', type=int, default=4, help='Number of parallel workers (default: 4)')
args = parser.parse_args()

ref_dir = args.ref_dir
query_dir = args.query_dir
meta_data = args.meta_data
out_dir = args.out_dir
threads = args.threads

os.makedirs(out_dir, exist_ok=True)

for sig in (signal.SIGINT, signal.SIGTERM):
  signal.signal(sig, lambda s, f: (cleanup(out_dir), sys.exit(1)))

fastq_files = [os.path.join(query_dir, f) for f in os.listdir(query_dir) if f.endswith('.fastq')]

node_to_family = defaultdict(str)
family_to_nodes = defaultdict(list)
with open(meta_data, 'r') as f:
  f.readline()
  for line in f:
    fields = line.strip().split('\t')
    node = fields[0]
    family = fields[4]
    node_to_family[node] = family
    family_to_nodes[family].append(node)

with ThreadPoolExecutor(max_workers=threads) as executor:
  futures = {
    executor.submit(align_family, family_to_nodes[os.path.basename(f).split('.')[0]], f, ref_dir, out_dir): f
    for f in fastq_files
  }
  for future in as_completed(futures):
    fastq_file = futures[future]
    try:
      future.result()
      print(f'Finished: {fastq_file}')
    except Exception as e:
      print(f'Error processing {fastq_file}: {e}')
      cleanup(out_dir)
      sys.exit(1)