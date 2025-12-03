import argparse
import os
import subprocess
import random
import shutil



parser = argparse.ArgumentParser(description='Run PanMAMA benchmark with configurable parameters')
parser.add_argument('panman_path', help='Path to PanMAN executable')
parser.add_argument('pmai_path', help='Path to PMAI executable')
parser.add_argument('num_haps', type=int, help='Number of haplotypes')
parser.add_argument('--num-reads', type=int, default=1000, help='Number of reads (default: 1000)')
parser.add_argument('--num-snps', type=int, default=0, help='Number of SNPs (default: 0)')
parser.add_argument('--percent-mutated', type=float, default=0, help='Percentage mutated (default: 0)')
parser.add_argument('--selection', type=str, default='overlap_coefficients', help='Probable node selection scheme (overlap_coefficients or node_scores) (default: coverlap coefficients)')
parser.add_argument('--prefix', type=str, help='Output prefix (default: auto-generated)')
parser.add_argument('--random-seed', type=str, default=None, help='Random seed (default: random)')
parser.add_argument('--out-dir', type=str, default='.', help='Output directory (default: current directory)')

args = parser.parse_args()

panman_path = os.path.abspath(args.panman_path)
pmai_path = os.path.abspath(args.pmai_path)
out_dir = os.path.abspath(args.out_dir)

random_seed = args.random_seed or random.randint(1, 2**31 - 1)

out_prefix = args.prefix or (
  f'{args.num_haps}_{args.num_snps}_{args.percent_mutated}_'
  f'{args.num_reads}_{random_seed}'
)

selection_scheme = args.selection
if selection_scheme not in ['overlap_coefficients', 'node_scores']:
  raise ValueError("Selection scheme must be 'overlap_coefficients' or 'node_scores'")

os.makedirs(out_dir, exist_ok=True)
tmp_dir = '/data/tmp/'
cmd = [ 'bash', 
  '/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/genreads.sh',
  '--seqtype', 'shotgun',
  '--numhap', str(args.num_haps),
  '--numsnp', str(args.num_snps),
  '--permut', str(args.percent_mutated),
  '--numreads', str(args.num_reads),
  '--rep', '1',
  '--cpus', '8',
  '--panmap', '/private/groups/corbettlab/alan/panmap/',
  '--panman', panman_path,
  '--pmi', pmai_path,
  '--random-seed', str(random_seed),
  '--swampy', '/private/home/bzhan146/tools/SWAMPy/src/simulate_metagenome.py',
  '--reference-primer-bed-file',
    '/private/home/bzhan146/tools/SWAMPy/primer_sets/nimagenV2.bed',
  '--reference-fasta-file',
    '/private/home/bzhan146/tools/SWAMPy/ref/MN908947.3.fasta',
  '--jvarkit', '/private/home/bzhan146/tools/jvarkit/dist/jvarkit.jar',
  '--out-prefix', os.path.join(tmp_dir, out_prefix),
]
subprocess.run(cmd, check=True)

shutil.move(os.path.join(tmp_dir, f'{out_prefix}.mutation_info.txt'), out_dir)
shutil.move(os.path.join(tmp_dir, f'{out_prefix}.abundance.txt'), out_dir)

panmap_path = '/private/groups/corbettlab/alan/panmap/'
panman_dir = os.path.dirname(panman_path)
pmai_dir = os.path.dirname(pmai_path)

readpath1 = os.path.join(tmp_dir, f'{out_prefix}_R1.fastq')
readpath2 = os.path.join(tmp_dir, f'{out_prefix}_R2.fastq')
panman_name = os.path.basename(panman_path)
pmai_name = os.path.basename(pmai_path)

docker_cmd = [
  'docker', 'run', '--rm',
  '-v', f'{panmap_path}:/panmap',
  '-v', f'{panman_dir}:/panmans',
  '-v', f'{pmai_dir}:/pmais',
  '-v', f'{tmp_dir}:/data',
  '-v', f'{out_dir}:/output',
  '-w', '/panmap',
  '--user', f'{os.getuid()}:{os.getgid()}',
  'panmap-dev',
]

if selection_scheme == 'overlap_coefficients':
  docker_cmd.extend([
    'bash', '-c',
    f'/panmap/build/bin/panmap /panmans/{panman_name} /data/{os.path.basename(readpath1)} '
    f'/data/{os.path.basename(readpath2)} -m /pmais/{pmai_name} --read-scores '
    f'--prefix /output/{out_prefix}.readScores --cpus 8'
  ])
else:
  docker_cmd.extend([
    'bash', '-c',
    f'/panmap/build/bin/panmap /panmans/{panman_name} /data/{os.path.basename(readpath1)} '
    f'/data/{os.path.basename(readpath2)} -m /pmais/{pmai_name} --overlap-coefficients 1000 '
    f'--prefix /output/{out_prefix}.overlapCoefficients --cpus 8'
  ])

subprocess.run(docker_cmd, check=True)

os.remove(readpath1)
os.remove(readpath2)
os.remove(os.path.join(tmp_dir, f'{out_prefix}.genreads.log'))