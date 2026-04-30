from Bio import SeqIO
import subprocess
import shutil
import sys
import os

def run_nucmer(fa_file, seq, prefix):
  # Run nucmer and show-coords
  nucmer_cmd = ["nucmer", "--maxmatch", "-l", "500", "-c", "1000", "--nosimplify", "-p", prefix, fa_file, fa_file]
  subprocess.run(nucmer_cmd, stderr=subprocess.DEVNULL, check=True)

  coords_cmd = ["show-coords",  "-r", "-c", "-l", "-T", f"{prefix}.delta"]
  coords_outfile = f"{prefix}.coords"
  coords_data = []
  with open(coords_outfile, "w") as outfile:
      result = subprocess.run(coords_cmd, stdout=subprocess.PIPE, check=True, text=True)
      outfile.write(result.stdout)
      cols = result.stdout.splitlines()[3].split('\t')
      coords_data_vec = result.stdout.splitlines()[4:]
      for i, line in enumerate(coords_data_vec):
        coords_data.append({})
        for col, value in zip(cols, line.split()):
          col = col[1:-1]
          if col == 'TAGS': continue
          if '.' in value:
            coords_data[i][col] = float(value)
          else:
            coords_data[i][col] = int(value)

  # Plot alignment
  mummerplot_cmd = ["mummerplot", "--postscript", "-p", prefix, f"{prefix}.delta"]
  subprocess.run(mummerplot_cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, check=True)

  ps_file = f"{prefix}.ps"
  pdf_file = f"{prefix}.pdf"
  ps2pdf_cmd = ["ps2pdf", ps_file, pdf_file]
  subprocess.run(ps2pdf_cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, check=True)

  # Sanity check: first alignment must be match across the entire genome
  self_alignment = coords_data[0]
  if not (self_alignment['S1'] == 1 and self_alignment['E1'] == len(seq) and self_alignment['S2'] == 1 and self_alignment['E2'] == len(seq) and self_alignment['% IDY'] == 100.0):
    raise ValueError("First alignment is not a perfect match across the entire genome")
  return coords_data

def ir_near_edges(coords_data, seqlen):
  near_beginning_coord = -1
  near_end_coord = -1
  for alignment in coords_data[1:]:
    if (alignment['S2'] < alignment['E2']): continue
    min_coord = min(alignment['S1'], alignment['E1'], alignment['S2'], alignment['E2'])
    max_coord = max(alignment['S1'], alignment['E1'], alignment['S2'], alignment['E2'])
    if min_coord < 2000: near_beginning_coord = min_coord
    if max_coord > seqlen - 2000: near_end_coord = seqlen - max_coord
  return near_beginning_coord, near_end_coord

def rotate_seq(seq, n, direction):
  if direction == 'left':
    return seq[n:] + seq[:n]
  elif direction == 'right':
    return seq[len(seq) - n:] + seq[:len(seq) - n]
  else:
    raise ValueError(f"Invalid direction: {direction}")

def write_fa_file(fa_file, header, seq):
  with open(fa_file, "w") as outfile:
    outfile.write(f">{header}\n{seq}")

def identical_inversion_regions(align1, align2):
  return (
    align1.get('S1') == align2.get('E2') and
    align1.get('E1') == align2.get('S2') and
    align1.get('S2') == align2.get('E1') and
    align1.get('E2') == align2.get('S1')
  )

def write_regions_to_fastas(final_seq, header, align, prefix, log_fh):
  ir_a_start, ir_a_end = align['S1'], align['E1']
  ir_b_start, ir_b_end = align['E2'], align['S2']

  ir_first = final_seq[ir_a_start - 1:ir_a_end]
  ir_second = final_seq[ir_b_start - 1:ir_b_end]
  before = final_seq[:ir_a_start - 1]
  between = final_seq[ir_a_end:ir_b_start - 1]
  after = final_seq[ir_b_end:]

  split_region = after + before
  if len(split_region) >= len(between):
    lsc_seq = split_region
    ssc_seq = between
    irb_seq = ir_first
    ira_seq = ir_second
  else:
    lsc_seq = between
    ssc_seq = split_region
    irb_seq = ir_second
    ira_seq = ir_first

  write_fa_file(f'{prefix}.LSC.fa', f'{header}_LSC', lsc_seq)
  write_fa_file(f'{prefix}.IRb.fa', f'{header}_IRb', irb_seq)
  write_fa_file(f'{prefix}.SSC.fa', f'{header}_SSC', ssc_seq)
  write_fa_file(f'{prefix}.IRa.fa', f'{header}_IRa', ira_seq)




fa_file = sys.argv[1]
prefix = sys.argv[2]

dir_name = os.path.dirname(prefix)
if dir_name:
  os.makedirs(dir_name, exist_ok=True)


intermediates_dir = os.path.join(dir_name, prefix + "_intermediates")
os.makedirs(intermediates_dir, exist_ok=True)
intermediates_prefix = os.path.join(intermediates_dir, prefix)

record = next(SeqIO.parse(fa_file, "fasta"))
header = record.id
seq = str(record.seq)

# Run nucmer and get coords data
coords_data = run_nucmer(fa_file, seq, intermediates_prefix)

# slightly rotate the genome to find maximal IR regions
near_beginning_coord, near_end_coord = ir_near_edges(coords_data, len(seq))

rotated_seq = ''
rotated_dir = ''
if near_beginning_coord > -1 and near_end_coord > -1:
  rotated_dir = 'left'
elif near_beginning_coord > -1:
  rotated_dir = 'right'
elif near_end_coord > -1:
  rotated_dir = 'left'

rotated_seq = rotate_seq(seq, 40000, rotated_dir)


rotate_attempts = 0
rotated_fa_file = f'{intermediates_prefix}_rotated{rotate_attempts}.fa'
rotated_coords_data = None
if rotated_seq != '':
  write_fa_file(rotated_fa_file, header, rotated_seq)
  rotated_coords_data = run_nucmer(rotated_fa_file, rotated_seq, intermediates_prefix + f"_rotated{rotate_attempts}")
  near_beginning_coord_rotated, near_end_coord_rotated = ir_near_edges(rotated_coords_data, len(rotated_seq))
  if near_beginning_coord_rotated > -1 or near_end_coord_rotated > -1:
    rotate_attempts += 1
    rotated_seq = rotate_seq(rotated_seq, 20000, rotated_dir)
    rotated_fa_file = f'{intermediates_prefix}_rotated{rotate_attempts}.fa'
    write_fa_file(rotated_fa_file, header, rotated_seq)
    rotated_coords_data = run_nucmer(rotated_fa_file, rotated_seq, intermediates_prefix + f"_rotated{rotate_attempts}")
    near_beginning_coord_rotated, near_end_coord_rotated = ir_near_edges(rotated_coords_data, len(rotated_seq))
    if near_beginning_coord_rotated > -1 or near_end_coord_rotated > -1:
      print(f"Rotated genome has IR regions near the edges for {prefix} after {rotate_attempts} attempts, {fa_file}")

final_fa_file = f'{intermediates_prefix}_final.fa'
final_coords_file = f'{intermediates_prefix}_final.coords'
final_pdf_file = f'{intermediates_prefix}_final.pdf'
if rotated_seq != '':
  write_fa_file(final_fa_file, header, rotated_seq)
  shutil.copy(f'{intermediates_prefix}_rotated{rotate_attempts}.coords', final_coords_file)
  shutil.copy(f'{intermediates_prefix}_rotated{rotate_attempts}.pdf', final_pdf_file)
else:
  write_fa_file(final_fa_file, header, seq)
  shutil.copy(f'{intermediates_prefix}.coords', final_coords_file)
  shutil.copy(f'{intermediates_prefix}.pdf', final_pdf_file)

final_coords_data = coords_data if rotated_seq == '' else rotated_coords_data
final_coords_data = final_coords_data[1:]

log_file = f'{prefix}.log'
log_fh = open(log_file, "w")

if len(final_coords_data) == 2:
  if identical_inversion_regions(final_coords_data[0], final_coords_data[1]):
    final_seq = seq if rotated_seq == '' else rotated_seq
    write_regions_to_fastas(final_seq, header, final_coords_data[0], prefix, log_fh)
    log_fh.write(f"Single identical inversion regions pairs found in {fa_file} and written to {prefix}.*.fa\n")
  else:
    log_fh.write(f'Single inversion pairs but not identical found in {fa_file}\n')
    log_fh.write(f'Single inversion pairs but not identical found in {fa_file}: {final_coords_data[0]}\n')
    log_fh.write(f'Single inversion pairs but not identical found in {fa_file}: {final_coords_data[1]}\n')
else:
  log_fh.write(f'Found more than a single inversion pair in {fa_file}\n')
log_fh.close()


