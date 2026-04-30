import pysam
import re
import sys

MD_RE = re.compile(r'(\d+)|(\^[A-Za-z]+)|([A-Za-z])')

def alignment_stats(read):
  if read.is_unmapped or read.is_secondary or read.is_supplementary:
    return None
  md = read.get_tag('MD') if read.has_tag('MD') else None
  if md is None:
    return None

  matches = mismatches = 0
  for token_num, token_del, token_sub in MD_RE.findall(md):
    if token_num:
      matches += int(token_num)
    elif token_sub:
      mismatches += 1

  insertions = deletions = soft_clip = hard_clip = 0
  for op, length in read.cigartuples:
    if op == 1:
      insertions += length
    elif op == 2:
      deletions += length
    elif op == 4:
      soft_clip += length
    elif op == 5:
      hard_clip += length

  aligned_bases = matches + mismatches + insertions + deletions
  if aligned_bases == 0:
    return None

  identity_A = matches / aligned_bases
  identity_B = matches / (aligned_bases + soft_clip + hard_clip)

  return {
    'name': read.query_name,
    'matches': matches,
    'mismatches': mismatches,
    'ins': insertions,
    'del': deletions,
    'soft': soft_clip,
    'hard': hard_clip,
    'aln_len': aligned_bases,
    'query_len': aligned_bases - deletions + soft_clip + hard_clip,
    'id_A': identity_A,
    'id_B': identity_B,
  }

def process(bam_path, out_path, min_mapq=20):
  cols = ['name', 'matches', 'mismatches', 'ins', 'del',
          'soft', 'hard', 'aln_len', 'query_len', 'id_A', 'id_B']
  with pysam.AlignmentFile(bam_path, 'rb') as bam, open(out_path, 'w') as out:
    out.write('\t'.join(cols) + '\n')
    for read in bam:
      if read.mapping_quality < min_mapq:
        continue
      r = alignment_stats(read)
      if r is None:
        continue
      out.write('\t'.join(str(r[c]) for c in cols) + '\n')

if __name__ == "__main__":
  process(sys.argv[1], sys.argv[2])