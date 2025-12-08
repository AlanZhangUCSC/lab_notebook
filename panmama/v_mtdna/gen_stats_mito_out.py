import sys
import os
import pandas as pd

headers = ['sample', 'selection_scheme', 'post_filter_reads', 'assigned_to_mammuthus', 'assigned_elsewhere']

DIR = sys.argv[1]

mammoths = ['DQ188829.2', 'NC_007596.2', 'NC_015529.1', 'JF912199.1']
processed = set()
data = []
for file in os.listdir(DIR):
  base_name = '.'.join(file.split('.')[:2])
  prefix = os.path.join(DIR, base_name)
  sample = base_name.split('.')[0]
  selection_scheme = file.split('.')[2]
  
  if (sample, selection_scheme) in processed:
    continue
  processed.add((sample, selection_scheme))
  post_filter_reads = -1
  assigned_to_mammuthus = -1
  assigned_elsewhere = -1
  
  read_scores_info = f'{prefix}.{selection_scheme}.read_scores_info.tsv'
  if (os.path.exists(read_scores_info)):
    fh = open(read_scores_info, 'r')
    fh.readline()
    num_reads = 0
    for line in fh.readlines():
      num_reads += len(line.strip().split('\t')[-1].split(','))
    fh.close()
    post_filter_reads = num_reads
  
  assigned_reads = f'{prefix}.{selection_scheme}.mgsr.assignedReads.out'
  if (os.path.exists(assigned_reads)):
    abundance = f'{prefix}.{selection_scheme}.mgsr.abundance.out'
    if not (os.path.exists(abundance)):
      print(f'Abundance file missing for {assigned_reads}', file=sys.stderr)
      exit(1)
    mammuthus_assigned_line = []
    with open(abundance, 'r') as fh:
      for i, line in enumerate(fh.readlines()):
        haplotypes = line.strip().split()[0]
        for mammoth in mammoths:
          if mammoth in haplotypes:
            mammuthus_assigned_line.append(i)
            break
    assigned_reads_set = set()
    with open(assigned_reads) as fh:
      lines = fh.readlines()
      for index in mammuthus_assigned_line:
        if int(lines[index].strip().split()[1]) > 0:
          assigned_reads_set = assigned_reads_set.union(set(lines[index].strip().split()[-1].split(',')))  
    assigned_to_mammuthus = len(assigned_reads_set)
  
  if post_filter_reads != -1 and assigned_to_mammuthus != -1:
    assigned_elsewhere = post_filter_reads - assigned_to_mammuthus
  
  data.append([sample, selection_scheme, post_filter_reads, assigned_to_mammuthus, assigned_elsewhere])

df = pd.DataFrame(data, columns=headers)

# replace -1 with NA
df = df.replace({
  'post_filter_reads': -1,
  'assigned_to_mammuthus': -1,
  'assigned_elsewhere': -1
}, pd.NA)

# sort by sample and selection_scheme
df = df.sort_values(by=['sample', 'selection_scheme'], ascending=[True, False])

# drop rows with NA in any column except sample and selection_scheme
df = df.dropna(subset=['post_filter_reads'])

# replace na with .
df = df.fillna('.')

df.to_csv('mito_assignment_stats.tsv', sep='\t', index=False)

