import sys
import os
import uuid
import pandas as pd
import subprocess
from ete3 import Tree


DIR = sys.argv[1]
out = sys.argv[2]

tree = Tree('input_data/v_mtdna.new.panman.nwk', format=1)


headers = ['sample',
           'selection_scheme',
           'primigenius_mappable',
           'post_filter_reads',
           'post_filter_reads_union_mappable',
           'assigned_to_primigenius',
           'assigned_to_primigenius_mappable',
           'assigned_to_mammuthus',
           'assigned_to_mammuthus_mappable',
           'assigned_to_elephantidae',
           'assigned_to_elephantidae_mappable',
           'assigned_elsewhere',
           'assigned_elsewhere_mappable']


primigeniuses = ['DQ188829.2', 'NC_007596.2']
mammoths = ['DQ188829.2', 'NC_007596.2', 'NC_015529.1', 'JF912199.1']
elephantidae = ['node_14398']

elephantidae_root_node = tree.search_nodes(name='node_14398')[0]
for node in elephantidae_root_node.get_descendants():
  elephantidae.append(node.name)

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
  post_filter_reads = 0
  assigned_to_primigenius = 0
  assigned_to_mammuthus = 0
  assigned_elsewhere = 0
  assigned_to_elephantidae = 0
  assigned_to_primigenius_mappable_union = 0
  assigned_to_mammuthus_mappable_union = 0
  assigned_elsewhere_mappable_union = 0
  assigned_to_elephantidae_mappable_union = 0
  
  read_scores_info = f'{prefix}.{selection_scheme}.read_scores_info.tsv'
  if (os.path.exists(read_scores_info)):
    fh = open(read_scores_info, 'r')
    fh.readline()
    num_reads = 0
    for line in fh.readlines():
      num_reads += len(line.strip().split('\t')[-1].split(','))
    fh.close()
    post_filter_reads = num_reads
  
  woolly_mammoth_bam = f'final_analysis/data/reads/{base_name}.mapped.bam'
  mappedReadIds = set()
  if (os.path.exists(woolly_mammoth_bam)):
    cmd = f'samtools view {woolly_mammoth_bam} | cut -f 1'
    mappedReadIds = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    if mappedReadIds != '':
      mappedReadIds = set(mappedReadIds.split('\n'))
  
  assigned_reads = f'{prefix}.{selection_scheme}.mgsr.assignedReads.out'

  read_ids = []
  readIdtoIndex = {}
  if (os.path.exists(assigned_reads)):
    uuid_string = str(uuid.uuid4())
    fastq_path = f'final_analysis/data/reads/{base_name}.fq'
    fastq_filtered = f'/data/tmp/{uuid_string}.{base_name}.filtered.fq'
    subprocess.run(f'~/tools/BBTools/bbduk.sh in={fastq_path} out={fastq_filtered} entropy=0.7', shell=True, check=True)

    fastq_filtered_sorted = f'/data/tmp/{uuid_string}.{base_name}.filtered.sorted.fq'
    subprocess.run(f"cat {fastq_filtered} | paste - - - - | sort -k1,1 -S 3G | tr '\\t' '\\n' > {fastq_filtered_sorted}", shell=True, check=True)

    with open(fastq_filtered_sorted, 'r') as fh:
      for i, line in enumerate(fh.readlines()):
        if line.strip() == '': break
        if i % 4 == 0:
          read_id = line.strip().split()[0][1:]
          read_ids.append(read_id)
          readIdtoIndex[read_id] = len(read_ids) - 1
    
    os.remove(fastq_filtered)
    os.remove(fastq_filtered_sorted)
  
  mappedReadIndices = set()
  pseudo = -1
  for read_id in mappedReadIds:
    if read_id in readIdtoIndex:
      mappedReadIndices.add(readIdtoIndex[read_id])
    else:
      mappedReadIndices.add(pseudo)
      pseudo -= 1
    
  if (os.path.exists(assigned_reads)):
    abundance = f'{prefix}.{selection_scheme}.mgsr.abundance.out'
    if not (os.path.exists(abundance)):
      print(f'Abundance file missing for {assigned_reads}', file=sys.stderr)
      exit(1)
    primigeniuses_assigned_line = []
    mammuthus_assigned_line = []
    elephantidae_assigned_line = []
    elsewhere_assigned_line = []
    with open(abundance, 'r') as fh:
      for i, line in enumerate(fh.readlines()):
        haplotypes = line.strip().split()[0]
        assigned_elephantidae = False
        for primigenius in primigeniuses:
          if primigenius in haplotypes:
            primigeniuses_assigned_line.append(i)
            break
        for mammoth in mammoths:
          if mammoth in haplotypes:
            mammuthus_assigned_line.append(i)
            break
        for elephant in elephantidae:
          if elephant in haplotypes:
            assigned_elephantidae = True
            elephantidae_assigned_line.append(i)
            break
        if not assigned_elephantidae:
          elsewhere_assigned_line.append(i)
        
    assigned_primigenius_reads_set = set()
    assigned_mammuthus_reads_set = set()
    assigned_elephantidae_reads_set = set()
    assigned_elsewhere_reads_set = set()
    with open(assigned_reads) as fh:
      lines = fh.readlines()
      for index in primigeniuses_assigned_line:
        if int(lines[index].strip().split()[1]) > 0:
          assigned_primigenius_reads_set = assigned_primigenius_reads_set.union(set(map(int, lines[index].strip().split()[-1].split(','))))
      for index in mammuthus_assigned_line:
        if int(lines[index].strip().split()[1]) > 0:
          assigned_mammuthus_reads_set = assigned_mammuthus_reads_set.union(set(map(int, lines[index].strip().split()[-1].split(','))))
      for index in elephantidae_assigned_line:
        if int(lines[index].strip().split()[1]) > 0:
          assigned_elephantidae_reads_set = assigned_elephantidae_reads_set.union(set(map(int, lines[index].strip().split()[-1].split(','))))
      for index in elsewhere_assigned_line:
        if int(lines[index].strip().split()[1]) > 0:
          assigned_elsewhere_reads_set = assigned_elsewhere_reads_set.union(set(map(int, lines[index].strip().split()[-1].split(','))))
    
    assigned_to_primigenius = len(assigned_primigenius_reads_set)
    assigned_to_primigenius_mappable_union = len(assigned_primigenius_reads_set.intersection(mappedReadIndices))
    assigned_to_mammuthus = len(assigned_mammuthus_reads_set)
    assigned_to_mammuthus_mappable_union = len(assigned_mammuthus_reads_set.intersection(mappedReadIndices))
    assigned_to_elephantidae = len(assigned_elephantidae_reads_set)
    assigned_to_elephantidae_mappable_union = len(assigned_elephantidae_reads_set.intersection(mappedReadIndices))
    assigned_elsewhere = len(assigned_elsewhere_reads_set)
    assigned_elsewhere_mappable_union = len(assigned_elsewhere_reads_set.intersection(mappedReadIndices))
  

  data.append([sample,
               selection_scheme,
               len(mappedReadIndices),
               post_filter_reads,
               assigned_to_elephantidae_mappable_union + assigned_elsewhere_mappable_union,
               assigned_to_primigenius,
               assigned_to_primigenius_mappable_union,
               assigned_to_mammuthus,
               assigned_to_mammuthus_mappable_union,
               assigned_to_elephantidae,
               assigned_to_elephantidae_mappable_union,
               assigned_elsewhere,
               assigned_elsewhere_mappable_union])


df = pd.DataFrame(data, columns=headers)

# # replace -1 with NA
# df = df.replace({
#   'post_filter_reads': -1,
#   'assigned_to_mammuthus': -1,
#   'assigned_to_primigenius': -1,
#   'assigned_to_elephantidae': -1,
#   'assigned_elsewhere': -1
# }, pd.NA)

# sort by sample and selection_scheme
df = df.sort_values(by=['sample', 'selection_scheme'], ascending=[True, False])

# drop rows with NA in any column except sample and selection_scheme
df = df.dropna(subset=['post_filter_reads'])

# replace na with .
df = df.fillna('.')

df.to_csv(out, sep='\t', index=False)
