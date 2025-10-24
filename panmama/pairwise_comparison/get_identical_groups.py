import sys
from collections import defaultdict

def find_identical_groups(filename):
  identical_map = defaultdict(list)
  
  with open(filename, 'r') as f:
    next(f)
    
    for line in f:
      fields = line.strip().split('\t')
      if (len(fields) < 6):
        continue
      seq1, seq2 = fields[0], fields[1]
      snps, snps_amb, gaps, gap_edge = map(int, fields[2:6])
      
      if snps == 0 and snps_amb == 0 and gaps == 0 and gap_edge == 0:
        key = tuple(sorted([seq1, seq2]))
        identical_map[key].append((seq1, seq2))
  
  groups = defaultdict(set)
  visited = set()
  
  for seq1, seq2 in identical_map.keys():
    if seq1 not in visited or seq2 not in visited:
      group = set()
      stack = [seq1, seq2]
      
      while stack:
        node = stack.pop()
        if node in visited:
          continue
        visited.add(node)
        group.add(node)
        
        for key in identical_map.keys():
          if node in key:
            other = key[1] if key[0] == node else key[0]
            if other not in visited:
              stack.append(other)
      
      if len(group) > 1:
        groups[frozenset(group)] = group
  
  sorted_groups = sorted(groups.values(), key=len, reverse=True)
  
  print('seqID\tidentical_group_id')
  identical_group_id = 0
  for idx, group in enumerate(sorted_groups):
    # skip over groups that do not have at least 2 leaf nodes
    num_leaf = 0
    for seq in group:
      if not seq.startswith('node_'):
        num_leaf += 1
    if num_leaf < 2: continue
    
    for seq in sorted(group):
      print(f"{seq}\t#{identical_group_id}#")
    identical_group_id+=1

find_identical_groups(sys.argv[1])
