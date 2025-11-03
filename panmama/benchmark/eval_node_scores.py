from collections import defaultdict
import sys
import os

dir_path = sys.argv[1]

sample_type = defaultdict(lambda: defaultdict(dict))
for filename in os.listdir(dir_path):
  if not filename.endswith('Scores.tsv'):  continue

  file_path = os.path.join(dir_path, filename)

  fh = open(file_path)
  lines = fh.readlines()
  fh.close()

  num_true = 0
  num_true_selected = 0
  for line in lines[1:]:
    fields = line.strip().split('\t')
    if fields[6] == '1':
      num_true += 1
      if fields[5] == '1':
        num_true_selected += 1
  
  sample_type_string = filename.split('.')[0]
  if (filename.split('.')[1] == 'seedWeights'):
    sample_type[sample_type_string]['seedWeights'] = [num_true_selected, num_true]
  elif (filename.split('.')[1] == 'nodeScores'):
    sample_type[sample_type_string]['readWeights'] = [num_true_selected, num_true]
  
print("SampleType\tSeedWeightedScore\tReadWeightedScore\tTotalTrue")
for sample_type_str in sample_type.keys():
  seedWeights = list(sample_type[sample_type_str]['seedWeights']) if 'seedWeights' in sample_type[sample_type_str] else [0,0]
  readWeights = list(sample_type[sample_type_str]['readWeights']) if 'readWeights' in sample_type[sample_type_str] else [0,0]
  print(f"{sample_type_str}\t{seedWeights[0]}\t{readWeights[0]}\t{readWeights[1]}")
  
