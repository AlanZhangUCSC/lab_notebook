from collections import defaultdict
import sys
import os

dir_path = sys.argv[1]

sample_type = defaultdict(lambda: defaultdict(dict))
for filename in os.listdir(dir_path):
  if not filename.endswith('Scores.tsv'):  continue
  isPerfect = 'perfect' in filename

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
  
  if (isPerfect):
    if (filename.split('.')[1] == 'nodeScores'):
      sample_type_string = filename.split('.')[0].replace('_perfect','')
      sample_type[sample_type_string]['readWeightsPerfect'] = [num_true_selected, num_true]
    elif (filename.split('.')[1] == 'seedWeights_perfect'):
      sample_type_string = filename.split('.')[0]
      sample_type[sample_type_string]['seedWeightsPerfect'] = [num_true_selected, num_true]
  else:
    sample_type_string = filename.split('.')[0]
    if (filename.split('.')[1] == 'seedWeights'):
      sample_type[sample_type_string]['seedWeights'] = [num_true_selected, num_true]
    elif (filename.split('.')[1] == 'nodeScores'):
      sample_type[sample_type_string]['readWeights'] = [num_true_selected, num_true]
  
print("SampleType\tSeedWeightedScore\tReadWeightedScore\tTotalTrue\tSeedWeightedScore_Perfect\tReadWeightedScore_Perfect\tTotalTrue_Perfect")
for sample_type_str in sample_type.keys():
  seedWeights = list(sample_type[sample_type_str]['seedWeights']) if 'seedWeights' in sample_type[sample_type_str] else [0,0]
  readWeights = list(sample_type[sample_type_str]['readWeights']) if 'readWeights' in sample_type[sample_type_str] else [0,0]
  assert(seedWeights[1] == readWeights[1])

  seedWeightsPerfect = list(sample_type[sample_type_str]['seedWeightsPerfect']) if 'seedWeightsPerfect' in sample_type[sample_type_str] else [0,0]
  readWeightsPerfect = list(sample_type[sample_type_str]['readWeightsPerfect']) if 'readWeightsPerfect' in sample_type[sample_type_str] else [0,0]
  assert(seedWeightsPerfect[1] == readWeightsPerfect[1])
  
  print(f"{sample_type_str}\t{seedWeights[0]}\t{readWeights[0]}\t{readWeights[1]}\t{seedWeightsPerfect[0]}\t{readWeightsPerfect[0]}\t{readWeightsPerfect[1]}")
  
