from collections import defaultdict
import sys
import os

dir_path = sys.argv[1]

sample_type = defaultdict(lambda: defaultdict(dict))
for filename in os.listdir(dir_path):
  if not filename.endswith('Scores.tsv'):  continue
  isPerfect = 'perfect' in filename
  isErrorDetect = 'err_detect' in filename

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
    if ('nodeScores' in filename):
      sample_type_string = filename.split('.')[0].replace('_perfect','')
      sample_type[sample_type_string]['readWeightsPerfect'] = [num_true_selected, num_true]
    elif ('nodeSeedScores' in filename):
      sample_type_string = filename.split('.')[0]
      sample_type[sample_type_string]['seedWeightsPerfect'] = [num_true_selected, num_true]
    elif ('nodeReadSeedScores' in filename):
      sample_type_string = filename.split('.')[0]
      sample_type[sample_type_string]['readSeedWeightsPerfect'] = [num_true_selected, num_true]
  elif (isErrorDetect):
    sample_type_string = filename.split('.')[0]
    if ('nodeScores' in filename):
      sample_type[sample_type_string]['readWeightsErrorDetect'] = [num_true_selected, num_true]
    elif ('nodeSeedScores' in filename):
      sample_type[sample_type_string]['seedWeightsErrorDetect'] = [num_true_selected, num_true]
    elif ('nodeReadSeedScores' in filename):
      sample_type[sample_type_string]['readSeedWeightsErrorDetect'] = [num_true_selected, num_true]
  else:
    sample_type_string = filename.split('.')[0]
    if ('nodeScores' in filename):
      sample_type[sample_type_string]['readWeights'] = [num_true_selected, num_true]
    elif ('nodeSeedScores' in filename):
      sample_type[sample_type_string]['seedWeights'] = [num_true_selected, num_true]
    elif ('nodeReadSeedScores' in filename):
      sample_type[sample_type_string]['readSeedWeights'] = [num_true_selected, num_true]

  
print("SampleType\tSeedWeightedScore\tReadWeightedScore\tReadSeedWeightedScore\tTotalTrue\tSeedWeightedScore_ErrorDetect\treadWeightedScore_ErrorDetect\tReadSeedWeightedScore_ErrorDetect\tTotalTrue_ErrorDetect\tSeedWeightedScore_Perfect\tReadWeightedScore_Perfect\tReadSeedWeightedScore_Perfect\tTotalTrue_Perfect")
for sample_type_str in sample_type.keys():
  if (not sample_type_str.startswith('amplicon')): 
    continue
  seedWeights = list(sample_type[sample_type_str]['seedWeights']) if 'seedWeights' in sample_type[sample_type_str] else [0,0]
  readWeights = list(sample_type[sample_type_str]['readWeights']) if 'readWeights' in sample_type[sample_type_str] else [0,0]
  readSeedWeights = list(sample_type[sample_type_str]['readSeedWeights']) if 'readSeedWeights' in sample_type[sample_type_str] else [0,0]
  assert(seedWeights[1] == readWeights[1])
  assert(seedWeights[1] == readSeedWeights[1])

  seedWeightsPerfect = list(sample_type[sample_type_str]['seedWeightsPerfect']) if 'seedWeightsPerfect' in sample_type[sample_type_str] else [0,0]
  readWeightsPerfect = list(sample_type[sample_type_str]['readWeightsPerfect']) if 'readWeightsPerfect' in sample_type[sample_type_str] else [0,0]
  readSeedWeightsPerfect = list(sample_type[sample_type_str]['readSeedWeightsPerfect']) if 'readSeedWeightsPerfect' in sample_type[sample_type_str] else [0,0]
  assert(seedWeightsPerfect[1] == readWeightsPerfect[1])
  assert(seedWeightsPerfect[1] == readSeedWeightsPerfect[1])

  seedWeightsErrorDetect = list(sample_type[sample_type_str]['seedWeightsErrorDetect']) if 'seedWeightsErrorDetect' in sample_type[sample_type_str] else [0,0]
  readWeightsErrorDetect = list(sample_type[sample_type_str]['readWeightsErrorDetect']) if 'readWeightsErrorDetect' in sample_type[sample_type_str] else [0,0]
  readSeedWeightsErrorDetect = list(sample_type[sample_type_str]['readSeedWeightsErrorDetect']) if 'readSeedWeightsErrorDetect' in sample_type[sample_type_str] else [0,0]
  assert(seedWeightsErrorDetect[1] == readWeightsErrorDetect[1])
  assert(seedWeightsErrorDetect[1] == readSeedWeightsErrorDetect[1])
  
  print(f"{sample_type_str}\t{seedWeights[0]}\t{readWeights[0]}\t{readSeedWeights[0]}\t{readWeights[1]}\t{seedWeightsErrorDetect[0]}\t{readWeightsErrorDetect[0]}\t{readSeedWeightsErrorDetect[0]}\t{seedWeightsErrorDetect[1]}\t{seedWeightsPerfect[0]}\t{readWeightsPerfect[0]}\t{readSeedWeightsPerfect[0]}\t{readWeightsPerfect[1]}")
  
