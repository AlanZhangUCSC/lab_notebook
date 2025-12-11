import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

data_by_samples = {}
data_accumulated = {}

for file in os.listdir('.'):
  if not (file.startswith('mito_assignment_stats') and file.endswith('.tsv')):
    continue
  matches = re.findall(r'MPS0\.[0-9]+', file)
  mps = float(matches[0][3:])
  matches = re.findall(r'entropyk[0-9]+', file)
  entropy_k = int(matches[0][8:])
  if entropy_k != 5: continue
  cur_data_by_samples = []
  fh = open(file, 'r')
  TP_accumulated = 0
  FP_accumulated = 0
  FN_accumulated = 0
  for line in fh.readlines()[1:]:
    prim_mappable, panmama_filtered, intersection = list(map(int, line.strip().split()[2:5]))
    if prim_mappable == 0 and panmama_filtered == 0:
      continue
    TP = intersection
    FP = panmama_filtered - TP
    FN = prim_mappable - TP
    TP_accumulated += TP
    FP_accumulated += FP
    FN_accumulated += FN
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    F1 = 2 * TP / (2 * TP + FP + FN) if (2 * TP + FP + FN) > 0 else 0.0
    cur_data_by_samples.append([precision, recall, F1])
  fh.close()
  headers = ['Precision', 'Recall', 'F1']
  df = pd.DataFrame(cur_data_by_samples, columns=headers)
  data_by_samples[(mps, entropy_k)] = df
  precision_accumulated = TP_accumulated / (TP_accumulated + FP_accumulated) if (TP_accumulated + FP_accumulated) > 0 else 0.0
  recall_accumulated = TP_accumulated / (TP_accumulated + FN_accumulated) if (TP_accumulated + FN_accumulated) > 0 else 0.0
  F1_accumulated = 2 * TP_accumulated / (2 * TP_accumulated + FP_accumulated + FN_accumulated) if (2 * TP_accumulated + FP_accumulated + FN_accumulated) > 0 else 0.0
  data_accumulated[(mps, entropy_k)] = [precision_accumulated, recall_accumulated, F1_accumulated]

accumulated_data = list(data_accumulated.values())

for key in data_accumulated:
  print(f'MPS {key[0]} entropyk {key[1]}: {data_accumulated[key]}')
precisions = [point[0] for point in accumulated_data]
recalls = [point[1] for point in accumulated_data]

sorted_data = sorted(zip(recalls, precisions))
sorted_recalls, sorted_precisions = zip(*sorted_data)

plt.figure(figsize=(8, 6))
sns.set_style("whitegrid")
sns.lineplot(x=sorted_recalls, y=sorted_precisions, marker='o', markersize=8, linewidth=2, color="steelblue")
plt.xlabel('Recall', fontsize=12)
plt.ylabel('Precision', fontsize=12)
plt.title('Precision-Recall Curve (entropy_k=5)', fontsize=14)
plt.xlim([min(sorted_recalls) - 0.02, 1.05])
plt.ylim([-0.05, 1.05])
plt.tight_layout()

plt.savefig('mito_assignment_entropyk5_prc.png', dpi=600)