import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data_by_samples = {}
data_accumulated = {}

for file in os.listdir('.'):
  if not (file.startswith('mito_assignment_stats') and file.endswith('.tsv')):
    continue
  matches = re.findall(r'MPS0\.[0-9]+', file)
  mps = float(matches[0][3:])
  matches = re.findall(r'entropyk[0-9]+', file)
  entropy_k = int(matches[0][8:])
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

for key in data_accumulated:
  print(f'MPS {key[0]} entropyk {key[1]}: {data_accumulated[key]}')

mps_values = sorted(set(key[0] for key in data_accumulated.keys()))
entropy_k_values = sorted(set(key[1] for key in data_accumulated.keys()))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
colors = plt.cm.tab10(np.linspace(0, 1, len(entropy_k_values)))

for entropy_idx, entropy_k in enumerate(entropy_k_values):
  precision_values = []
  recall_values = []
  f1_values = []
  for mps in mps_values:
    key = (mps, entropy_k)
    if key in data_accumulated:
      precision_values.append(data_accumulated[key][0])
      recall_values.append(data_accumulated[key][1])
      f1_values.append(data_accumulated[key][2])
    else:
      precision_values.append(np.nan)
      recall_values.append(np.nan)
      f1_values.append(np.nan)
  
  axes[0].plot(mps_values, precision_values, marker='o', color=colors[entropy_idx], label=f'k={entropy_k}', linestyle='-')
  axes[0].plot(mps_values, recall_values, marker='s', color=colors[entropy_idx], linestyle='--')
  axes[1].plot(mps_values, f1_values, marker='o', color=colors[entropy_idx], label=f'k={entropy_k}')

axes[0].plot([], [], marker='o', color='black', linestyle='-', label='Precision')
axes[0].plot([], [], marker='s', color='black', linestyle='--', label='Recall')
axes[0].set_xticks(mps_values)
axes[0].set_xlabel('MPS')
axes[0].set_ylabel('Score')
axes[0].set_title('Precision and Recall')
axes[0].legend(loc='best')
axes[0].grid(True, alpha=0.3)
axes[0].set_ylim([0, 1])

axes[1].set_xticks(mps_values)
axes[1].set_xlabel('MPS')
axes[1].set_ylabel('F1 Score')
axes[1].set_title('F1 Score')
axes[1].legend()
axes[1].grid(True, alpha=0.3)
axes[1].set_ylim([0, 1])

plt.tight_layout()
plt.savefig('mito_assignment_stats_comparison.png', dpi=300, bbox_inches='tight')
