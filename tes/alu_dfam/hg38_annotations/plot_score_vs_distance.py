import sys
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

x_coords = {}
with open(sys.argv[1], 'r') as f:
  for line in f:
    data = json.loads(line)
    if "version" in data: continue
    if data['name'] == '': continue
    x_coords[data['name']] = data['x_dist']
print(f'Found {len(x_coords)} x_coords')

metadata = {}
with open(sys.argv[2], 'r') as f:
  for line in f:
    sample, target, model_acc, model_name, bitscore, evalue, hmm_start, hmm_end, hmm_len, strand, aln_start, aln_end, env_start, env_end, target_length = line.strip().split('\t')
    if sample not in x_coords: continue
    metadata[sample] = {
      'seqlen': int(hmm_end) - int(hmm_start),
      'bitscore': float(bitscore),
      'evalue': float(evalue)
    }
print(f'Loaded {len(metadata)} metadata')

df = pd.DataFrame([
  {'seqid': s, 'x': x_coords[s], **metadata[s]}
  for s in metadata
])

df['neg_log_evalue'] = -np.log10(df['evalue'].clip(lower=1e-300))

def fit_label(x, y):
  slope, intercept, r, p, _ = stats.linregress(x, y)
  return f'y = {slope:.3g}x + {intercept:.3g}\nR² = {r**2:.3f}, p = {p:.2e}'

def annotate(ax, x, y):
  ax.text(0.05, 0.95, fit_label(x, y),
          transform=ax.transAxes, va='top',
          bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.8})

def panel(ax, xcol, ycol, xlabel, ylabel, title):
  sns.regplot(
    data=df, x=xcol, y=ycol,
    scatter_kws={'s': 20, 'alpha': 0.6},
    line_kws={'color': 'red'},
    ax=ax,
  )
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_title(title)
  annotate(ax, df[xcol], df[ycol])

sns.set_theme(style='whitegrid')
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
ax1, ax2, ax3 = axes

panel(ax1, 'bitscore', 'x', 'Bit score', 'x coordinate', 'Bit score vs x_coords')
panel(ax2, 'seqlen', 'x', 'Sequence length', 'x coordinate', 'Alignment length vs x_coords')
panel(ax3, 'seqlen', 'bitscore', 'Sequence length', 'Bit score', 'Alignment length vs bit score')

plt.tight_layout()
plt.savefig('dotplots.png', dpi=150)
plt.show()