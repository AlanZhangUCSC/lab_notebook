import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

results_dir = sys.argv[1]

data = []
for directory in os.listdir(results_dir):
  cur_dir = os.path.join(results_dir, directory)
  if not os.path.isdir(cur_dir):
    continue
  panmap_out_files = [f for f in os.listdir(cur_dir) if f.endswith('.abundance.out')]
  for panmap_out_file in panmap_out_files:
    panmap_out_file_path = os.path.join(cur_dir, panmap_out_file)
    haplotype = panmap_out_file.split('_')[0]
    num_reads = int(panmap_out_file.split('_')[1])
    damaged = 'nodmg' not in panmap_out_file
    maskends = 'maskends' in panmap_out_file
    estimated_proportion = 0.0
    is_top_haplotype = False
    with open(panmap_out_file_path, 'r') as f:
      for index, line in enumerate(f):
        fields = line.strip().split('\t')
        if haplotype in fields[0]:
          if index == 0:
            is_top_haplotype = True
          estimated_proportion = float(fields[1])
          break

    weighted_peak_distance = None
    if damaged:
      weighted_peak_distance = 0.0
      with open(panmap_out_file_path.replace('.abundance.out', '.abundance.with_distance.out'), 'r') as f:
        for line in f:
          fields = line.strip().split('\t')
          abundance = float(fields[1])
          distance = float(fields[2])
          weighted_peak_distance += abundance * distance
    data.append([
      num_reads,
      damaged,
      maskends,
      haplotype,
      estimated_proportion,
      is_top_haplotype,
      weighted_peak_distance,
    ])

df = pd.DataFrame(
  data,
  columns=[
    'num_reads',
    'damaged',
    'maskends',
    'haplotype',
    'estimated_proportion',
    'is_top_haplotype',
    'weighted_peak_distance',
  ],
)
# None for undamaged rows makes the column object dtype; coerce so the KDE gets floats.
df['weighted_peak_distance'] = pd.to_numeric(df['weighted_peak_distance'], errors='coerce')

UNDAMAGED = 'Undamaged'
DAMAGED = 'Damaged'
DAMAGED_MASKED = 'Damaged + masked ends'

def label_condition(row):
  if not row['damaged']:
    return UNDAMAGED
  return DAMAGED_MASKED if row['maskends'] else DAMAGED

# maskends is only meaningful for damaged reads; flag any file that violates this
# so a naming typo doesn't silently vanish into the Undamaged group.
impossible = df[(~df['damaged']) & df['maskends']]
if not impossible.empty:
  print(
    'Warning: {} undamaged file(s) also matched "maskends"'.format(len(impossible)),
    file=sys.stderr,
  )

df['condition'] = df.apply(label_condition, axis=1)

read_order = sorted(df['num_reads'].unique())
condition_order = [DAMAGED, DAMAGED_MASKED, UNDAMAGED]
condition_order = [c for c in condition_order if c in set(df['condition'])]
palette = dict(zip([UNDAMAGED, DAMAGED, DAMAGED_MASKED], sns.color_palette('Set2', 3)))

sns.set_theme(style='whitegrid')

def violin_with_points(data, y, hue_order, ax):
  sns.violinplot(
    data=data,
    x='num_reads',
    y=y,
    hue='condition',
    order=read_order,
    hue_order=hue_order,
    palette=palette,
    inner=None,
    cut=0,
    linewidth=1,
    ax=ax,
  )
  # Capture before the strip layer duplicates the hue entries, and lighten only
  # the violin bodies.
  handles, labels = ax.get_legend_handles_labels()
  for body in ax.collections:
    body.set_alpha(0.3)
  sns.stripplot(
    data=data,
    x='num_reads',
    y=y,
    hue='condition',
    order=read_order,
    hue_order=hue_order,
    palette=palette,
    dodge=True,
    jitter=0.12,
    size=3,
    alpha=0.75,
    linewidth=0,
    ax=ax,
  )
  ax.legend(handles, labels, title='', loc='best')

fig, ax = plt.subplots(figsize=(max(9, len(read_order) * 2.2), 6))
violin_with_points(df, 'estimated_proportion', condition_order, ax)
ax.set_xlabel('Number of reads')
ax.set_ylabel('Estimated proportion')
ax.set_title('Estimated proportion by read count and damage treatment')
fig.tight_layout()
fig.savefig('proportion_violins.png', dpi=200)

top_fraction = (
  df.groupby(['num_reads', 'condition'])['is_top_haplotype'].agg(['mean', 'size']).reset_index()
)
top_fraction = top_fraction.rename(columns={'mean': 'fraction_top', 'size': 'n'})

fig2, ax2 = plt.subplots(figsize=(8, 6))
sns.lineplot(
  data=top_fraction,
  x='num_reads',
  y='fraction_top',
  hue='condition',
  hue_order=condition_order,
  palette=palette,
  marker='o',
  ax=ax2,
)
ax2.set_xscale('log')
ax2.set_ylim(0, 1.05)
ax2.set_xlabel('Number of reads')
ax2.set_ylabel('Fraction of replicates with true haplotype ranked top')
ax2.set_title('Top-haplotype recovery rate by read count')
ax2.legend(title='', loc='best')
fig2.tight_layout()
fig2.savefig('top_haplotype_fraction.png', dpi=200)

damaged_df = df[df['damaged'] & df['weighted_peak_distance'].notna()]
damaged_order = [c for c in [DAMAGED, DAMAGED_MASKED] if c in set(damaged_df['condition'])]

if damaged_df.empty:
  print('No damaged samples with a weighted peak distance; skipping plot 3.', file=sys.stderr)
else:
  fig3, ax3 = plt.subplots(figsize=(max(9, len(read_order) * 1.8), 6))
  violin_with_points(damaged_df, 'weighted_peak_distance', damaged_order, ax3)
  ax3.set_ylim(bottom=0)
  ax3.set_xlabel('Number of reads')
  ax3.set_ylabel('Weighted peak distance')
  ax3.set_title('Weighted peak distance by read count (damaged samples)')
  fig3.tight_layout()
  fig3.savefig('weighted_peak_distance_violins.png', dpi=200)

print(top_fraction.to_string(index=False))
print(
  damaged_df.groupby(['num_reads', 'condition'])['weighted_peak_distance']
  .agg(['median', 'mean', 'size'])
  .to_string()
)