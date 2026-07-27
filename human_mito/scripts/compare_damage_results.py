import sys
import os
import bisect
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

results_dir = sys.argv[1]

TOLERANCES = [0, 1, 3, 5, 7, 9]
PEAK_TOLERANCE_COLUMNS = ['peak_tolerance_{}'.format(t) for t in TOLERANCES]
GENOME_TOLERANCE_COLUMNS = ['genome_tolerance_{}'.format(t) for t in TOLERANCES]

def _cumulative(buckets):
  running = 0.0
  out = []
  for i in range(len(TOLERANCES)):
    running += buckets[i]
    out.append(running)
  return out

def proportions_by_tolerance(distance_path):
  # peak/genome buckets[i] collect abundance whose distance first qualifies at
  # TOLERANCES[i]; the trailing slot absorbs nodes beyond the largest tolerance.
  peak_buckets = [0.0] * (len(TOLERANCES) + 1)
  genome_buckets = [0.0] * (len(TOLERANCES) + 1)
  weighted_peak_distance = 0.0
  with open(distance_path, 'r') as f:
    for line in f:
      fields = line.strip().split('\t')
      if len(fields) < 4:
        continue
      abundance = float(fields[1])
      peak_distance = float(fields[2])
      genome_distance = float(fields[3])
      weighted_peak_distance += abundance * peak_distance
      peak_buckets[bisect.bisect_left(TOLERANCES, peak_distance)] += abundance
      genome_buckets[bisect.bisect_left(TOLERANCES, genome_distance)] += abundance
  return _cumulative(peak_buckets), _cumulative(genome_buckets), weighted_peak_distance

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

    distance_path = panmap_out_file_path.replace('.abundance.out', '.abundance.with_distance.out')
    assert os.path.exists(distance_path), 'missing distance file: {}'.format(distance_path)
    peak_proportions, genome_proportions, weighted_peak_distance = proportions_by_tolerance(distance_path)
    data.append([
      num_reads,
      damaged,
      maskends,
      haplotype,
      estimated_proportion,
      is_top_haplotype,
      weighted_peak_distance,
    ] + peak_proportions + genome_proportions)

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
  ] + PEAK_TOLERANCE_COLUMNS + GENOME_TOLERANCE_COLUMNS,
)

UNDAMAGED = 'Undamaged'
DAMAGED = 'Damaged'
DAMAGED_MASKED = 'Damaged + masked ends'

def label_condition(row):
  if not row['damaged']:
    return UNDAMAGED
  return DAMAGED_MASKED if row['maskends'] else DAMAGED

impossible = df[(~df['damaged']) & df['maskends']]
if not impossible.empty:
  print(
    'Warning: {} undamaged file(s) also matched "maskends"'.format(len(impossible)),
    file=sys.stderr,
  )

df['condition'] = df.apply(label_condition, axis=1)

# peak tolerance 0 sums every node at peak-distance 0, which need not equal the
# single name-matched row; a large gap means the two are not measuring the same thing.
tolerance_zero_gap = (df['peak_tolerance_0'] - df['estimated_proportion']).abs()
if tolerance_zero_gap.max() > 1e-9:
  print(
    'Note: tolerance 0 differs from estimated_proportion (max gap {:.4g}, {} of {} rows)'.format(
      tolerance_zero_gap.max(), int((tolerance_zero_gap > 1e-9).sum()), len(df)
    ),
    file=sys.stderr,
  )

read_order = sorted(df['num_reads'].unique())
condition_order = [DAMAGED, DAMAGED_MASKED, UNDAMAGED]
condition_order = [c for c in condition_order if c in set(df['condition'])]
palette = dict(zip([UNDAMAGED, DAMAGED, DAMAGED_MASKED], sns.color_palette('Set2', 3)))

tolerance_order = [str(t) for t in TOLERANCES]
tolerance_palette = dict(zip(tolerance_order, sns.color_palette('viridis', len(TOLERANCES))))

sns.set_theme(style='whitegrid')

def violin_with_points(data, y, hue, hue_order, hue_palette, ax, jitter=0.12, size=3, legend=True):
  sns.violinplot(
    data=data,
    x='num_reads',
    y=y,
    hue=hue,
    order=read_order,
    hue_order=hue_order,
    palette=hue_palette,
    inner=None,
    cut=0,
    density_norm='width',
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
    hue=hue,
    order=read_order,
    hue_order=hue_order,
    palette=hue_palette,
    dodge=True,
    jitter=jitter,
    size=size,
    alpha=0.75,
    linewidth=0,
    ax=ax,
  )
  if legend:
    ax.legend(handles, labels, title='', loc='best')
  elif ax.get_legend() is not None:
    ax.get_legend().remove()
  return handles, labels

fig, ax = plt.subplots(figsize=(max(9, len(read_order) * 2.2), 6))
violin_with_points(df, 'estimated_proportion', 'condition', condition_order, palette, ax)
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

fig3, ax3 = plt.subplots(figsize=(max(9, len(read_order) * 2.2), 6))
violin_with_points(df, 'weighted_peak_distance', 'condition', condition_order, palette, ax3)
ax3.set_ylim(bottom=0)
ax3.set_xlabel('Number of reads')
ax3.set_ylabel('Weighted peak distance')
ax3.set_title('Weighted peak distance by read count and damage treatment')
fig3.tight_layout()
fig3.savefig('weighted_peak_distance_violins.png', dpi=200)

def tolerance_panels(tolerance_columns, prefix, title, outfile):
  long_df = df.melt(
    id_vars=['num_reads', 'condition', 'haplotype'],
    value_vars=tolerance_columns,
    var_name='tolerance',
    value_name='proportion',
  )
  long_df['tolerance'] = long_df['tolerance'].str.replace(prefix, '', regex=False)

  fig, axes = plt.subplots(
    len(condition_order),
    1,
    figsize=(max(9.0, len(read_order) * 2.2), 4.5 * len(condition_order)),
    sharey=True,
    sharex=True,
  )
  axes = [axes] if len(condition_order) == 1 else list(axes)

  handles = labels = None
  for panel_ax, condition in zip(axes, condition_order):
    subset = long_df[long_df['condition'] == condition]
    handles, labels = violin_with_points(
      subset,
      'proportion',
      'tolerance',
      tolerance_order,
      tolerance_palette,
      panel_ax,
      jitter=0.08,
      size=2.0,
      legend=False,
    )
    panel_ax.set_title(condition)
    panel_ax.set_xlabel('')
    panel_ax.set_ylabel('Estimated proportion')
    panel_ax.set_ylim(bottom=0)

  axes[-1].set_xlabel('Number of reads')
  fig.legend(
    handles,
    labels,
    title='Tolerance (bases)',
    loc='lower center',
    ncol=len(tolerance_order),
    frameon=True,
  )
  fig.suptitle(title)
  fig.tight_layout(rect=[0, 0.05, 1, 0.98])
  fig.savefig(outfile, dpi=200)
  return long_df

peak_long = tolerance_panels(
  PEAK_TOLERANCE_COLUMNS,
  'peak_tolerance_',
  'Estimated proportion by read count and read-region distance tolerance',
  'proportion_by_peak_tolerance_panels.png',
)
genome_long = tolerance_panels(
  GENOME_TOLERANCE_COLUMNS,
  'genome_tolerance_',
  'Estimated proportion by read count and whole-genome distance tolerance',
  'proportion_by_genome_tolerance_panels.png',
)

print(top_fraction.to_string(index=False))
print(
  df.groupby(['num_reads', 'condition'])['weighted_peak_distance']
  .agg(['median', 'mean', 'size'])
  .to_string()
)
for name, long_df in [('read-region', peak_long), ('whole-genome', genome_long)]:
  print('\nmedian proportion by {} tolerance:'.format(name))
  print(
    long_df.groupby(['condition', 'num_reads', 'tolerance'])['proportion']
    .median()
    .unstack('tolerance')
    .to_string()
  )
print(
  '\nWrote proportion_violins.png, top_haplotype_fraction.png, '
  'weighted_peak_distance_violins.png, proportion_by_peak_tolerance_panels.png, '
  'proportion_by_genome_tolerance_panels.png'
)