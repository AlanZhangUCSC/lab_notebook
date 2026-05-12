import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict

parser = argparse.ArgumentParser(description="Calculate family distribution from mgsr input.")
parser.add_argument("mgsr_in", help="Input mgsr file")
parser.add_argument("output_prefix", help="Prefix for output files")
args = parser.parse_args()
mgsr_in = args.mgsr_in
output_prefix = args.output_prefix

family_counts = defaultdict(int)
family_counts["."] = 0
with open(mgsr_in, 'r') as fh:
  for line in fh:
    nodes_string, family, count, read_indices = line.strip().split('\t')
    family_counts[family] += int(count)

family_count_list = list(family_counts.items())
family_count_list.sort(key=lambda x: x[1], reverse=True)

total_counts = sum(count for _, count in family_count_list)

with open(output_prefix + '.family_distribution.tsv', 'w') as fh:
  fh.write('family\tcount\tproporiton(%%)\n')
  for family, count in family_count_list:
    if family == ".": continue
    fh.write(f"{family}\t{count}\t{count/total_counts*100:.5f}\n")
  fh.write(f".\t{family_counts['.']}\t{family_counts['.']/total_counts*100:.5f}\n")

plot_data = sorted(
  family_count_list,
  key=lambda x: (x[0] == ".", x[0].lower())
)

sns.set_theme(style="whitegrid")
labels = [f if f != "." else "Unknown" for f, _ in plot_data]
counts = np.array([c for _, c in plot_data], dtype=float)
total = counts.sum()
percentages = counts / total * 100

n_known = sum(1 for f, _ in plot_data if f != ".")
palette = sns.color_palette("husl", n_colors=max(n_known, 1))
colors = []
known_idx = 0
for f, _ in plot_data:
  if f == ".":
    colors.append((0.6, 0.6, 0.6))
  else:
    colors.append(palette[known_idx])
    known_idx += 1

INSIDE_THRESHOLD = 3.0

fig, ax = plt.subplots(figsize=(12, 9))
wedges, _ = ax.pie(
  counts,
  colors=colors,
  startangle=90,
  wedgeprops={'edgecolor': 'white', 'linewidth': 1}
)

for i, wedge in enumerate(wedges):
  angle = (wedge.theta2 + wedge.theta1) / 2.0
  pct = percentages[i]
  label_text = f"{labels[i]} ({pct:.1f}%)"

  if pct >= INSIDE_THRESHOLD:
    x = 0.65 * np.cos(np.deg2rad(angle))
    y = 0.65 * np.sin(np.deg2rad(angle))
    ax.text(
      x, y, f"{labels[i]}\n{pct:.1f}%",
      ha='center', va='center',
      fontsize=10, fontweight='bold', color='white'
    )
  else:
    x_outer = np.cos(np.deg2rad(angle))
    y_outer = np.sin(np.deg2rad(angle))
    x_text = 1.35 * x_outer
    y_text = 1.35 * y_outer
    ha = 'left' if x_outer >= 0 else 'right'
    connectionstyle = f"angle,angleA=0,angleB={angle}"
    ax.annotate(
      label_text,
      xy=(x_outer, y_outer),
      xytext=(x_text, y_text),
      ha=ha, va='center',
      fontsize=9,
      arrowprops=dict(
        arrowstyle='-',
        connectionstyle=connectionstyle,
        color='gray',
        lw=0.7
      )
    )

ax.set_title('Alu Family Distribution', fontsize=14, fontweight='bold', pad=20)
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.6, 1.6)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(output_prefix + '.family_distribution.png', dpi=300, bbox_inches='tight')
plt.close()