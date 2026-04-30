import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('out/block_mutations_optimized.tsv', sep='\t')

mutation_types = ['block_translocations', 'block_insertions', 'block_deletions', 
                  'nuc_substitutions', 'nuc_insertions', 'nuc_deletions']

reduction_data = []
for _, row in df.iterrows():
  tree = row['tree']
  for mut_type in mutation_types:
    original = row[mut_type]
    optimized = row[f'{mut_type}_optimized']
    reduction_pct = ((original - optimized) / original) * 100 if original > 0 else 0
    reduction_data.append({
      'tree': tree.upper(),
      'mutation_type': mut_type.replace('_', ' ').title(),
      'reduction_percent': reduction_pct
    })

reduction_df = pd.DataFrame(reduction_data)
trees = ['RSV', 'SARS', 'HIV']

color_palettes = {
  'palette1': ['#E63946', '#457B9D', '#2A9D8F'],
  'palette2': ['#8338EC', '#FB5607', '#3A86FF'],
  'palette3': ['#D62828', '#F77F00', '#06A77D'],
  'palette4': ['#780000', '#C1121F', '#FFC300'],
  'palette5': ['#4A4E69', '#9A8C98', '#C9ADA7']
}

for palette_name, colors in color_palettes.items():
  fig, ax = plt.subplots(figsize=(14, 8))
  
  pivot_reduction = reduction_df.pivot(index='mutation_type', columns='tree', values='reduction_percent')
  pivot_reduction = pivot_reduction.reindex(['Block Translocations', 'Block Insertions', 'Block Deletions',
                                             'Nuc Substitutions', 'Nuc Insertions', 'Nuc Deletions'])
  
  x = np.arange(len(pivot_reduction.index))
  width = 0.25
  
  for i, tree in enumerate(trees):
    bars = ax.bar(x + i*width, pivot_reduction[tree], width, label=tree, color=colors[i], edgecolor='black', linewidth=0.5)
    for bar in bars:
      height = bar.get_height()
      ax.text(bar.get_x() + bar.get_width()/2., height,
              f'{height:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
  
  ax.set_ylabel('Mutation Reduction (%)', fontsize=13, fontweight='bold')
  ax.set_xlabel('Mutation Category', fontsize=13, fontweight='bold')
  ax.set_title('Tree Optimization: Mutation Reduction Across Categories', fontsize=15, fontweight='bold', pad=20)
  ax.set_xticks(x + width)
  ax.set_xticklabels(pivot_reduction.index, rotation=45, ha='right', fontsize=11)
  ax.legend(title='Tree', fontsize=11, title_fontsize=12)
  ax.grid(axis='y', alpha=0.3, linestyle='--')
  ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
  ax.set_ylim(bottom=min(0, pivot_reduction.min().min() - 5))
  
  plt.tight_layout()
  plt.savefig(f'out/mutation_reduction_{palette_name}.png', dpi=300, bbox_inches='tight')
  plt.close()
  print(f"Saved {palette_name}")

print("\nColor Palette Details:")
print("palette1: Red, Blue, Teal (high contrast)")
print("palette2: Purple, Orange, Bright Blue (vibrant)")
print("palette3: Dark Red, Orange, Green (warm to cool)")
print("palette4: Maroon, Crimson, Gold (warm gradient)")
print("palette5: Gray-Blue, Mauve, Tan (muted/professional)")