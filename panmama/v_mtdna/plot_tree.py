import matplotlib.pyplot as plt
from Bio import Phylo

tree_files = ["v_mtdna.new.nwk"]

fig, axes = plt.subplots(1, len(tree_files), figsize=(15, 7), sharex=True, sharey=True)

if len(tree_files) == 1:
  axes = [axes]

trees = [Phylo.read(tree_file, "newick") for tree_file in tree_files]

max_tree_depth = max(max(tree.depths().values()) for tree in trees)

for tree, ax in zip(trees, axes):
  Phylo.draw(tree, axes=ax, do_show=False)
  ax.set_xlim(0, max_tree_depth)
  ax.set_ylim(-5, 70)
  
  for text in ax.texts:
    text.set_visible(False)
  
  for line in ax.get_lines():
    line.set_linewidth(3)
  
  ax.axis("off")

scale_length = 0.5
scale_y_position = 40
bar_x_start = 0.05 * max_tree_depth
bar_x_end = bar_x_start + scale_length

axes[0].hlines(y=scale_y_position, xmin=bar_x_start, xmax=bar_x_end, color="black", linewidth=3)
axes[0].text((bar_x_start + bar_x_end) / 2, scale_y_position + (0.03 * max_tree_depth), "0.5", 
  ha="center", va="bottom", fontsize=12)

plt.tight_layout()
plt.savefig("panel_trees_scaled_11.11.25.png", dpi=600, bbox_inches='tight')
plt.show()