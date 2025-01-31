from ete3 import Tree
import sys

# Read the tree from file
tree = Tree(sys.argv[1], format=1)  # format=1 is for standard Newick

# Print parent-child relationships
for node in tree.traverse("preorder"):  # preorder ensures depth-first traversal
  if not node.is_root():
    print(node.up.name, node.name)
