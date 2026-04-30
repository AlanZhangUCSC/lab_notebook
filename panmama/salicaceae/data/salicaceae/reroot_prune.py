#!/usr/bin/env python3
import sys
import argparse
import dendropy

def main():
  parser = argparse.ArgumentParser(description="Reroot on outgroup, prune it, write Newick to stdout.")
  parser.add_argument("tree", help="Input Newick file")
  parser.add_argument("outgroup", help="Outgroup taxon label")
  args = parser.parse_args()

  tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)

  outgroup_node = tree.find_node_with_taxon_label(args.outgroup)
  if outgroup_node is None:
    sys.exit(f"Error: outgroup '{args.outgroup}' not found in tree.")

  edge_len = outgroup_node.edge.length
  half = (edge_len / 2.0) if edge_len is not None else None
  tree.reroot_at_edge(outgroup_node.edge, length1=half, length2=half, update_bipartitions=False)

  tree.prune_taxa([outgroup_node.taxon])
  tree.purge_taxon_namespace()

  sys.stdout.write(tree.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True))

if __name__ == "__main__":
  main()