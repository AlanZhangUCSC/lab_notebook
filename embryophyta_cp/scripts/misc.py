def sister_leaf_names(leaf):
  parent = leaf.up
  if parent is None:
    return []
  return [l.name for l in parent.get_leaves() if l is not leaf]

def closest_relative_names(leaf):
  sisters = sister_leaf_names(leaf)
  if len(sisters) >= 2:
    return sisters
  parent = leaf.up
  if parent is None or parent.up is None:
    return sisters
  extended = list(sisters)
  for sib in parent.get_sisters():
    extended.extend(l.name for l in sib.get_leaves())
  return extended