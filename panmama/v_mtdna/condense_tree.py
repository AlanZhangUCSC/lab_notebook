#!/usr/bin/env python3
"""
Condense a phylogenetic tree by taxonomic family and subsample using greedy
phylogenetic diversity (PD) maximization.

Algorithm:
1. Parse Newick tree and family annotations
2. Build condensed tree with one tip per family using UPGMA on pairwise distances
3. Greedily select families that maximize total branch length coverage

Time complexity: O(n) parsing + O(f²) UPGMA + O(f²k) greedy selection
  where n=samples, f=families, k=target count
Memory: O(n + f²) for tree structure and distance matrix
"""

import sys
import csv
from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class Node:
    name: str = ""
    branch_length: float = 0.0
    children: list = field(default_factory=list)
    parent: Optional['Node'] = None
    family: str = ""
    
    @property
    def is_tip(self) -> bool:
        return len(self.children) == 0


class NewickParser:
    def __init__(self, newick: str):
        self.s = newick.strip()
        self.pos = 0
    
    def parse(self) -> Node:
        root = self._parse_node()
        if self.pos < len(self.s) and self.s[self.pos] == ';':
            self.pos += 1
        return root
    
    def _parse_node(self) -> Node:
        node = Node()
        
        if self.pos < len(self.s) and self.s[self.pos] == '(':
            self.pos += 1
            while True:
                child = self._parse_node()
                child.parent = node
                node.children.append(child)
                if self.pos < len(self.s) and self.s[self.pos] == ',':
                    self.pos += 1
                else:
                    break
            if self.pos < len(self.s) and self.s[self.pos] == ')':
                self.pos += 1
        
        label = []
        while self.pos < len(self.s) and self.s[self.pos] not in ':,);(':
            label.append(self.s[self.pos])
            self.pos += 1
        node.name = ''.join(label)
        
        if self.pos < len(self.s) and self.s[self.pos] == ':':
            self.pos += 1
            len_chars = []
            while self.pos < len(self.s) and self.s[self.pos] in '0123456789.eE-+':
                len_chars.append(self.s[self.pos])
                self.pos += 1
            if len_chars:
                node.branch_length = float(''.join(len_chars))
        
        return node


def collect_tips(node: Node) -> list[Node]:
    if node.is_tip:
        return [node]
    tips = []
    for child in node.children:
        tips.extend(collect_tips(child))
    return tips


def to_newick(node: Node) -> str:
    if node.is_tip:
        return f"{node.name}:{node.branch_length}"
    children_str = ','.join(to_newick(c) for c in node.children)
    return f"({children_str}){node.name}:{node.branch_length}"


def load_annotations(path: str) -> dict[str, str]:
    sample_to_family = {}
    with open(path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 2:
                sample_to_family[row[0]] = row[1]
    return sample_to_family


def assign_families(node: Node, annotations: dict[str, str]):
    if node.is_tip:
        node.family = annotations.get(node.name, f"Unknown_{node.name}")
    for child in node.children:
        assign_families(child, annotations)


def compute_pairwise_distances(tips: list[Node], families: list[str]) -> np.ndarray:
    """
    Compute pairwise distances between family representatives.
    Uses path-to-root intersection to find MRCA distance efficiently.
    """
    n = len(families)
    family_to_idx = {f: i for i, f in enumerate(families)}
    
    family_tips = {}
    for tip in tips:
        if tip.family not in family_tips:
            family_tips[tip.family] = tip
    
    reps = [family_tips[f] for f in families]
    
    paths = []
    for rep in reps:
        path = {}
        dist = 0.0
        curr = rep
        while curr is not None:
            path[id(curr)] = dist
            dist += curr.branch_length
            curr = curr.parent
        paths.append(path)
    
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            min_dist = float('inf')
            for node_id, d_i in paths[i].items():
                if node_id in paths[j]:
                    min_dist = min(min_dist, d_i + paths[j][node_id])
            dist_matrix[i, j] = dist_matrix[j, i] = min_dist
    
    return dist_matrix


def build_upgma_tree(families: list[str], dist_matrix: np.ndarray) -> Node:
    """
    Build UPGMA tree from distance matrix.
    UPGMA produces ultrametric trees suitable for visualization.
    """
    n = len(families)
    
    nodes = [Node(name=f, family=f) for f in families]
    heights = np.zeros(n)
    cluster_size = np.ones(n, dtype=int)
    cluster_id = np.arange(n)
    
    dist = dist_matrix.copy()
    np.fill_diagonal(dist, np.inf)
    
    for _ in range(n - 1):
        active = np.where(cluster_id == np.arange(n))[0]
        if len(active) < 2:
            break
        
        min_val = np.inf
        mi, mj = -1, -1
        for i in active:
            for j in active:
                if i < j and dist[i, j] < min_val:
                    min_val = dist[i, j]
                    mi, mj = i, j
        
        if mi < 0:
            break
        
        new_height = min_val / 2.0
        
        merged = Node()
        nodes[mi].branch_length = new_height - heights[mi]
        nodes[mj].branch_length = new_height - heights[mj]
        nodes[mi].parent = merged
        nodes[mj].parent = merged
        merged.children = [nodes[mi], nodes[mj]]
        
        nodes[mi] = merged
        heights[mi] = new_height
        
        si, sj = cluster_size[mi], cluster_size[mj]
        cluster_size[mi] = si + sj
        
        for k in range(n):
            if cluster_id[k] == k and k != mi:
                dist[mi, k] = dist[k, mi] = (dist[mi, k] * si + dist[mj, k] * sj) / (si + sj)
        
        cluster_id[mj] = mi
        dist[mj, :] = np.inf
        dist[:, mj] = np.inf
    
    root_idx = np.where(cluster_id == np.arange(n))[0][0]
    return nodes[root_idx]


def greedy_pd_subsample(root: Node, required_family: str, target_count: int) -> list[str]:
    """
    Greedy phylogenetic diversity maximization.
    
    At each step, select the family that adds the most unique branch length
    to the current selection. This achieves a (1 - 1/e) approximation of
    optimal PD, which is the best polynomial-time guarantee for this NP-hard problem.
    """
    tips = collect_tips(root)
    n = len(tips)
    name_to_idx = {tip.name: i for i, tip in enumerate(tips)}
    
    paths = []
    for tip in tips:
        path = []
        curr = tip
        dist = 0.0
        while curr is not None:
            path.append((id(curr), curr.branch_length, dist))
            dist += curr.branch_length
            curr = curr.parent
        paths.append(path)
    
    selected = []
    covered_edges = set()
    is_selected = [False] * n
    
    if required_family in name_to_idx:
        idx = name_to_idx[required_family]
        is_selected[idx] = True
        selected.append(required_family)
        for node_id, _, _ in paths[idx]:
            covered_edges.add(node_id)
    
    while len(selected) < target_count and len(selected) < n:
        best_gain = -1.0
        best_idx = -1
        
        for i in range(n):
            if is_selected[i]:
                continue
            
            gain = sum(bl for node_id, bl, _ in paths[i] if node_id not in covered_edges)
            
            if gain > best_gain:
                best_gain = gain
                best_idx = i
        
        if best_idx < 0:
            break
        
        is_selected[best_idx] = True
        selected.append(tips[best_idx].name)
        for node_id, _, _ in paths[best_idx]:
            covered_edges.add(node_id)
    
    return selected


def prune_tree(node: Node, keep_families: set[str]) -> Optional[Node]:
    """Prune tree to keep only specified families, collapsing single-child nodes."""
    if node.is_tip:
        if node.name in keep_families:
            return Node(name=node.name, branch_length=node.branch_length, family=node.family)
        return None
    
    kept_children = []
    for child in node.children:
        pruned = prune_tree(child, keep_families)
        if pruned is not None:
            kept_children.append(pruned)
    
    if not kept_children:
        return None
    
    if len(kept_children) == 1:
        kept_children[0].branch_length += node.branch_length
        return kept_children[0]
    
    new_node = Node(name=node.name, branch_length=node.branch_length)
    for child in kept_children:
        child.parent = new_node
    new_node.children = kept_children
    return new_node


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <tree.nwk> <annotations.tsv> <output.nwk> <num_families> [required_family]")
        print("\nArguments:")
        print("  tree.nwk        Input Newick tree file")
        print("  annotations.tsv TSV with columns: sample_name, family")
        print("  output.nwk      Output Newick file")
        print("  num_families    Number of families to keep (e.g., 100)")
        print("  required_family Family that must be included (default: Elephantidae)")
        sys.exit(1)
    
    tree_path = sys.argv[1]
    annot_path = sys.argv[2]
    output_path = sys.argv[3]
    target_families = int(sys.argv[4])
    required_family = sys.argv[5] if len(sys.argv) > 5 else "Elephantidae"
    
    print(f"Loading tree from {tree_path}...", file=sys.stderr)
    with open(tree_path, 'r') as f:
        newick = f.read()
    
    parser = NewickParser(newick)
    root = parser.parse()
    
    tips = collect_tips(root)
    print(f"Loaded tree with {len(tips)} tips", file=sys.stderr)
    
    print(f"Loading annotations from {annot_path}...", file=sys.stderr)
    annotations = load_annotations(annot_path)
    print(f"Loaded {len(annotations)} sample annotations", file=sys.stderr)
    
    assign_families(root, annotations)
    
    family_tips = {}
    for tip in tips:
        if tip.family not in family_tips:
            family_tips[tip.family] = []
        family_tips[tip.family].append(tip)
    
    families = list(family_tips.keys())
    print(f"Found {len(families)} unique families", file=sys.stderr)
    
    print("Computing pairwise distances...", file=sys.stderr)
    dist_matrix = compute_pairwise_distances(tips, families)
    
    print("Building condensed UPGMA tree...", file=sys.stderr)
    condensed_root = build_upgma_tree(families, dist_matrix)
    
    family_nodes = collect_tips(condensed_root)
    print(f"Condensed tree has {len(family_nodes)} family tips", file=sys.stderr)
    
    has_required = any(n.name == required_family for n in family_nodes)
    if not has_required:
        print(f"Warning: Required family '{required_family}' not found in tree", file=sys.stderr)
    
    print(f"Subsampling to {target_families} families using greedy PD...", file=sys.stderr)
    selected = greedy_pd_subsample(condensed_root, required_family, target_families)
    print(f"Selected {len(selected)} families", file=sys.stderr)
    
    keep_set = set(selected)
    final_tree = prune_tree(condensed_root, keep_set)
    
    output_newick = to_newick(final_tree) + ";"
    with open(output_path, 'w') as f:
        f.write(output_newick)
    print(f"Wrote subsampled tree to {output_path}", file=sys.stderr)
    
    families_path = output_path + ".families.txt"
    with open(families_path, 'w') as f:
        for fam in selected:
            f.write(fam + "\n")
    print(f"Wrote selected families list to {families_path}", file=sys.stderr)


if __name__ == "__main__":
    main()