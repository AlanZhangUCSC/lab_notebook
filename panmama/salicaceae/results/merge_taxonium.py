#!/usr/bin/env python3
"""
Merge two taxoniumtools JSONL outputs that share topology but differ in
branch lengths and internal-node labels.

Output keeps everything from --primary (including x_dist), but rewrites
internal node names to "{secondary_label}_{primary_label}".
"""
import argparse
import json
import sys
from pathlib import Path


def load_jsonl(path):
  header = None
  nodes = []
  with open(path) as f:
    for i, line in enumerate(f):
      line = line.strip()
      if not line:
        continue
      obj = json.loads(line)
      if i == 0 and "total_nodes" in obj:
        header = obj
      else:
        nodes.append(obj)
  return header, nodes


def build_children_map(nodes):
  children = {n["node_id"]: [] for n in nodes}
  root = None
  for n in nodes:
    pid = n["parent_id"]
    nid = n["node_id"]
    if pid == nid:
      root = nid
    else:
      children[pid].append(nid)
  if root is None:
    raise ValueError("No root found (no node with parent_id == node_id).")
  return children, root


def compute_tip_signatures(nodes, children, root):
  by_id = {n["node_id"]: n for n in nodes}
  sig = {}
  order = []
  stack = [(root, False)]
  while stack:
    nid, processed = stack.pop()
    if processed:
      order.append(nid)
      continue
    stack.append((nid, True))
    for c in children[nid]:
      stack.append((c, False))

  for nid in order:
    node = by_id[nid]
    if node["is_tip"]:
      sig[nid] = (node["name"],)
    else:
      merged = []
      for c in children[nid]:
        merged.extend(sig[c])
      merged.sort()
      sig[nid] = tuple(merged)
  return sig, by_id


def merge(primary_path, secondary_path, out_path):
  hdr_p, nodes_p = load_jsonl(primary_path)
  hdr_s, nodes_s = load_jsonl(secondary_path)

  if len(nodes_p) != len(nodes_s):
    print(f"warning: node counts differ ({len(nodes_p)} vs {len(nodes_s)})",
          file=sys.stderr)

  ch_p, root_p = build_children_map(nodes_p)
  ch_s, root_s = build_children_map(nodes_s)
  sig_p, by_p = compute_tip_signatures(nodes_p, ch_p, root_p)
  sig_s, by_s = compute_tip_signatures(nodes_s, ch_s, root_s)

  s_index = {}
  for nid, s in sig_s.items():
    if not by_s[nid]["is_tip"]:
      s_index[s] = by_s[nid]["name"]

  unmatched = 0
  for n in nodes_p:
    if n["is_tip"]:
      continue
    s = sig_p[n["node_id"]]
    sec_label = s_index.get(s)
    if sec_label is None:
      unmatched += 1
      continue
    n["name"] = f"{sec_label}_{n['name']}"

  if unmatched:
    print(f"warning: {unmatched} internal nodes had no topology match",
          file=sys.stderr)

  with open(out_path, "w") as f:
    if hdr_p is not None:
      f.write(json.dumps(hdr_p) + "\n")
    for n in nodes_p:
      f.write(json.dumps(n) + "\n")
  print(f"wrote {out_path} ({len(nodes_p)} nodes, "
        f"{len(nodes_p) - sum(1 for n in nodes_p if n['is_tip'])} internal, "
        f"{unmatched} unmatched)", file=sys.stderr)


if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument("--primary", required=True,
                  help="JSONL whose x_dist and confidence labels are kept")
  ap.add_argument("--secondary", required=True,
                  help="JSONL whose internal node labels (node_N) are merged in")
  ap.add_argument("-o", "--output", required=True)
  args = ap.parse_args()
  merge(args.primary, args.secondary, args.output)