import urllib.request, tarfile, io, csv

url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
with urllib.request.urlopen(url) as r:
  buf = io.BytesIO(r.read())

parent = {}
with tarfile.open(fileobj=buf, mode="r:gz") as tar:
  with tar.extractfile("nodes.dmp") as f:
    for line in f:
      parts = line.decode().split("\t|\t")
      parent[int(parts[0])] = int(parts[1])

def is_primate(tid, root=9443):
  seen = set()
  while tid in parent and tid not in seen:
    if tid == root: return True
    seen.add(tid)
    nxt = parent[tid]
    if nxt == tid: return False
    tid = nxt
  return False

with open("all_ucsc_assemblies.tsv") as fin, \
     open("ucsc_primate_assemblies.tsv", "w") as fout:
  for row in csv.reader(fin, delimiter="\t"):
    if len(row) >= 4 and row[3].isdigit() and is_primate(int(row[3])):
      fout.write("\t".join(row) + "\n")