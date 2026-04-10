#!/usr/bin/env python3
import csv
from collections import defaultdict, deque
from pathlib import Path

pairs_path = Path("/home/irvin/data/BAPS/outputs/baps_ecoli_pos_pairs_in_fasta.tsv")

st_path = Path("/mnt/d/baps-genophi-generalisation/host_to_st_full.tsv")

cluster_path = Path("/mnt/d/baps-genophi-generalisation/phage_contig_to_cluster.tsv")

out_dir = Path("/mnt/d/baps-genophi-generalisation/results")
out_dir.mkdir(parents=True, exist_ok=True)

# host -> ST
host_to_st = {}
with st_path.open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        host_to_st[r["host_accession"]] = r["ST"]

# phage -> cluster
phage_to_cluster = {}
with cluster_path.open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        phage_to_cluster[r["phage_contig"]] = r["cluster_id"]

# read pairs
rows = []
with pairs_path.open() as f:
    header = f.readline().strip().split()
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        rows.append(dict(zip(header, parts[:3])))

annot = []
missing_st = 0
missing_cl = 0

for r in rows:
    h = r["host_accession"]
    p = r["phage_contig"]
    st = host_to_st.get(h, "UNKNOWN")
    cl = phage_to_cluster.get(p)

    if st == "UNKNOWN":
        missing_st += 1
    if cl is None:
        missing_cl += 1
        continue

    rr = dict(r)
    rr["ST"] = st
    rr["cluster_id"] = cl

    annot.append(rr)
annot = [r for r in annot if r["ST"] != "UNKNOWN"]
print("annot_rows_after_drop_UNKNOWN_ST =", len(annot))
print("input_rows =", len(rows))
print("annot_rows =", len(annot))
print("missing_ST_rows =", missing_st)
print("missing_cluster_rows_dropped =", missing_cl)

# bipartite graph ST <-> cluster
adj = defaultdict(set)
edge_rows = defaultdict(list)

for r in annot:
    s = "ST::" + r["ST"]
    c = "CL::" + r["cluster_id"]
    adj[s].add(c)
    adj[c].add(s)
    edge_rows[(s, c)].append(r)

# connected components
seen = set()
components = []

for node in adj:
    if node in seen:
        continue
    q = deque([node])
    seen.add(node)
    comp_rows = []

    while q:
        x = q.popleft()
        for y in adj[x]:
            edge = (x, y) if x.startswith("ST::") else (y, x)
            comp_rows.extend(edge_rows.get(edge, []))
            if y not in seen:
                seen.add(y)
                q.append(y)

    uniq = []
    used = set()
    for r in comp_rows:
        key = (r["host_accession"], r["phage_contig"], r["interaction"], r["ST"], r["cluster_id"])
        if key not in used:
            used.add(key)
            uniq.append(r)

    components.append(uniq)

components.sort(key=len, reverse=True)

print("n_components =", len(components))
print("largest_component_rows =", len(components[0]) if components else 0)

# greedy split
target = {"train": 0.70, "val": 0.15, "test": 0.15}
total = sum(len(c) for c in components)
desired = {k: total * v for k, v in target.items()}

splits = {"train": [], "val": [], "test": []}
sizes = {"train": 0, "val": 0, "test": 0}

for comp in components:
    split = min(
        sizes,
        key=lambda k: sizes[k] / desired[k] if desired[k] > 0 else 1e9
    )
    splits[split].extend(comp)
    sizes[split] += len(comp)

# write outputs
for split, data in splits.items():
    out = out_dir / f"baps_full_pos_{split}.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["host_accession", "phage_contig", "interaction", "ST", "cluster_id"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(data)
    print(split, "rows =", len(data), "->", out)

# leakage checks
split_st = {k: {r["ST"] for r in v} for k, v in splits.items()}
split_cl = {k: {r["cluster_id"] for r in v} for k, v in splits.items()}

for a, b in [("train", "val"), ("train", "test"), ("val", "test")]:
    print(f"ST_overlap_{a}_{b} =", len(split_st[a] & split_st[b]))
    print(f"CL_overlap_{a}_{b} =", len(split_cl[a] & split_cl[b]))
