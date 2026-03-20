#!/usr/bin/env python3
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

pairs_path = Path("/home/irvin/data/BAPS/outputs/baps_ecoli_pos_pairs_in_fasta.tsv")
st_path = Path("/mnt/d/baps-genophi-generalisation/host_to_st_full.tsv")
cluster_path = Path("/mnt/d/baps-genophi-generalisation/phage_contig_to_cluster.tsv")
out_dir = Path("/mnt/d/baps-genophi-generalisation/results")
out_dir.mkdir(parents=True, exist_ok=True)

random.seed(42)

host_to_st = {}
with st_path.open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        host_to_st[r["host_accession"]] = r["ST"]

phage_to_cluster = {}
with cluster_path.open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        phage_to_cluster[r["phage_contig"]] = r["cluster_id"]

rows = []
with pairs_path.open() as f:
    header = f.readline().strip().split()
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        r = dict(zip(header, parts[:3]))
        st = host_to_st.get(r["host_accession"], "UNKNOWN")
        cl = phage_to_cluster.get(r["phage_contig"])
        if st == "UNKNOWN" or cl is None:
            continue
        r["ST"] = st
        r["cluster_id"] = cl
        rows.append(r)

print("usable_rows =", len(rows))

st_counts = Counter(r["ST"] for r in rows)
sts = list(st_counts)
random.shuffle(sts)

target = {"train": 0.7, "val": 0.15, "test": 0.15}
total = sum(st_counts.values())
desired = {k: total * v for k, v in target.items()}
current = {"train": 0, "val": 0, "test": 0}
st_split = {}

for st in sorted(sts, key=lambda x: st_counts[x], reverse=True):
    split = min(current, key=lambda s: current[s] / desired[s] if desired[s] > 0 else 1e9)
    st_split[st] = split
    current[split] += st_counts[st]

print("ST_assigned_rows =", current)

cluster_votes = defaultdict(Counter)
for r in rows:
    cluster_votes[r["cluster_id"]][st_split[r["ST"]]] += 1

cluster_split = {}
for cl, cnt in cluster_votes.items():
    best = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    cluster_split[cl] = best

splits = {"train": [], "val": [], "test": []}
dropped = 0

for r in rows:
    a = st_split[r["ST"]]
    b = cluster_split[r["cluster_id"]]
    if a == b:
        splits[a].append(r)
    else:
        dropped += 1

print("dropped_conflicting_rows =", dropped)

for split, data in splits.items():
    out = out_dir / f"baps_full_pos_{split}_majority.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["host_accession", "phage_contig", "interaction", "ST", "cluster_id"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(data)
    print(split, "rows =", len(data), "->", out)

split_st = {k: {r["ST"] for r in v} for k, v in splits.items()}
split_cl = {k: {r["cluster_id"] for r in v} for k, v in splits.items()}

for a, b in [("train", "val"), ("train", "test"), ("val", "test")]:
    print(f"ST_overlap_{a}_{b} =", len(split_st[a] & split_st[b]))
    print(f"CL_overlap_{a}_{b} =", len(split_cl[a] & split_cl[b]))
