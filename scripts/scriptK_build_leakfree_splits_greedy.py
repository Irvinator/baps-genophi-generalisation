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
        r = dict(zip(header, parts[:3]))
        st = host_to_st.get(r["host_accession"], "UNKNOWN")
        cl = phage_to_cluster.get(r["phage_contig"])
        if st == "UNKNOWN" or cl is None:
            continue
        r["ST"] = st
        r["cluster_id"] = cl
        rows.append(r)

print("usable_rows =", len(rows))

# count rows per ST and cluster
st_counts = Counter(r["ST"] for r in rows)
cl_counts = Counter(r["cluster_id"] for r in rows)

sts = list(st_counts.keys())
cls = list(cl_counts.keys())

random.shuffle(sts)
random.shuffle(cls)

def assign_groups(groups, counts, target=(0.7, 0.15, 0.15)):
    total = sum(counts[g] for g in groups)
    desired = {
        "train": total * target[0],
        "val": total * target[1],
        "test": total * target[2],
    }
    current = {"train": 0, "val": 0, "test": 0}
    assign = {}

    for g in sorted(groups, key=lambda x: counts[x], reverse=True):
        split = min(current, key=lambda s: current[s] / desired[s] if desired[s] > 0 else 1e9)
        assign[g] = split
        current[split] += counts[g]

    return assign, current

st_assign, st_sizes = assign_groups(sts, st_counts)
cl_assign, cl_sizes = assign_groups(cls, cl_counts)

print("ST target sizes =", st_sizes)
print("CL target sizes =", cl_sizes)

# keep only rows where ST split == cluster split
splits = {"train": [], "val": [], "test": []}
dropped = 0

for r in rows:
    s1 = st_assign[r["ST"]]
    s2 = cl_assign[r["cluster_id"]]
    if s1 == s2:
        splits[s1].append(r)
    else:
        dropped += 1

print("dropped_conflicting_rows =", dropped)

for split, data in splits.items():
    out = out_dir / f"baps_full_pos_{split}_greedy.tsv"
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
