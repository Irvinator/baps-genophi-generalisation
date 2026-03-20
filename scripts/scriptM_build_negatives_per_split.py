#!/usr/bin/env python3
import csv
import random
from pathlib import Path

random.seed(42)

splits = {
    "train": Path("/mnt/d/baps-genophi-generalisation/results/baps_full_pos_train_FINAL.tsv"),
    "val":   Path("/mnt/d/baps-genophi-generalisation/results/baps_full_pos_val_FINAL.tsv"),
    "test":  Path("/mnt/d/baps-genophi-generalisation/results/baps_full_pos_test_FINAL.tsv"),
}

out_dir = Path("/mnt/d/baps-genophi-generalisation/results")
neg_per_pos = 1  # change later if needed

for split, path in splits.items():
    rows = []
    hosts = set()
    phages = set()
    pos_pairs = set()

    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)
            h = row["host_accession"]
            p = row["phage_contig"]
            hosts.add(h)
            phages.add(p)
            pos_pairs.add((h, p))

    hosts = sorted(hosts)
    phages = sorted(phages)

    pos_by_host = {}
    for h in hosts:
        pos_by_host[h] = set()

    for h, p in pos_pairs:
        pos_by_host[h].add(p)

    neg_rows = []
    for h in hosts:
        pos_set = pos_by_host[h]
        n_pos = len(pos_set)
        n_neg = neg_per_pos * n_pos
        candidates = [p for p in phages if p not in pos_set]
        if not candidates or n_neg == 0:
            continue
        chosen = random.sample(candidates, k=min(n_neg, len(candidates)))
        for p in chosen:
            neg_rows.append({
                "host_accession": h,
                "phage_contig": p,
                "interaction": "0"
            })

    pos_out = []
    for row in rows:
        pos_out.append({
            "host_accession": row["host_accession"],
            "phage_contig": row["phage_contig"],
            "interaction": "1"
        })

    final = pos_out + neg_rows
    random.shuffle(final)

    out_path = out_dir / f"baps_full_{split}_posneg_FINAL.tsv"
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["host_accession", "phage_contig", "interaction"], delimiter="\t")
        w.writeheader()
        w.writerows(final)

    print(split, "pos=", len(pos_out), "neg=", len(neg_rows), "total=", len(final), "->", out_path)
