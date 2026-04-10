#!/usr/bin/env python3
import csv
from pathlib import Path

files = {
    "train": Path("/mnt/d/baps-genophi-generalisation/results/baps_full_train_posneg_FINAL.tsv"),
    "val":   Path("/mnt/d/baps-genophi-generalisation/results/baps_full_val_posneg_FINAL.tsv"),
}

out_dir = Path("/mnt/d/baps-genophi-generalisation/retrain_baseline/inputs")
out_dir.mkdir(parents=True, exist_ok=True)

for split, path in files.items():
    hosts = []
    phages = []
    seen_h = set()
    seen_p = set()

    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            h = row["host_accession"]
            p = row["phage_contig"]
            if h not in seen_h:
                seen_h.add(h)
                hosts.append(h)
            if p not in seen_p:
                seen_p.add(p)
                phages.append(p)

    with (out_dir / f"{split}_hosts.txt").open("w") as g:
        for h in hosts:
            g.write(h + "\n")

    with (out_dir / f"{split}_phages.txt").open("w") as g:
        for p in phages:
            g.write(p + "\n")

    print(split, "hosts=", len(hosts), "phages=", len(phages))
