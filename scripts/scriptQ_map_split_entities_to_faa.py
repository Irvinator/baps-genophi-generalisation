#!/usr/bin/env python3
import csv
from pathlib import Path

base = Path("/mnt/d/baps-genophi-generalisation/retrain_baseline/inputs")

host_dir = Path("/mnt/d/baps-genophi-generalisation/host_AAs_baps")
phage_dir = Path("/home/irvin/data/BAPS/phage_AAs_baps")

host_map = {p.stem: str(p) for p in host_dir.glob("*.faa")}
phage_map = {p.stem: str(p) for p in phage_dir.glob("*.faa")}

for split in ["train", "val"]:
    hosts = [x.strip() for x in (base / f"{split}_hosts.txt").read_text().splitlines() if x.strip()]
    phages = [x.strip() for x in (base / f"{split}_phages.txt").read_text().splitlines() if x.strip()]

    host_out = base / f"{split}_host_files.tsv"
    phage_out = base / f"{split}_phage_files.tsv"

    n_host = 0
    with host_out.open("w", newline="") as g:
        g.write("strain\tfaa_path\n")
        for h in hosts:
            if h in host_map:
                g.write(f"{h}\t{host_map[h]}\n")
                n_host += 1

    n_phage = 0
    with phage_out.open("w", newline="") as g:
        g.write("phage\tfaa_path\n")
        for p in phages:
            if p in phage_map:
                g.write(f"{p}\t{phage_map[p]}\n")
                n_phage += 1

    print(split, "host_files_found=", n_host, "of", len(hosts), "phage_files_found=", n_phage, "of", len(phages))
