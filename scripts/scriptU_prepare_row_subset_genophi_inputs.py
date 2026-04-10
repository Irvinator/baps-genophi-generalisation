#!/usr/bin/env python3
import os
from pathlib import Path

base = Path("/mnt/d/baps-genophi-generalisation")
subset_dir = base / "retrain_reduced_rows1500"
out_dir = subset_dir / "genophi_inputs"
(out_dir / "train_host_AAs").mkdir(parents=True, exist_ok=True)
(out_dir / "train_phage_AAs").mkdir(parents=True, exist_ok=True)

host_src = base / "host_AAs_baps"
phage_src = Path("/mnt/d/BAPS/phage_AAs_baps")
hosts = [x.strip() for x in (subset_dir / "selected_hosts.txt").read_text().splitlines() if x.strip()]
phages = [x.strip() for x in (subset_dir / "selected_phages.txt").read_text().splitlines() if x.strip()]

host_n = 0
for h in hosts:
    src = host_src / f"{h}.faa"
    dst = out_dir / "train_host_AAs" / f"{h}.faa"
    if src.exists():
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.symlink(src, dst)
        host_n += 1

phage_n = 0
for p in phages:
    src = phage_src / f"{p}.faa"
    dst = out_dir / "train_phage_AAs" / f"{p}.faa"
    if src.exists():
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.symlink(src, dst)
        phage_n += 1

print("host_links =", host_n)
print("phage_links =", phage_n)
