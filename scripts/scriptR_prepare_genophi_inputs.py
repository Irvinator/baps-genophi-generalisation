#!/usr/bin/env python3
import csv
import os
from pathlib import Path

base = Path("/mnt/d/baps-genophi-generalisation")
inp = base / "retrain_baseline" / "inputs"
out = base / "retrain_baseline" / "genophi_inputs"

(out / "train_host_AAs").mkdir(parents=True, exist_ok=True)
(out / "train_phage_AAs").mkdir(parents=True, exist_ok=True)
(out / "val_host_AAs").mkdir(parents=True, exist_ok=True)
(out / "val_phage_AAs").mkdir(parents=True, exist_ok=True)

def link_from_tsv(tsv_path, out_dir):
    n = 0
    with open(tsv_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            key = list(row.keys())[0]
            src_key = list(row.keys())[1]
            name = row[key]
            src = Path(row[src_key])
            dst = out_dir / f"{name}.faa"
            if not dst.exists():
                os.symlink(src, dst)
            n += 1
    return n

def write_pheno(in_tsv, out_csv):
    rows = 0
    with open(in_tsv, newline="") as f, open(out_csv, "w", newline="") as g:
        r = csv.DictReader(f, delimiter="\t")
        w = csv.DictWriter(g, fieldnames=["strain", "phage", "interaction"])
        w.writeheader()
        for row in r:
            w.writerow({
                "strain": row["host_accession"],
                "phage": row["phage_contig"],
                "interaction": row["interaction"],
            })
            rows += 1
    return rows

n1 = link_from_tsv(inp / "train_host_files.tsv", out / "train_host_AAs")
n2 = link_from_tsv(inp / "train_phage_files.tsv", out / "train_phage_AAs")
n3 = link_from_tsv(inp / "val_host_files.tsv", out / "val_host_AAs")
n4 = link_from_tsv(inp / "val_phage_files.tsv", out / "val_phage_AAs")

r1 = write_pheno(base / "results" / "baps_full_train_posneg_FINAL.tsv",
                 out / "train_phenotype_matrix.csv")
r2 = write_pheno(base / "results" / "baps_full_val_posneg_FINAL.tsv",
                 out / "val_phenotype_matrix.csv")

print("train_host_links =", n1)
print("train_phage_links =", n2)
print("val_host_links =", n3)
print("val_phage_links =", n4)
print("train_rows =", r1)
print("val_rows =", r2)
