#!/usr/bin/env python3
import csv
from collections import Counter
from pathlib import Path

base = Path("/mnt/d/baps-genophi-generalisation")
train_path = base / "results" / "baps_full_train_posneg_FINAL.tsv"
out_dir = base / "retrain_reduced_rows1500"
out_dir.mkdir(parents=True, exist_ok=True)

TARGET_ROWS = 1500

rows = []
host_counts = Counter()
phage_counts = Counter()

with train_path.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    fieldnames = r.fieldnames
    for row in r:
        h = row["host_accession"]
        p = row["phage_contig"]
        rows.append(row)
        host_counts[h] += 1
        phage_counts[p] += 1

# score rows by how connected their host/phage are
def row_score(row):
    h = row["host_accession"]
    p = row["phage_contig"]
    return (host_counts[h] * phage_counts[p], host_counts[h], phage_counts[p])

rows_sorted = sorted(rows, key=row_score, reverse=True)
subset_rows = rows_sorted[:TARGET_ROWS]

final_hosts = sorted(set(r["host_accession"] for r in subset_rows))
final_phages = sorted(set(r["phage_contig"] for r in subset_rows))

subset_tsv = out_dir / "train_subset_posneg.tsv"
with subset_tsv.open("w", newline="") as g:
    w = csv.DictWriter(g, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(subset_rows)

phenotype_csv = out_dir / "train_subset_phenotype_matrix.csv"
with phenotype_csv.open("w", newline="") as g:
    w = csv.DictWriter(g, fieldnames=["strain", "phage", "interaction"])
    w.writeheader()
    for row in subset_rows:
        w.writerow({
            "strain": row["host_accession"],
            "phage": row["phage_contig"],
            "interaction": row["interaction"],
        })

(out_dir / "selected_hosts.txt").write_text("".join(h + "\n" for h in final_hosts))
(out_dir / "selected_phages.txt").write_text("".join(p + "\n" for p in final_phages))

print("target_rows =", TARGET_ROWS)
print("final_rows =", len(subset_rows))
print("final_hosts =", len(final_hosts))
print("final_phages =", len(final_phages))
print("wrote =", subset_tsv)
print("wrote =", phenotype_csv)
