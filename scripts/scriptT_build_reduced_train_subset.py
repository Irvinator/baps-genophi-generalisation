#!/usr/bin/env python3
import csv
from collections import Counter, defaultdict
from pathlib import Path

base = Path("/mnt/d/baps-genophi-generalisation")
train_path = base / "results" / "baps_full_train_posneg_FINAL.tsv"
cluster_map_path = base / "phage_contig_to_cluster.tsv"
st_map_path = base / "host_to_st_full.tsv"
out_dir = base / "retrain_reduced_dense"
out_dir.mkdir(parents=True, exist_ok=True)

TARGET_HOSTS = 300
TARGET_PHAGES = 300

# host -> ST
host_to_st = {}
with st_map_path.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        host_to_st[row["host_accession"]] = row["ST"]

# phage -> cluster
phage_to_cluster = {}
with cluster_map_path.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        phage_to_cluster[row["phage_contig"]] = row["cluster_id"]

rows = []
host_counts = Counter()
phage_counts = Counter()
st_to_hosts = defaultdict(set)

with train_path.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    fieldnames = r.fieldnames
    for row in r:
        h = row["host_accession"]
        p = row["phage_contig"]
        st = host_to_st.get(h, "UNKNOWN")
        cl = phage_to_cluster.get(p)
        if st == "UNKNOWN" or cl is None:
            continue
        rows.append(row)
        host_counts[h] += 1
        phage_counts[p] += 1
        st_to_hosts[st].add(h)

# choose hosts: one per ST first, then most-connected hosts
selected_hosts = set()
for st, hs in st_to_hosts.items():
    best_h = max(hs, key=lambda h: host_counts[h])
    selected_hosts.add(best_h)

remaining_hosts = sorted(
    [h for h in host_counts if h not in selected_hosts],
    key=lambda h: (-host_counts[h], h)
)
for h in remaining_hosts:
    if len(selected_hosts) >= TARGET_HOSTS:
        break
    selected_hosts.add(h)

rows_h = [row for row in rows if row["host_accession"] in selected_hosts]

# choose most-connected phages among those hosts
phage_counts_h = Counter(row["phage_contig"] for row in rows_h)
selected_phages = set(p for p, _ in phage_counts_h.most_common(TARGET_PHAGES))

subset_rows = [
    row for row in rows_h
    if row["phage_contig"] in selected_phages
]

final_hosts = sorted(set(row["host_accession"] for row in subset_rows))
final_phages = sorted(set(row["phage_contig"] for row in subset_rows))

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

print("target_hosts =", TARGET_HOSTS)
print("target_phages =", TARGET_PHAGES)
print("final_hosts =", len(final_hosts))
print("final_phages =", len(final_phages))
print("subset_rows =", len(subset_rows))
print("wrote =", subset_tsv)
print("wrote =", phenotype_csv)
