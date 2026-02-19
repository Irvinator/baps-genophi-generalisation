import pandas as pd
import re
from pathlib import Path

IN_TSV = Path.home() / "data/BAPS/baps_ecoli.tsv"
OUT_DIR = Path.home() / "data/BAPS/outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# baps_ecoli.tsv is tab-separated and begins with a leading tab in the header
df = pd.read_csv(IN_TSV, sep="\t")

# Columns expected from your file: org, sample, contig, n_genes, n_anns, genes
# sample looks like: "GCA_018015695.1_-_PDT000998088.1" etc.
def extract_gca(sample: str):
    m = re.search(r"(GCA_\d+\.\d+)", str(sample))
    return m.group(1) if m else None

df["host_accession"] = df["sample"].apply(extract_gca)

pos = df[["host_accession", "contig"]].dropna().drop_duplicates()
pos = pos.rename(columns={"contig": "phage_contig"})
pos["interaction"] = 1

pos_path = OUT_DIR / "baps_ecoli_pos_pairs.tsv"
pos.to_csv(pos_path, sep="\t", index=False)

acc_path = OUT_DIR / "baps_ecoli_accessions.txt"
pos["host_accession"].drop_duplicates().to_csv(acc_path, index=False, header=False)

print(f"Wrote positives: {pos_path} (rows={len(pos)})")
print(f"Wrote accessions: {acc_path} (unique={pos['host_accession'].nunique()})")
print("Example rows:")
print(pos.head(5).to_string(index=False))
