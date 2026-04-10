#!/usr/bin/env python3
import os
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import numpy as np

BASIC_CSV = "/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_basic_features.csv"
PROTEIN_DIR = Path("/mnt/d/baps-genophi-generalisation/aim2_ecoli96/prodigal_proteins")
OUTPUT = "/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_feature_table.csv"

basic = pd.read_csv(BASIC_CSV)

rows = []

for _, row in basic.iterrows():
    phage_id = row["phage_id"]
    genome_length = row["genome_length"]

    faa = PROTEIN_DIR / f"{phage_id}_proteins.faa"

    if not faa.exists():
        rows.append({
            "phage_id": phage_id,
            "num_predicted_proteins": np.nan,
            "mean_protein_length_aa": np.nan,
            "median_protein_length_aa": np.nan,
            "total_protein_aa": np.nan,
            "coding_density_proxy": np.nan,
            "proteins_per_10kb": np.nan
        })
        continue

    lengths = [len(rec.seq) for rec in SeqIO.parse(str(faa), "fasta")]

    if len(lengths) == 0:
        rows.append({
            "phage_id": phage_id,
            "num_predicted_proteins": 0,
            "mean_protein_length_aa": 0,
            "median_protein_length_aa": 0,
            "total_protein_aa": 0,
            "coding_density_proxy": 0,
            "proteins_per_10kb": 0
        })
        continue

    total_protein_aa = sum(lengths)

    rows.append({
        "phage_id": phage_id,
        "num_predicted_proteins": len(lengths),
        "mean_protein_length_aa": np.mean(lengths),
        "median_protein_length_aa": np.median(lengths),
        "total_protein_aa": total_protein_aa,
        # convert aa to nt proxy by multiplying by 3
        "coding_density_proxy": (total_protein_aa * 3) / genome_length if genome_length > 0 else np.nan,
        "proteins_per_10kb": len(lengths) / (genome_length / 10000) if genome_length > 0 else np.nan
    })

protein_df = pd.DataFrame(rows)

merged = basic.merge(protein_df, on="phage_id", how="left")
merged.to_csv(OUTPUT, index=False)

print("Saved to:", OUTPUT)
print()
print(merged.head())
print()
print("Shape:", merged.shape)
