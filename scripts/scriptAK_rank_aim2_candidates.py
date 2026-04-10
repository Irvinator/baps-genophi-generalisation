#!/usr/bin/env python3
import pandas as pd
import numpy as np

INPUT = "/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_feature_table.csv"
OUTPUT = "/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_ranked_candidates.csv"

df = pd.read_csv(INPUT)

# z-score normalisation
for col in ["coding_density_proxy", "num_predicted_proteins", "proteins_per_10kb"]:
    df[f"{col}_z"] = (df[col] - df[col].mean()) / df[col].std()

# simple composite score
df["candidate_priority_score"] = (
    df["coding_density_proxy_z"] +
    df["num_predicted_proteins_z"] +
    df["proteins_per_10kb_z"]
)

# highest score = higher priority
df = df.sort_values("candidate_priority_score", ascending=False).reset_index(drop=True)
df["priority_rank"] = df.index + 1

# reorder useful columns
cols = [
    "priority_rank",
    "phage_id",
    "candidate_priority_score",
    "genome_length",
    "gc_content",
    "num_contigs",
    "num_predicted_proteins",
    "mean_protein_length_aa",
    "median_protein_length_aa",
    "total_protein_aa",
    "coding_density_proxy",
    "proteins_per_10kb"
]
df = df[cols]

df.to_csv(OUTPUT, index=False)

print("Saved ranked candidates to:", OUTPUT)
print()
print(df.head(15))
