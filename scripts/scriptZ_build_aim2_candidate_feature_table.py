import pandas as pd
import numpy as np

# ==========================
# INPUT / OUTPUT
# ==========================
INPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/train_feature_build/phage/feature_table.csv"
OUTPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_rebootability_candidate_features.csv"

# ==========================
# LOAD
# ==========================
df = pd.read_csv(INPUT)

print("Input shape:", df.shape)
print("Columns:", df.columns[:10].tolist())

# Identify phage feature columns
pc_cols = [c for c in df.columns if c.startswith("pc_")]

print("Number of phage feature columns:", len(pc_cols))

# ==========================
# BUILD SUMMARY FEATURES
# ==========================
out = pd.DataFrame()
out["phage"] = df["phage"]

# Number of present protein-family features
out["n_pc_features_present"] = (df[pc_cols] > 0).sum(axis=1)

# Total feature count (same as above if binary, but keep explicit)
out["total_pc_feature_count"] = df[pc_cols].sum(axis=1)

# Density across full feature space
out["pc_feature_density"] = out["n_pc_features_present"] / len(pc_cols)

# Optional simple ranking variables
out["pc_feature_richness_rank"] = out["n_pc_features_present"].rank(method="average", ascending=False)

# ==========================
# SAVE
# ==========================
out.to_csv(OUTPUT, index=False)

print("\nSaved candidate Aim 2 feature table to:")
print(OUTPUT)

print("\nPreview:")
print(out.head(10))

print("\nSummary:")
print(out.describe(include="all"))
