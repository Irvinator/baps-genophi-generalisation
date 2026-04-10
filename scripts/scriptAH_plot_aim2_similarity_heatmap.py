#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

INPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/train_feature_build/phage/feature_table.csv"
OUTPUT_FIG = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_phage_similarity_heatmap.png"
OUTPUT_CSV = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_phage_similarity_matrix.csv"

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(INPUT)
pc_cols = [c for c in df.columns if c.startswith("pc_")]
X = df[pc_cols].astype(int)

print("Input shape:", df.shape)
print("Number of phages:", len(df))
print("Number of pc_ columns:", len(pc_cols))

# -----------------------------
# Compute pairwise Jaccard distance
# For binary data, Jaccard is a sensible choice
# similarity = 1 - distance
# -----------------------------
dist_vec = pdist(X.values, metric="jaccard")
sim_vec = 1 - dist_vec
sim_mat = squareform(sim_vec)

# Put diagonal back to 1
np.fill_diagonal(sim_mat, 1.0)

# -----------------------------
# Hierarchical clustering
# -----------------------------
Z = linkage(dist_vec, method="average")
order = leaves_list(Z)

ordered_sim = sim_mat[np.ix_(order, order)]
ordered_names = df["phage"].iloc[order].tolist()

# Save matrix
sim_df = pd.DataFrame(ordered_sim, index=ordered_names, columns=ordered_names)
sim_df.to_csv(OUTPUT_CSV)

print("Saved similarity matrix to:", OUTPUT_CSV)

# -----------------------------
# Plot heatmap
# -----------------------------
plt.figure(figsize=(8, 7))
plt.imshow(ordered_sim, aspect="auto")
plt.colorbar(label="Jaccard similarity")
plt.title("Hierarchical clustering of candidate phage feature space")
plt.xlabel("Phages (clustered)")
plt.ylabel("Phages (clustered)")
plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=300)

print("Saved heatmap to:", OUTPUT_FIG)
