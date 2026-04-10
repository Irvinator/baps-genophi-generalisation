#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

INPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/train_feature_build/phage/feature_table.csv"
OUTPUT_FIG = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_phage_pca.png"
OUTPUT_CSV = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_phage_pca_coordinates.csv"

df = pd.read_csv(INPUT)

# Keep phage id + pc_ features only
pc_cols = [c for c in df.columns if c.startswith("pc_")]
X = df[pc_cols].copy()

print("Input shape:", df.shape)
print("Number of pc_ columns:", len(pc_cols))

# PCA
pca = PCA(n_components=2)
coords = pca.fit_transform(X)

pca_df = pd.DataFrame({
    "phage": df["phage"],
    "PC1": coords[:, 0],
    "PC2": coords[:, 1]
})

pca_df.to_csv(OUTPUT_CSV, index=False)

print("Explained variance ratio:", pca.explained_variance_ratio_)
print("Saved PCA coordinates to:", OUTPUT_CSV)

# Plot
plt.figure(figsize=(7.5, 6))
plt.scatter(pca_df["PC1"], pca_df["PC2"], alpha=0.8)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
plt.title("PCA of candidate phage feature space")
plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=300)

print("Saved PCA figure to:", OUTPUT_FIG)
