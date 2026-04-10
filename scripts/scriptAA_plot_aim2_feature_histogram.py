import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Input / output
# -----------------------------
INPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_rebootability_candidate_features.csv"
OUTPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_phage_feature_histogram.png"

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(INPUT)

print("Loaded:", df.shape)
print(df.head())

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(8, 5))
plt.hist(df["n_pc_features_present"], bins=20, edgecolor="black")
plt.xlabel("Number of phage protein-family features present")
plt.ylabel("Number of phages")
plt.title("Distribution of candidate phage feature richness")
plt.tight_layout()

plt.savefig(OUTPUT, dpi=300)
print("Saved to:", OUTPUT)
