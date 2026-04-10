import pandas as pd
import matplotlib.pyplot as plt

# Load GenoPHI phage feature table
df = pd.read_csv("/home/irvin/projects/GenoPHI/results/baps_assign_phage/phage_combined_feature_table_eval.csv")

# Identify protein family columns
cluster_cols = [c for c in df.columns if c.startswith("pc_")]

# Count number of protein families present per phage
df["families_present"] = (df[cluster_cols] > 0).sum(axis=1)

# Plot histogram
plt.figure(figsize=(8,6))
plt.hist(df["families_present"], bins=20)

plt.xlabel("Number of GenoPHI protein families detected per BAPS phage")
plt.ylabel("Number of phages")
plt.title("Distribution of GenoPHI protein families present in BAPS phage genomes")

plt.tight_layout()

plt.savefig("phage_family_overlap_correct.png", dpi=300)

print("Saved phage_family_overlap_correct.png")
