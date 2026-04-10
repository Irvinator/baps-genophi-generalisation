#!/usr/bin/env python3
import pandas as pd

INPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_rebootability_candidate_features.csv"
OUTPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_feature_summary_table.csv"

df = pd.read_csv(INPUT)

cols = [
    "n_pc_features_present",
    "total_pc_feature_count",
    "pc_feature_density"
]

summary = pd.DataFrame({
    "mean": df[cols].mean(),
    "median": df[cols].median(),
    "std": df[cols].std(),
    "min": df[cols].min(),
    "q1": df[cols].quantile(0.25),
    "q3": df[cols].quantile(0.75),
    "max": df[cols].max()
})

summary = summary.round(4)
summary.to_csv(OUTPUT)

print("Saved to:", OUTPUT)
print()
print(summary)
