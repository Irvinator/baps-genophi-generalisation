import pandas as pd

df = pd.read_csv("/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_feature_table.csv")

cols = [
    "genome_length",
    "gc_content",
    "num_predicted_proteins",
    "coding_density_proxy",
    "proteins_per_10kb"
]

summary = df[cols].agg(['mean', 'median', 'std', 'min', 'max']).T
summary["q1"] = df[cols].quantile(0.25)
summary["q3"] = df[cols].quantile(0.75)

summary = summary[["mean", "median", "std", "min", "q1", "q3", "max"]]

print(summary.round(4))
