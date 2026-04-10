import pandas as pd
import matplotlib.pyplot as plt

# BAPS phage feature table
baps = pd.read_csv("/home/irvin/projects/GenoPHI/results/baps_assign_phage/phage_combined_feature_table_eval.csv")
baps_cols = [c for c in baps.columns if c.startswith("pc_")]
baps["families_present"] = (baps[baps_cols] > 0).sum(axis=1)

# GenoPHI training phage feature table
train = pd.read_csv("/home/irvin/projects/GenoPHI/results/ecoli_pf_quick/phage/features/feature_table.csv")
train_cols = [c for c in train.columns if c.startswith("pc_")]
train["families_present"] = (train[train_cols] > 0).sum(axis=1)

plt.figure(figsize=(8, 6))

bins = range(
    0,
    int(max(train["families_present"].max(), baps["families_present"].max())) + 2
)

plt.hist(
    train["families_present"],
    bins=bins,
    alpha=0.6,
    label="GenoPHI training phages"
)

plt.hist(
    baps["families_present"],
    bins=bins,
    alpha=0.6,
    label="BAPS test phages"
)

plt.xlabel("Number of GenoPHI protein families detected per phage")
plt.ylabel("Number of phages")
plt.title("Comparison of protein family content in GenoPHI training and BAPS test phages")
plt.legend()
plt.tight_layout()
plt.savefig("training_vs_baps_phage_family_overlap.png", dpi=300)
print("Saved training_vs_baps_phage_family_overlap.png")
