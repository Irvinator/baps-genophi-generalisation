#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# =========================
# Paths
# =========================
BASE = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150")

importance_fp = BASE / "pilot_train_val_outputs" / "pilot_feature_importance.csv"

host_selected_fp = BASE / "train_feature_build" / "strain" / "features_gca_min2" / "selected_features.csv"
phage_selected_fp = BASE / "train_feature_build" / "phage_only_tmp" / "strain" / "features" / "selected_features.csv"

out_fp = BASE / "pilot_train_val_outputs" / "pilot_top_feature_mapping.csv"

# =========================
# Load data
# =========================
print("Loading feature importance...")
imp = pd.read_csv(importance_fp)

print("Loading selected feature mappings...")
host_sel = pd.read_csv(host_selected_fp)
phage_sel = pd.read_csv(phage_selected_fp)

# =========================
# Build summary tables
# =========================
def summarise_selected_features(df, prefix_name):
    """
    For each retained feature (e.g. sc_25817 / pc_124), summarise:
    - number of original cluster labels collapsed into it
    - first few example labels
    """
    grouped = (
        df.groupby("Feature")["Cluster_Label"]
        .apply(list)
        .reset_index()
    )

    grouped["n_original_proteins"] = grouped["Cluster_Label"].apply(len)
    grouped["example_cluster_label"] = grouped["Cluster_Label"].apply(lambda x: x[0] if len(x) > 0 else "")
    grouped["all_cluster_labels"] = grouped["Cluster_Label"].apply(lambda x: "; ".join(map(str, x[:10])))
    grouped["type"] = prefix_name

    return grouped[["Feature", "type", "n_original_proteins", "example_cluster_label", "all_cluster_labels"]]

host_summary = summarise_selected_features(host_sel, "host")
phage_summary = summarise_selected_features(phage_sel, "phage")

# Combine
feature_map = pd.concat([host_summary, phage_summary], ignore_index=True)

# =========================
# Keep only top 20 importances
# =========================
top_imp = imp.head(20).copy()

merged = top_imp.merge(feature_map, left_on="feature", right_on="Feature", how="left")

# Clean columns
merged = merged[[
    "feature",
    "importance",
    "type",
    "n_original_proteins",
    "example_cluster_label",
    "all_cluster_labels"
]]

# =========================
# Save
# =========================
merged.to_csv(out_fp, index=False)

print("\nSaved mapping table to:")
print(out_fp)

print("\nTop feature mapping preview:")
print(merged.head(20).to_string(index=False))
