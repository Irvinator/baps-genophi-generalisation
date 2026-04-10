#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split

DATA = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/final_merged_feature_table.csv")
IMPORTANCE = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/pilot_train_val_outputs/pilot_feature_importance.csv")
OUTDIR = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/pilot_train_val_outputs")

def main():
    df = pd.read_csv(DATA)
    imp = pd.read_csv(IMPORTANCE)

    meta_cols = ["strain", "phage"]
    target_col = "interaction"
    feature_cols = [c for c in df.columns if c not in meta_cols + [target_col]]

    X = df[feature_cols]
    y = df[target_col].astype(int)

    # Reproduce the exact pilot split
    X_train, X_val, y_train, y_val = train_test_split(
        X, y,
        test_size=0.2,
        random_state=42,
        stratify=y
    )

    print("Train shape:", X_train.shape)
    print("Val shape:", X_val.shape)

    # Binary presence across samples
    train_present = (X_train > 0).sum(axis=0)
    val_present   = (X_val > 0).sum(axis=0)

    coverage = pd.DataFrame({
        "feature": feature_cols,
        "train_count": train_present.values,
        "val_count": val_present.values
    })

    coverage["train_present"] = coverage["train_count"] > 0
    coverage["val_present"] = coverage["val_count"] > 0
    coverage["missing_in_val"] = coverage["train_present"] & (~coverage["val_present"])

    n_train_present = int(coverage["train_present"].sum())
    n_val_present = int(coverage["val_present"].sum())
    n_missing_in_val = int(coverage["missing_in_val"].sum())

    print("\nOverall feature coverage")
    print("Features present in train:", n_train_present)
    print("Features present in val:", n_val_present)
    print("Features present in train but absent in val:", n_missing_in_val)
    print("Fraction of train-present features absent in val:",
          round(n_missing_in_val / n_train_present, 4) if n_train_present else "NA")

    # Join with feature importance
    merged = imp.merge(coverage, on="feature", how="left")

    # Top 20 predictive features
    top20 = merged.head(20).copy()
    print("\nTop 20 feature coverage")
    print(top20[["feature", "importance", "train_count", "val_count", "missing_in_val"]].to_string(index=False))

    # Save
    coverage.to_csv(OUTDIR / "pilot_feature_coverage_all.csv", index=False)
    top20.to_csv(OUTDIR / "pilot_top20_feature_coverage.csv", index=False)

    print("\nSaved:")
    print(OUTDIR / "pilot_feature_coverage_all.csv")
    print(OUTDIR / "pilot_top20_feature_coverage.csv")

if __name__ == "__main__":
    main()
