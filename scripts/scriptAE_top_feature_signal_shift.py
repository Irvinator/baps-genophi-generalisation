#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split

DATA = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/final_merged_feature_table.csv")
IMPORTANCE = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/pilot_train_val_outputs/pilot_feature_importance.csv")
OUTDIR = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/pilot_train_val_outputs")

def frac_present(df, feature):
    return (df[feature] > 0).mean() if feature in df.columns else 0.0

def main():
    df = pd.read_csv(DATA)
    imp = pd.read_csv(IMPORTANCE)

    meta_cols = ["strain", "phage"]
    target_col = "interaction"
    feature_cols = [c for c in df.columns if c not in meta_cols + [target_col]]

    X = df[feature_cols]
    y = df[target_col].astype(int)

    # Reproduce pilot split
    X_train, X_val, y_train, y_val = train_test_split(
        X, y,
        test_size=0.2,
        random_state=42,
        stratify=y
    )

    train_df = X_train.copy()
    train_df["interaction"] = y_train.values
    val_df = X_val.copy()
    val_df["interaction"] = y_val.values

    train_pos = train_df[train_df["interaction"] == 1]
    train_neg = train_df[train_df["interaction"] == 0]
    val_pos = val_df[val_df["interaction"] == 1]
    val_neg = val_df[val_df["interaction"] == 0]

    rows = []
    top20 = imp.head(20)

    for _, r in top20.iterrows():
        feat = r["feature"]
        rows.append({
            "feature": feat,
            "importance": r["importance"],
            "train_pos_frac": frac_present(train_pos, feat),
            "train_neg_frac": frac_present(train_neg, feat),
            "val_pos_frac": frac_present(val_pos, feat),
            "val_neg_frac": frac_present(val_neg, feat),
        })

    out = pd.DataFrame(rows)
    out["train_signal"] = out["train_pos_frac"] - out["train_neg_frac"]
    out["val_signal"] = out["val_pos_frac"] - out["val_neg_frac"]
    out["signal_shift"] = out["val_signal"] - out["train_signal"]

    out.to_csv(OUTDIR / "pilot_top20_signal_shift.csv", index=False)

    print(out.to_string(index=False))
    print("\nSaved:")
    print(OUTDIR / "pilot_top20_signal_shift.csv")

if __name__ == "__main__":
    main()
