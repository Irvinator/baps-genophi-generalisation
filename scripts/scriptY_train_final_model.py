#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
from catboost import CatBoostClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


DATA = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/final_merged_feature_table.csv")
OUTDIR = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/final_model_outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print(f"Reading: {DATA}")
    df = pd.read_csv(DATA)

    print("Input shape:", df.shape)
    print("Columns:", df.columns[:10].tolist())

    meta_cols = ["strain", "phage"]
    target_col = "interaction"
    feature_cols = [c for c in df.columns if c not in meta_cols + [target_col]]

    X = df[feature_cols]
    y = df[target_col].astype(int)

    print("Feature matrix shape:", X.shape)
    print("Positive labels:", int(y.sum()))
    print("Negative labels:", int((1 - y).sum()))

    model = CatBoostClassifier(
        iterations=500,
        depth=6,
        learning_rate=0.05,
        loss_function="Logloss",
        eval_metric="AUC",
        verbose=100,
        random_seed=42,
    )

    model.fit(X, y)

    probs = model.predict_proba(X)[:, 1]
    preds = (probs >= 0.5).astype(int)

    metrics = {
        "roc_auc": roc_auc_score(y, probs),
        "pr_auc": average_precision_score(y, probs),
        "accuracy": accuracy_score(y, preds),
        "precision": precision_score(y, preds, zero_division=0),
        "recall": recall_score(y, preds, zero_division=0),
        "f1": f1_score(y, preds, zero_division=0),
        "mcc": matthews_corrcoef(y, preds),
    }

    print("\nFinal metrics")
    for k, v in metrics.items():
        print(f"{k}: {v:.4f}")

    cm = confusion_matrix(y, preds)
    print("\nConfusion matrix:")
    print(cm)

    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(OUTDIR / "final_metrics.csv", index=False)

    pd.DataFrame(
        {
            "strain": df["strain"],
            "phage": df["phage"],
            "interaction": y,
            "pred_prob": probs,
            "pred_label": preds,
        }
    ).to_csv(OUTDIR / "final_predictions.csv", index=False)

    cm_df = pd.DataFrame(cm, index=["true_0", "true_1"], columns=["pred_0", "pred_1"])
    cm_df.to_csv(OUTDIR / "confusion_matrix.csv")

    importance = pd.DataFrame(
        {
            "feature": feature_cols,
            "importance": model.get_feature_importance(),
        }
    ).sort_values("importance", ascending=False)
    importance.to_csv(OUTDIR / "feature_importance.csv", index=False)

    model.save_model(str(OUTDIR / "final_catboost_model.cbm"))

    print(f"\nSaved outputs to: {OUTDIR}")


if __name__ == "__main__":
    main()
