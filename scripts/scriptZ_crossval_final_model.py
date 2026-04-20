#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
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
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "reduced" / "final_merged_feature_table.csv"
OUTDIR = ROOT / "outputs" / "final_model_outputs_cv"
OUTDIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print(f"Reading: {DATA}")
    df = pd.read_csv(DATA)

    meta_cols = ["strain", "phage"]
    target_col = "interaction"
    feature_cols = [c for c in df.columns if c not in meta_cols + [target_col]]

    X = df[feature_cols]
    y = df[target_col].astype(int)

    print("Input shape:", df.shape)
    print("Feature matrix shape:", X.shape)
    print("Positive labels:", int(y.sum()))
    print("Negative labels:", int((1 - y).sum()))

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    fold_metrics = []
    all_preds = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), start=1):
        print(f"\n=== Fold {fold} ===")

        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        model = CatBoostClassifier(
            iterations=500,
            depth=6,
            learning_rate=0.05,
            loss_function="Logloss",
            eval_metric="AUC",
            verbose=0,
            random_seed=42,
        )

        model.fit(X_train, y_train)

        probs = model.predict_proba(X_test)[:, 1]
        preds = (probs >= 0.5).astype(int)

        metrics = {
            "fold": fold,
            "n_test": len(test_idx),
            "roc_auc": roc_auc_score(y_test, probs),
            "pr_auc": average_precision_score(y_test, probs),
            "accuracy": accuracy_score(y_test, preds),
            "precision": precision_score(y_test, preds, zero_division=0),
            "recall": recall_score(y_test, preds, zero_division=0),
            "f1": f1_score(y_test, preds, zero_division=0),
            "mcc": matthews_corrcoef(y_test, preds),
        }
        fold_metrics.append(metrics)

        print({k: round(v, 4) if isinstance(v, float) else v for k, v in metrics.items()})

        fold_pred_df = pd.DataFrame({
            "fold": fold,
            "strain": df.iloc[test_idx]["strain"].values,
            "phage": df.iloc[test_idx]["phage"].values,
            "interaction": y_test.values,
            "pred_prob": probs,
            "pred_label": preds,
        })
        all_preds.append(fold_pred_df)

    metrics_df = pd.DataFrame(fold_metrics)
    preds_df = pd.concat(all_preds, ignore_index=True)

    summary = metrics_df.drop(columns=["fold", "n_test"]).agg(["mean", "std"]).T
    print("\n=== CV Summary ===")
    print(summary)

    # overall confusion matrix from out-of-fold predictions
    cm = confusion_matrix(preds_df["interaction"], preds_df["pred_label"])
    print("\nOut-of-fold confusion matrix:")
    print(cm)

    metrics_df.to_csv(OUTDIR / "cv_fold_metrics.csv", index=False)
    summary.to_csv(OUTDIR / "cv_summary_metrics.csv")
    preds_df.to_csv(OUTDIR / "cv_predictions.csv", index=False)
    pd.DataFrame(cm, index=["true_0", "true_1"], columns=["pred_0", "pred_1"]).to_csv(
        OUTDIR / "cv_confusion_matrix.csv"
    )

    print(f"\nSaved outputs to: {OUTDIR}")


if __name__ == "__main__":
    main()
