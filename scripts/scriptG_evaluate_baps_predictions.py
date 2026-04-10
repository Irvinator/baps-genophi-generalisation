#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
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


def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Evaluate GenoPHI predictions against BAPS truth labels."
    )
    ap.add_argument(
        "--truth_csv",
        required=True,
        help="CSV with truth labels, e.g. columns: strain, phage, interaction",
    )
    ap.add_argument(
        "--pred_csv",
        required=True,
        help="CSV with GenoPHI predictions, e.g. columns: strain, phage, Confidence, Final_Prediction",
    )
    ap.add_argument(
        "--out_csv",
        required=True,
        help="Where to save merged truth + predictions CSV",
    )
    args = ap.parse_args()

    truth = pd.read_csv(args.truth_csv)
    pred = pd.read_csv(args.pred_csv)

    # Resolve truth columns
    strain_t = find_col(truth, ["strain"])
    phage_t = find_col(truth, ["phage"])
    y_t = find_col(truth, ["interaction", "label", "y_true"])

    if not all([strain_t, phage_t, y_t]):
        raise ValueError(
            f"Truth file missing required columns. Got: {truth.columns.tolist()}"
        )

    # Resolve prediction columns
    strain_p = find_col(pred, ["strain"])
    phage_p = find_col(pred, ["phage"])
    score_p = find_col(
        pred,
        [
            "Confidence",
            "probability",
            "prediction_probability",
            "score",
            "pred_score",
            "predicted_probability",
            "positive_class_probability",
            "median_prediction",
            "prediction",
            "pred",
        ],
    )
    predclass_p = find_col(
        pred,
        [
            "Final_Prediction",
            "predicted_class",
            "class",
        ],
    )

    if not all([strain_p, phage_p, score_p]):
        raise ValueError(
            f"Prediction file missing required columns. Got: {pred.columns.tolist()}"
        )

    # Standardize truth
    truth = truth[[strain_t, phage_t, y_t]].copy()
    truth.columns = ["strain", "phage", "interaction"]

    # Standardize prediction
    keep_cols = [strain_p, phage_p, score_p]
    if predclass_p:
        keep_cols.append(predclass_p)

    pred = pred[keep_cols].copy()

    rename_map = {
        strain_p: "strain",
        phage_p: "phage",
        score_p: "score",
    }
    if predclass_p:
        rename_map[predclass_p] = "pred_class"

    pred = pred.rename(columns=rename_map)

    # Normalize types
    truth["strain"] = truth["strain"].astype(str)
    truth["phage"] = truth["phage"].astype(str)
    truth["interaction"] = truth["interaction"].astype(int)

    pred["strain"] = pred["strain"].astype(str)
    pred["phage"] = pred["phage"].astype(str)
    pred["score"] = pred["score"].astype(float)

    if "pred_class" in pred.columns:
        # Handle common string/bool/int class encodings
        pred["pred_class"] = pred["pred_class"].replace(
            {
                "True": 1,
                "False": 0,
                "true": 1,
                "false": 0,
                True: 1,
                False: 0,
            }
        )
        pred["pred_class"] = pred["pred_class"].astype(int)

    merged = truth.merge(pred, on=["strain", "phage"], how="left")

    missing = merged["score"].isna().sum()
    print(f"Merged rows: {len(merged)}")
    print(f"Missing predictions: {missing}")

    if missing > 0:
        print("Dropping rows with missing predictions for metric calculation.")
        merged = merged.dropna(subset=["score"]).copy()

    if len(merged) == 0:
        raise ValueError("No rows left after merging/dropping missing predictions.")

    y_true = merged["interaction"].astype(int)
    y_score = merged["score"].astype(float)

    if "pred_class" in merged.columns:
        y_pred = merged["pred_class"].astype(int)
    else:
        y_pred = (y_score >= 0.5).astype(int)

    metrics = {
        "n_pairs": len(merged),
        "roc_auc": roc_auc_score(y_true, y_score),
        "pr_auc": average_precision_score(y_true, y_score),
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred),
    }

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    metrics.update(
        {
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn),
            "tp": int(tp),
        }
    )

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)

    print("\nMetrics:")
    for k, v in metrics.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
