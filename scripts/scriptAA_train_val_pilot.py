#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from catboost import CatBoostClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split


DATA = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/final_merged_feature_table.csv")
OUTDIR = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/pilot_train_val_outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)


def plot_confusion_matrix(cm, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 5.6), dpi=300)
    row_sums = cm.sum(axis=1, keepdims=True)
    pct = cm / row_sums * 100

    labels = [
        [f"{cm[0,0]}\n({pct[0,0]:.1f}%)", f"{cm[0,1]}\n({pct[0,1]:.1f}%)"],
        [f"{cm[1,0]}\n({pct[1,0]:.1f}%)", f"{cm[1,1]}\n({pct[1,1]:.1f}%)"],
    ]

    im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Predicted\nNon-interaction", "Predicted\nInteraction"], fontsize=11)
    ax.set_yticklabels(["Actual\nNon-interaction", "Actual\nInteraction"], fontsize=11)

    threshold = cm.max() / 2.0
    for i in range(2):
        for j in range(2):
            ax.text(
                j, i, labels[i][j],
                ha="center", va="center",
                color="white" if cm[i, j] > threshold else "black",
                fontsize=12,
            )

    ax.set_title("Pilot train/validation confusion matrix", fontsize=14, pad=12)
    ax.set_xlabel("Predicted class", fontsize=12)
    ax.set_ylabel("True class", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def plot_roc(y_true, y_prob, out_png: Path) -> float:
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    auc = roc_auc_score(y_true, y_prob)

    plt.figure(figsize=(6.2, 5.0), dpi=300)
    plt.plot(fpr, tpr, label=f"ROC-AUC = {auc:.3f}")
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("Pilot validation ROC curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    return auc


def plot_pr(y_true, y_prob, out_png: Path) -> float:
    precision, recall, _ = precision_recall_curve(y_true, y_prob)
    ap = average_precision_score(y_true, y_prob)

    plt.figure(figsize=(6.2, 5.0), dpi=300)
    plt.plot(recall, precision, label=f"PR-AUC = {ap:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Pilot validation precision-recall curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    return ap


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

    X_train, X_val, y_train, y_val, meta_train, meta_val = train_test_split(
        X,
        y,
        df[meta_cols],
        test_size=0.2,
        random_state=42,
        stratify=y,
    )

    print("Train shape:", X_train.shape)
    print("Validation shape:", X_val.shape)
    print("Validation positives:", int(y_val.sum()))
    print("Validation negatives:", int((1 - y_val).sum()))

    model = CatBoostClassifier(
        iterations=500,
        depth=6,
        learning_rate=0.05,
        loss_function="Logloss",
        eval_metric="AUC",
        verbose=100,
        random_seed=42,
    )

    model.fit(X_train, y_train)

    val_probs = model.predict_proba(X_val)[:, 1]
    val_preds = (val_probs >= 0.5).astype(int)

    metrics = {
        "roc_auc": roc_auc_score(y_val, val_probs),
        "pr_auc": average_precision_score(y_val, val_probs),
        "accuracy": accuracy_score(y_val, val_preds),
        "precision": precision_score(y_val, val_preds, zero_division=0),
        "recall": recall_score(y_val, val_preds, zero_division=0),
        "f1": f1_score(y_val, val_preds, zero_division=0),
        "mcc": matthews_corrcoef(y_val, val_preds),
    }

    print("\nPilot validation metrics")
    for k, v in metrics.items():
        print(f"{k}: {v:.4f}")

    cm = confusion_matrix(y_val, val_preds)
    print("\nValidation confusion matrix:")
    print(cm)

    pd.DataFrame([metrics]).to_csv(OUTDIR / "pilot_validation_metrics.csv", index=False)
    pd.DataFrame(cm, index=["true_0", "true_1"], columns=["pred_0", "pred_1"]).to_csv(
        OUTDIR / "pilot_validation_confusion_matrix.csv"
    )

    pred_df = pd.DataFrame({
        "strain": meta_val["strain"].values,
        "phage": meta_val["phage"].values,
        "interaction": y_val.values,
        "pred_prob": val_probs,
        "pred_label": val_preds,
    })
    pred_df.to_csv(OUTDIR / "pilot_validation_predictions.csv", index=False)

    importance = pd.DataFrame({
        "feature": feature_cols,
        "importance": model.get_feature_importance(),
    }).sort_values("importance", ascending=False)
    importance.to_csv(OUTDIR / "pilot_feature_importance.csv", index=False)

    plot_confusion_matrix(cm, OUTDIR / "pilot_validation_confusion_matrix.png")
    plot_roc(y_val, val_probs, OUTDIR / "pilot_validation_roc_curve.png")
    plot_pr(y_val, val_probs, OUTDIR / "pilot_validation_pr_curve.png")

    model.save_model(str(OUTDIR / "pilot_catboost_model.cbm"))
    print(f"\nSaved outputs to: {OUTDIR}")


if __name__ == "__main__":
    main()
