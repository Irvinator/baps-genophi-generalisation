#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


INFILE = Path("/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/pilot_train_val_outputs/pilot_feature_importance.csv")
OUTDIR = Path("/mnt/d/baps-genophi-generalisation/report_figures")
OUTDIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    df = pd.read_csv(INFILE)

    # top 20
    top = df.sort_values("importance", ascending=False).head(20).copy()

    # add feature type
    top["feature_type"] = top["feature"].apply(
        lambda x: "Host (sc_)" if str(x).startswith("sc_") else "Phage (pc_)"
    )

    # reverse for plotting
    top = top.iloc[::-1]

    plt.figure(figsize=(9, 7), dpi=300)

    # Use default matplotlib colors automatically by grouping
    host = top[top["feature_type"] == "Host (sc_)"]
    phage = top[top["feature_type"] == "Phage (pc_)"]

    plt.barh(host["feature"], host["importance"], label="Host features")
    plt.barh(phage["feature"], phage["importance"], label="Phage features")

    plt.xlabel("Feature importance")
    plt.ylabel("Feature")
    plt.title("Top 20 CatBoost feature importances from the pilot retraining model")
    plt.legend()
    plt.tight_layout()

    out_png = OUTDIR / "pilot_top20_feature_importance.png"
    out_pdf = OUTDIR / "pilot_top20_feature_importance.pdf"

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"Saved to: {out_png}")
    print(f"Saved to: {out_pdf}")


if __name__ == "__main__":
    main()
