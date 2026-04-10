import os
from pathlib import Path
import pandas as pd
import numpy as np

# ================================
# INPUT: change this if needed
# ================================
PHAGE_DIR = "/mnt/d/baps-genophi-generalisation/retrain_baseline/genophi_inputs/train_phage_AAs"

OUTPUT = "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/aim2_phage_feature_summary.csv"


def parse_fasta_lengths(filepath):
    lengths = []
    seq = ""

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    lengths.append(len(seq))
                    seq = ""
            else:
                seq += line.strip()

        if seq:
            lengths.append(len(seq))

    return lengths


def extract_features(phage_dir):
    rows = []

    files = list(Path(phage_dir).glob("*.faa"))
    print(f"Found {len(files)} phage files")

    for i, fp in enumerate(files):
        if i % 50 == 0:
            print(f"Processing {i}/{len(files)}")

        lengths = parse_fasta_lengths(fp)

        if len(lengths) == 0:
            continue

        lengths = np.array(lengths)

        row = {
            "phage": fp.stem,
            "n_proteins": len(lengths),
            "total_aa_length": lengths.sum(),
            "avg_protein_length": lengths.mean(),
            "std_protein_length": lengths.std(),
            "pct_long_proteins": np.mean(lengths > 300),
            "pct_short_proteins": np.mean(lengths < 100),
        }

        rows.append(row)

    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = extract_features(PHAGE_DIR)

    print("Final shape:", df.shape)
    print(df.head())

    df.to_csv(OUTPUT, index=False)
    print("Saved to:", OUTPUT)
