import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import random


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pos_pairs_tsv", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--n_hosts", type=int, default=200)
    parser.add_argument("--max_pos_per_host", type=int, default=10)
    parser.add_argument("--neg_ratio", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    # Load positives
    df = pd.read_csv(args.pos_pairs_tsv, sep="\t")
    print(f"[scriptB] Loaded positives: {len(df):,}")

    print(f"[scriptB] Unique hosts: {df['host_accession'].nunique():,}")
    print(f"[scriptB] Unique phage contigs: {df['phage_contig'].nunique():,}")

    # Sample hosts
    hosts = df["host_accession"].unique()
    sampled_hosts = np.random.choice(
        hosts,
        size=min(args.n_hosts, len(hosts)),
        replace=False
    )

    df_hosts = df[df["host_accession"].isin(sampled_hosts)].copy()

    # Cap positives per host
    df_pos = (
        df_hosts.groupby("host_accession", group_keys=False)
        .apply(lambda x: x.sample(
            n=min(len(x), args.max_pos_per_host),
            random_state=args.seed
        ))
        .reset_index(drop=True)
    )

    df_pos["interaction"] = 1

    print(f"[scriptB] Positives kept after cap: {len(df_pos):,}")

    # Build negatives
    all_contigs = df["phage_contig"].unique()
    neg_rows = []

    for host in sampled_hosts:
        host_pos = df_pos[df_pos["host_accession"] == host]["phage_contig"].values
        candidate_neg = list(set(all_contigs) - set(host_pos))

        n_neg = args.neg_ratio * len(host_pos)

        if len(candidate_neg) == 0 or n_neg == 0:
            continue

        chosen_neg = np.random.choice(
            candidate_neg,
            size=min(n_neg, len(candidate_neg)),
            replace=False
        )

        for contig in chosen_neg:
            neg_rows.append((host, contig, 0))

    neg_df = pd.DataFrame(
        neg_rows,
        columns=["host_accession", "phage_contig", "interaction"]
    )

    print(f"[scriptB] Negatives created: {len(neg_df):,}")

    final_df = pd.concat([df_pos, neg_df], ignore_index=True)

    final_df.to_csv(args.out_csv, index=False)

    print(f"[scriptB] Wrote dataset: {args.out_csv}")
    print(final_df["interaction"].value_counts())


if __name__ == "__main__":
    main()
