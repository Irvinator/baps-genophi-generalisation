#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def load_contig_list(path: Path) -> np.ndarray:
    contigs = pd.read_csv(path, header=None)[0].astype(str).values
    return contigs


def main():
    ap = argparse.ArgumentParser(description="Build balanced (pos/neg) evaluation pairs from BAPS positives.")
    ap.add_argument("--pos_pairs_tsv", required=True, help="TSV with columns: host_accession, phage_contig, interaction (1)")
    ap.add_argument("--allowed_phages", required=True, help="Text file: one phage contig per line (e.g. from BAPS FASTA)")
    ap.add_argument("--out_csv", required=True, help="Output CSV path")
    ap.add_argument("--n_hosts", type=int, default=200)
    ap.add_argument("--max_pos_per_host", type=int, default=10)
    ap.add_argument("--neg_per_pos", type=int, default=1, help="Number of negatives per positive")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out_hosts_txt", default=None, help="Optional: save sampled hosts list here")
    args = ap.parse_args()

    pos_path = Path(args.pos_pairs_tsv)
    phage_list_path = Path(args.allowed_phages)
    out_csv = Path(args.out_csv)

    if not pos_path.exists():
        raise FileNotFoundError(f"pos_pairs_tsv not found: {pos_path}")
    if not phage_list_path.exists():
        raise FileNotFoundError(f"allowed_phages not found: {phage_list_path}")

    rng = np.random.default_rng(args.seed)

    # Read positives
    df = pd.read_csv(pos_path, sep="\t")
    required = {"host_accession", "phage_contig"}
    if not required.issubset(df.columns):
        raise ValueError(f"Expected columns {required} in {pos_path}, got: {list(df.columns)}")

    df["host_accession"] = df["host_accession"].astype(str)
    df["phage_contig"] = df["phage_contig"].astype(str)

    # Ensure interaction column exists + is 1 for positives
    if "interaction" not in df.columns:
        df["interaction"] = 1
    else:
        df["interaction"] = 1

    # Allowed phage universe (e.g. all contigs in FASTA)
    all_phages = load_contig_list(phage_list_path)
    all_phage_set = set(all_phages.tolist())

    # Filter any positives that somehow aren’t in allowed list (extra safety)
    before = len(df)
    df = df[df["phage_contig"].isin(all_phage_set)].copy()
    after = len(df)

    print(f"[scriptF] Loaded positives: {before:,}")
    print(f"[scriptF] Positives in allowed phage list: {after:,}")
    print(f"[scriptF] Unique hosts in positives: {df['host_accession'].nunique():,}")
    print(f"[scriptF] Allowed phage contigs: {len(all_phages):,}")

    # Sample hosts
    hosts = df["host_accession"].unique()
    n_hosts = min(args.n_hosts, len(hosts))
    sampled_hosts = rng.choice(hosts, size=n_hosts, replace=False)

    if args.out_hosts_txt:
        out_hosts = Path(args.out_hosts_txt)
        out_hosts.parent.mkdir(parents=True, exist_ok=True)
        pd.Series(sampled_hosts).to_csv(out_hosts, index=False, header=False)
        print(f"[scriptF] Wrote sampled host list: {out_hosts} (n={len(sampled_hosts)})")

    # Restrict positives to sampled hosts
    df = df[df["host_accession"].isin(sampled_hosts)].copy()

    # Cap positives per host
    df_pos = (
        df.groupby("host_accession", group_keys=False)
          .apply(lambda x: x.sample(n=min(len(x), args.max_pos_per_host), random_state=args.seed))
          .reset_index(drop=True)
    )
    df_pos["interaction"] = 1

    print(f"[scriptF] Positives kept after cap: {len(df_pos):,}")

    # Build negatives: for each host, sample phages not in its positive set
    neg_rows = []
    # Precompute host -> set(positive phages)
    host_to_pos = df_pos.groupby("host_accession")["phage_contig"].apply(set).to_dict()

    for host, pos_set in host_to_pos.items():
        n_pos = len(pos_set)
        n_neg = args.neg_per_pos * n_pos
        if n_neg <= 0:
            continue

        # Candidate negatives = allowed phages minus positives
        # For speed: use boolean mask on numpy array
        # (still fine at this scale)
        candidates = [p for p in all_phages if p not in pos_set]
        if len(candidates) == 0:
            continue

        chosen = rng.choice(candidates, size=min(n_neg, len(candidates)), replace=False)
        for phage in chosen:
            neg_rows.append((host, phage, 0))

    df_neg = pd.DataFrame(neg_rows, columns=["host_accession", "phage_contig", "interaction"])
    print(f"[scriptF] Negatives created: {len(df_neg):,}")

    out = pd.concat([df_pos[["host_accession","phage_contig","interaction"]], df_neg], ignore_index=True)

    # Shuffle rows so positives/negatives are mixed
    out = out.sample(frac=1, random_state=args.seed).reset_index(drop=True)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)

    print(f"[scriptF] Wrote eval dataset: {out_csv}")
    print(out["interaction"].value_counts())


if __name__ == "__main__":
    main()
