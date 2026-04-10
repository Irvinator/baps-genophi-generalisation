#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser(description="Build BAPS200 evaluation pairs (positives + sampled negatives).")
    ap.add_argument("--pos_pairs_tsv", required=True, help="TSV with header: host_accession, phage_contig, interaction=1")
    ap.add_argument("--hosts_list", required=True, help="Text file with 200 host accessions (one per line)")
    ap.add_argument("--contig_universe", required=True, help="Text file with all valid contigs (one per line), e.g. from FASTA")
    ap.add_argument("--neg_ratio", type=int, default=1, help="Negatives per positive (per host)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    hosts = set(pd.read_csv(args.hosts_list, header=None)[0].astype(str))
    contigs = pd.read_csv(args.contig_universe, header=None)[0].astype(str).tolist()
    contigs_set = set(contigs)

    df = pd.read_csv(args.pos_pairs_tsv, sep="\t")
    # normalize expected column names
    df.columns = [c.strip() for c in df.columns]

    if "host_accession" not in df.columns or "phage_contig" not in df.columns:
        raise ValueError(f"Expected columns host_accession, phage_contig. Got: {df.columns.tolist()}")

    df = df[df["host_accession"].astype(str).isin(hosts)].copy()
    df["host_accession"] = df["host_accession"].astype(str)
    df["phage_contig"] = df["phage_contig"].astype(str)

    # ensure positives are in contig universe
    df = df[df["phage_contig"].isin(contigs_set)].copy()
    df["interaction"] = 1
    df = df.drop_duplicates(["host_accession", "phage_contig"])

    print(f"[scriptF] Positives after filtering: {len(df)}")
    print(f"[scriptF] Hosts in positives: {df['host_accession'].nunique()} / {len(hosts)}")

    neg_rows = []
    for host, sub in df.groupby("host_accession"):
        pos = set(sub["phage_contig"].tolist())
        # negatives must be valid contigs and not in positives
        candidates = list(contigs_set - pos)
        n_neg = args.neg_ratio * len(pos)
        if n_neg <= 0 or len(candidates) == 0:
            continue
        chosen = rng.choice(candidates, size=min(n_neg, len(candidates)), replace=False)
        for c in chosen:
            neg_rows.append((host, c, 0))

    neg_df = pd.DataFrame(neg_rows, columns=["host_accession", "phage_contig", "interaction"])

    out = pd.concat([df[["host_accession","phage_contig","interaction"]], neg_df], ignore_index=True)
    out.to_csv(args.out_csv, index=False)

    print(f"[scriptF] Negatives created: {len(neg_df)}")
    print(f"[scriptF] Wrote: {args.out_csv}")
    print(out["interaction"].value_counts())

if __name__ == "__main__":
    main()
