#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from pathlib import Path


def extract_gca(sample: str):
    """
    Extract NCBI assembly accession from strings like:
      "GCA_018015695.1_-_PDT000998088.1"
    Returns e.g. "GCA_018015695.1"
    """
    m = re.search(r"(GCA_\d+\.\d+)", str(sample))
    return m.group(1) if m else None


def main():
    parser = argparse.ArgumentParser(
        description="Extract E. coli positive host–phage contig pairs from BAPS annotations table."
    )
    parser.add_argument(
        "--baps_annotations",
        required=True,
        help="Path to BAPS annotations TSV (e.g. current_BAPS_anns.tsv OR your pre-filtered baps_ecoli.tsv).",
    )
    parser.add_argument(
        "--out_pairs",
        required=True,
        help="Output TSV path for positive pairs (host_accession, phage_contig, interaction).",
    )
    parser.add_argument(
        "--out_hosts",
        required=True,
        help="Output TXT path for unique host accessions (one per line).",
    )
    parser.add_argument(
        "--ecoli_only",
        action="store_true",
        help="If set, filter rows to Escherichia coli using available taxonomy columns (best-effort).",
    )
    args = parser.parse_args()

    in_path = Path(args.baps_annotations).expanduser()
    out_pairs = Path(args.out_pairs).expanduser()
    out_hosts = Path(args.out_hosts).expanduser()

    if not in_path.exists():
        raise FileNotFoundError(f"Input TSV not found: {in_path}")

    out_pairs.parent.mkdir(parents=True, exist_ok=True)
    out_hosts.parent.mkdir(parents=True, exist_ok=True)

    # Read TSV robustly. Some BAPS exports have a leading tab in the header.
    # pandas handles it, but we strip accidental whitespace from column names.
    df = pd.read_csv(in_path, sep="\t", low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    # Validate expected columns
    required_cols = {"sample", "contig"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns {sorted(missing)} in {in_path}. "
            f"Found columns: {list(df.columns)[:30]} ..."
        )

    # Optional: filter to E. coli if the TSV is not already E. coli-only.
    # We do best-effort depending on which column(s) exist.
    if args.ecoli_only:
        ecoli_mask = None

        # Try common column names that might contain taxonomy/organism strings
        candidate_cols = [
            "org", "organism", "organism_name", "host", "host_name",
            "host_taxonomy", "taxonomy", "lineage"
        ]
        for col in candidate_cols:
            if col in df.columns:
                m = df[col].astype(str).str.contains(r"Escherichia coli", case=False, na=False)
                ecoli_mask = m if ecoli_mask is None else (ecoli_mask | m)

        if ecoli_mask is None:
            print("[scriptA] --ecoli_only requested but no taxonomy-like columns found; skipping filter.")
        else:
            before = len(df)
            df = df[ecoli_mask].copy()
            print(f"[scriptA] Filtered to E. coli: {len(df)} / {before} rows kept")

    # Extract host accession from sample field
    df["host_accession"] = df["sample"].apply(extract_gca)

    # Build positive pairs table
    pos = df[["host_accession", "contig"]].dropna().drop_duplicates()
    pos = pos.rename(columns={"contig": "phage_contig"})
    pos["interaction"] = 1

    # Write outputs
    pos.to_csv(out_pairs, sep="\t", index=False)
    pos["host_accession"].drop_duplicates().to_csv(out_hosts, index=False, header=False)

    print(f"[scriptA] Wrote positives: {out_pairs} (rows={len(pos)})")
    print(f"[scriptA] Wrote accessions: {out_hosts} (unique={pos['host_accession'].nunique()})")
    print("[scriptA] Example rows:")
    print(pos.head(5).to_string(index=False))


if __name__ == "__main__":
    main()


