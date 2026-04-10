#!/usr/bin/env python3
"""
Compute min Mash distance per BAPS assembly accession (GCA/GCF) to a set of GenoPHI genomes.

Input: mash dist TSV with 5 columns:
  query_path<TAB>ref_path<TAB>distance<TAB>pvalue<TAB>shared_hashes

We:
  - extract BAPS accession (GCA_... or GCF_...) from ref_path
  - compute min(distance) across all query genomes for each BAPS accession
  - write CSV: baps_acc, min_dist_to_genophi_train, approx_ani

Example:
  python scripts/scriptE_compute_min_mash_distance.py \
    --mash_tsv outputs/mash_dist_genophiTrain_vs_baps200.abs.tsv \
    --out_csv outputs/baps200_min_dist_to_genophi_train.csv
"""

from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path
from math import inf


ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def extract_acc_from_path(s: str) -> str | None:
    m = ACC_RE.search(s)
    return m.group(1) if m else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mash_tsv", required=True, type=Path, help="Mash dist TSV (5 columns)")
    ap.add_argument("--out_csv", required=True, type=Path, help="Output CSV file")
    ap.add_argument("--expect_baps_list", type=Path, default=None,
                    help="Optional file of allowed BAPS accessions (one per line); useful to ignore noise")
    args = ap.parse_args()

    allowed: set[str] | None = None
    if args.expect_baps_list:
        allowed = set(x.strip() for x in args.expect_baps_list.read_text().splitlines() if x.strip())

    # min distance per baps acc
    mins: dict[str, float] = {}
    rows = 0
    skipped_no_acc = 0
    skipped_not_allowed = 0

    with args.mash_tsv.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            rows += 1

            ref_path = row[1]
            dist_str = row[2]

            baps_acc = extract_acc_from_path(ref_path)
            if not baps_acc:
                skipped_no_acc += 1
                continue
            if allowed is not None and baps_acc not in allowed:
                skipped_not_allowed += 1
                continue

            try:
                d = float(dist_str)
            except ValueError:
                continue

            cur = mins.get(baps_acc, inf)
            if d < cur:
                mins[baps_acc] = d

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)

    # write
    with args.out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["baps_acc", "min_dist_to_genophi_train", "approx_ani"])
        for acc in sorted(mins.keys()):
            d = mins[acc]
            w.writerow([acc, f"{d:.6f}", f"{1.0 - d:.6f}"])

    # report
    print(f"[scriptE] Read rows: {rows:,}")
    print(f"[scriptE] Unique BAPS accessions: {len(mins):,}")
    if allowed is not None:
        print(f"[scriptE] Allowed list size: {len(allowed):,}")
        print(f"[scriptE] Skipped (not in allowed list): {skipped_not_allowed:,}")
    print(f"[scriptE] Skipped (no accession found in ref path): {skipped_no_acc:,}")
    print(f"[scriptE] Wrote: {args.out_csv}")

    # quick stats
    vals = list(mins.values())
    if vals:
        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        mean = sum(vals_sorted) / n
        p50 = vals_sorted[n // 2]
        p25 = vals_sorted[int(0.25 * (n - 1))]
        p75 = vals_sorted[int(0.75 * (n - 1))]
        print("[scriptE] Distance summary:")
        print(f"  min={vals_sorted[0]:.6f}  p25={p25:.6f}  p50={p50:.6f}  p75={p75:.6f}  max={vals_sorted[-1]:.6f}  mean={mean:.6f}")


if __name__ == "__main__":
    main()
