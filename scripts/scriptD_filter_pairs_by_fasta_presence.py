#!/usr/bin/env python3
"""
Filter BAPS host–phage positive pairs to those whose phage contig IDs
exist in the BAPS lytic phage FASTA headers.

Input pairs TSV expected columns:
  host_accession<TAB>phage_contig<TAB>interaction

Two options for contig presence:
  (A) Provide a precomputed contig list (one contig per line)
  (B) Provide the gz FASTA and we will extract contigs from headers

Example:
  python scripts/scriptD_filter_pairs_by_fasta_presence.py \
    --pairs_tsv outputs/baps_ecoli_pos_pairs.tsv \
    --contig_list outputs/baps_all_phage_contigs_from_fasta.txt \
    --out_tsv outputs/baps_ecoli_pos_pairs_in_fasta.tsv
"""

from __future__ import annotations
import argparse
import gzip
import re
from pathlib import Path


HEADER_CONTIG_RE = re.compile(r"__([A-Z0-9]+\.\d+)__")  # e.g. __NAFV01000136.1__


def extract_contig_from_header(header_line: str) -> str | None:
    """
    Extract contig accession from a FASTA header line.
    Expected format contains: __<CONTIG>__ e.g. __NAFV01000136.1__
    """
    m = HEADER_CONTIG_RE.search(header_line)
    return m.group(1) if m else None


def load_contigs_from_list(path: Path) -> set[str]:
    contigs: set[str] = set()
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            c = line.strip()
            if c:
                contigs.add(c)
    return contigs


def load_contigs_from_fasta_gz(fasta_gz: Path, max_headers: int | None = None) -> set[str]:
    """
    Stream through gz FASTA and collect contigs from header lines.
    If max_headers is set, stops after that many header lines (debug only).
    """
    contigs: set[str] = set()
    headers_seen = 0
    with gzip.open(fasta_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line:
                continue
            if line[0] == ">":
                headers_seen += 1
                c = extract_contig_from_header(line)
                if c:
                    contigs.add(c)
                if max_headers is not None and headers_seen >= max_headers:
                    break
    return contigs


def filter_pairs(pairs_tsv: Path, contigs_ok: set[str], out_tsv: Path) -> tuple[int, int]:
    """
    Write filtered TSV and return (rows_in, rows_out) excluding header.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    rows_in = 0
    rows_out = 0

    with pairs_tsv.open("r", encoding="utf-8") as fin, out_tsv.open("w", encoding="utf-8") as fout:
        header = fin.readline()
        if not header:
            raise RuntimeError(f"Empty file: {pairs_tsv}")
        fout.write(header)

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            rows_in += 1
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            phage_contig = parts[1]
            if phage_contig in contigs_ok:
                fout.write(line + "\n")
                rows_out += 1

    return rows_in, rows_out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs_tsv", required=True, type=Path, help="BAPS host-phage pairs TSV")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--contig_list", type=Path, help="One contig per line (fast mode)")
    group.add_argument("--fasta_gz", type=Path, help="BAPS_lytic_phage_0.1.fasta.gz (direct mode)")
    ap.add_argument("--out_tsv", required=True, type=Path, help="Filtered output TSV")
    ap.add_argument("--debug_max_headers", type=int, default=None, help="Only for --fasta_gz debug")
    args = ap.parse_args()

    if args.contig_list:
        contigs_ok = load_contigs_from_list(args.contig_list)
        source_desc = f"contig list: {args.contig_list}"
    else:
        contigs_ok = load_contigs_from_fasta_gz(args.fasta_gz, max_headers=args.debug_max_headers)
        source_desc = f"FASTA gz: {args.fasta_gz}"

    rows_in, rows_out = filter_pairs(args.pairs_tsv, contigs_ok, args.out_tsv)

    missing = rows_in - rows_out
    print(f"[scriptD] Loaded contigs from {source_desc}")
    print(f"[scriptD] Contigs available: {len(contigs_ok):,}")
    print(f"[scriptD] Pairs in (excluding header): {rows_in:,}")
    print(f"[scriptD] Pairs kept: {rows_out:,}")
    print(f"[scriptD] Pairs removed (missing contig): {missing:,}")
    print(f"[scriptD] Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()
