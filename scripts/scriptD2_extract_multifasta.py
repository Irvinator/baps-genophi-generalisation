#!/usr/bin/env python3
import argparse
import gzip
import sys

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baps_fasta_gz", required=True, help="BAPS_lytic_phage_0.1.fasta.gz")
    ap.add_argument("--contig_list", required=True, help="One contig per line (e.g. NAFV01000136.1)")
    ap.add_argument("--out_fasta_gz", required=True, help="Output multi-fasta .fa.gz")
    return ap.parse_args()

def main():
    args = parse_args()

    wanted = set()
    with open(args.contig_list, "r") as f:
        for line in f:
            c = line.strip()
            if c:
                wanted.add(c)

    print(f"[scriptD2] Wanted contigs: {len(wanted):,}", file=sys.stderr)

    kept = 0
    total = 0
    write = False

    with gzip.open(args.baps_fasta_gz, "rt") as fin, gzip.open(args.out_fasta_gz, "wt") as fout:
        for line in fin:
            if line.startswith(">"):
                total += 1
                # header format contains "__<contig>__"
                # example: >Escherichia_coli__GCA_002099625.1_-_ASM209962v1__NAFV01000136.1__259__562 ...
                hdr = line.strip()
                contig = None
                parts = hdr.split("__")
                if len(parts) >= 3:
                    contig = parts[2]  # the contig token
                write = (contig in wanted)
                if write:
                    kept += 1
                    fout.write(line)
                continue

            if write:
                fout.write(line)

    print(f"[scriptD2] FASTA records scanned: {total:,}", file=sys.stderr)
    print(f"[scriptD2] FASTA records kept:   {kept:,}", file=sys.stderr)
    print(f"[scriptD2] Wrote: {args.out_fasta_gz}", file=sys.stderr)

if __name__ == "__main__":
    main()
