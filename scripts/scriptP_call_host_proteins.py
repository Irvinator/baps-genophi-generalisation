#!/usr/bin/env python3
import csv
import subprocess
from pathlib import Path
from multiprocessing import Pool

map_tsv = Path("/mnt/d/baps-genophi-generalisation/retrain_baseline/inputs/all_host_fna_map.tsv")
out_dir = Path("/mnt/d/baps-genophi-generalisation/host_AAs_baps")
out_dir.mkdir(parents=True, exist_ok=True)

NPROC = 10  # leave a couple of cores free

jobs = []
with map_tsv.open(newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        acc = row["host_accession"]
        fna = row["genome_fna"]
        faa = out_dir / f"{acc}.faa"
        jobs.append((acc, fna, str(faa)))

def run_one(job):
    acc, fna, faa = job
    out = Path(faa)
    if out.exists() and out.stat().st_size > 0:
        return ("skipped", acc)

    cmd = [
        "prodigal",
        "-i", fna,
        "-a", faa,
        "-p", "single",
        "-q"
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return ("done", acc)
    except subprocess.CalledProcessError:
        return ("failed", acc)

if __name__ == "__main__":
    done = skipped = failed = 0
    with Pool(NPROC) as pool:
        for i, (status, acc) in enumerate(pool.imap_unordered(run_one, jobs), start=1):
            if status == "done":
                done += 1
            elif status == "skipped":
                skipped += 1
            else:
                failed += 1

            if i % 200 == 0:
                print(f"processed={i} done={done} skipped={skipped} failed={failed}")

    print(f"finished done={done} skipped={skipped} failed={failed}")
