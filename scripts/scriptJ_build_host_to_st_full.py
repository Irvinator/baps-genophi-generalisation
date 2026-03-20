#!/usr/bin/env python3
import re
from pathlib import Path

inp = Path("/mnt/d/baps-genophi-generalisation/mlst_results_bulk.tsv")
out = Path("/mnt/d/baps-genophi-generalisation/host_to_st_full.tsv")

seen = set()

with inp.open() as f, out.open("w") as g:
    g.write("host_accession\tST\n")

    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue

        path, _, st = parts[0], parts[1], parts[2]

        m = re.search(r"(GCA_\d+\.\d+)", path)
        if not m:
            continue

        host = m.group(1)

        st = st.replace("ecoli", "").strip()
        if st in {"", "-", "NA"}:
            st = "UNKNOWN"

        if host not in seen:
            seen.add(host)
            g.write(f"{host}\t{st}\n")

print("hosts_typed =", len(seen))
print("wrote =", out)
