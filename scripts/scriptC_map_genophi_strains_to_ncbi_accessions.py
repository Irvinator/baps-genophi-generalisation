#!/usr/bin/env python3
import json
import subprocess
from pathlib import Path

INFILE  = Path("/home/irvin/data/BAPS/outputs/genophi_ecoli_strain_names.txt")
OUTFILE = Path("/home/irvin/data/BAPS/outputs/genophi_ecoli_strain_to_accession.tsv")

TAXON_ECOLI = "562"

def run_datasets_search(term: str, max_hits: int = 5):
    """
    Return up to max_hits genome assembly records from NCBI Datasets for a search term.
    """
    cmd = [
        "datasets", "summary", "genome", "taxon", TAXON_ECOLI,
        "--search", term,
        "--as-json-lines"
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        return [], p.stderr.strip()

    records = []
    for line in p.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            continue

    # Some versions return one big JSON rather than json-lines; handle that too
    if len(records) == 0:
        try:
            obj = json.loads(p.stdout)
            # best effort: wrap in list
            records = [obj]
        except Exception:
            pass

    return records[:max_hits], ""

def pick_best_accession(records):
    """
    Best-effort selection:
    prefer reference/representative if present; otherwise take first accession found.
    """
    # Different Datasets versions nest these fields differently.
    # We'll scan for plausible assembly accession keys.
    def extract_accessions(obj):
        accs = set()

        # Common: obj["assemblies"][i]["assembly"]["accession"]
        if isinstance(obj, dict):
            assemblies = obj.get("assemblies") or obj.get("reports") or []
            if isinstance(assemblies, list):
                for a in assemblies:
                    if isinstance(a, dict):
                        # try several nestings
                        for path in [
                            ("assembly", "accession"),
                            ("accession",),
                            ("assembly_accession",),
                        ]:
                            cur = a
                            ok = True
                            for k in path:
                                if isinstance(cur, dict) and k in cur:
                                    cur = cur[k]
                                else:
                                    ok = False
                                    break
                            if ok and isinstance(cur, str) and (cur.startswith("GCA_") or cur.startswith("GCF_")):
                                accs.add(cur)

                        # sometimes accession is in top level of a
                        for k, v in a.items():
                            if isinstance(v, str) and (v.startswith("GCA_") or v.startswith("GCF_")):
                                accs.add(v)

            # also scan whole dict shallowly
            for k, v in obj.items():
                if isinstance(v, str) and (v.startswith("GCA_") or v.startswith("GCF_")):
                    accs.add(v)

        return sorted(accs)

    # try to prefer records containing "representative" / "reference"
    # (exact fields vary; so we just string-search JSON as a fallback)
    scored = []
    for r in records:
        accs = extract_accessions(r)
        blob = json.dumps(r).lower()
        score = 0
        if "reference" in blob:
            score += 2
        if "representative" in blob:
            score += 1
        scored.append((score, accs, r))

    scored.sort(reverse=True, key=lambda x: x[0])

    for score, accs, r in scored:
        if accs:
            return accs[0], score, accs

    return None, 0, []

def main():
    strains = [s.strip() for s in INFILE.read_text().splitlines() if s.strip()]
    OUTFILE.parent.mkdir(parents=True, exist_ok=True)

    with OUTFILE.open("w") as f:
        f.write("strain_name\tassembly_accession\tmatch_score\tall_accessions_found\n")

        for s in strains:
            term = s  # you can also do f'"{s}"' but datasets handles plain strings fine
            print(f"Searching: {s}", flush=True)

            recs, err = run_datasets_search(term)
            if err:
                f.write(f"{s}\t\t0\tERROR: {err}\n")
                continue

            best, score, all_accs = pick_best_accession(recs)
            if best is None:
                f.write(f"{s}\t\t0\tNO_MATCH\n")
            else:
                f.write(f"{s}\t{best}\t{score}\t{';'.join(all_accs)}\n")

    print(f"Wrote mapping: {OUTFILE}")

if __name__ == "__main__":
    main()
