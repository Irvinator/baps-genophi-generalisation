import pandas as pd
import numpy as np
from pathlib import Path

IN_POS = Path.home() / "data/BAPS/outputs/baps_ecoli_pos_pairs.tsv"
OUT_DIR = Path.home() / "data/BAPS/outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- Controls (edit these) ----
N_HOSTS = 200                  # number of bacterial hosts (accessions)
MAX_POS_PER_HOST = 10          # positives sampled per host (cap)
NEG_PER_POS = 1                # how many negatives per positive (1 => balanced)
SEED = 42
# --------------------------------

rng = np.random.default_rng(SEED)

# Load positives
pos = pd.read_csv(IN_POS, sep="\t")
pos = pos.dropna().drop_duplicates()

# Unique hosts and phages
hosts = pos["host_accession"].unique()
phages = pos["phage_contig"].unique()

print(f"Loaded positives: {len(pos):,}")
print(f"Unique hosts: {len(hosts):,}")
print(f"Unique phage contigs: {len(phages):,}")

# Sample hosts
if N_HOSTS > len(hosts):
    N_HOSTS = len(hosts)
sampled_hosts = rng.choice(hosts, size=N_HOSTS, replace=False)

# Build positives subset
pos_sub = (
    pos_filtered.groupby("host_accession", group_keys=False)
    .apply(lambda df: df.sample(n=min(len(df), MAX_POS_PER_HOST), random_state=SEED))
    .reset_index(drop=True)
))

pos_sub = pos_sub.rename(columns={"host_accession": "strain", "phage_contig": "phage"})
pos_sub["interaction"] = 1

# Build negatives
neg_rows = []
pos_by_host = (
    pos[pos["host_accession"].isin(sampled_hosts)]
    .groupby("host_accession")["phage_contig"]
    .apply(set)
    .to_dict()
)

for host in sampled_hosts:
    host_pos = pos_by_host.get(host, set())
    if not host_pos:
        continue

    # candidates are all phages not in host positives
    candidates = np.array(list(set(phages) - host_pos))
    if len(candidates) == 0:
        continue

    # number of negatives to sample
    n_pos = min(len(host_pos), MAX_POS_PER_HOST)
    n_neg = n_pos * NEG_PER_POS

    chosen = rng.choice(candidates, size=min(n_neg, len(candidates)), replace=False)
    for ph in chosen:
        neg_rows.append((host, ph, 0))

neg = pd.DataFrame(neg_rows, columns=["strain", "phage", "interaction"])

# Combine
df = pd.concat([pos_sub[["strain","phage","interaction"]], neg], ignore_index=True)

# Shuffle rows
df = df.sample(frac=1, random_state=SEED).reset_index(drop=True)

out_path = OUT_DIR / f"baps_ecoli_posneg_hosts{N_HOSTS}_pos{MAX_POS_PER_HOST}_neg{NEG_PER_POS}_seed{SEED}.csv"
df.to_csv(out_path, index=False)

print("\nOutput:", out_path)
print(df["interaction"].value_counts())
print("\nExample:")
print(df.head(10).to_string(index=False))
