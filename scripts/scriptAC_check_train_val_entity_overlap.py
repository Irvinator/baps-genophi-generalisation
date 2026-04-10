#!/usr/bin/env python3
from pathlib import Path

BASE = Path("/mnt/d/baps-genophi-generalisation/retrain_baseline/inputs")

def read_list(fp: Path):
    with open(fp, "r") as f:
        return {line.strip() for line in f if line.strip()}

train_hosts = read_list(BASE / "train_hosts.txt")
val_hosts   = read_list(BASE / "val_hosts.txt")
train_phages = read_list(BASE / "train_phages.txt")
val_phages   = read_list(BASE / "val_phages.txt")

shared_hosts = train_hosts & val_hosts
shared_phages = train_phages & val_phages

print("Train hosts:", len(train_hosts))
print("Val hosts:", len(val_hosts))
print("Shared hosts:", len(shared_hosts))
if shared_hosts:
    print("Example shared hosts:", list(sorted(shared_hosts))[:10])

print()

print("Train phages:", len(train_phages))
print("Val phages:", len(val_phages))
print("Shared phages:", len(shared_phages))
if shared_phages:
    print("Example shared phages:", list(sorted(shared_phages))[:10])
