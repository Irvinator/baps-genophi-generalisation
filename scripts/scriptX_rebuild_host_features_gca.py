#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from genophi.mmseqs2_clustering import run_feature_assignment


FAA_DIR = Path(
    "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/genophi_inputs/train_host_AAs"
)
BEST_HITS = Path(
    "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/train_feature_build/strain/best_hits.tsv"
)
OUT_DIR = Path(
    "/mnt/d/baps-genophi-generalisation/retrain_reduced_rows1500_ph150/train_feature_build/strain"
)

OUT_PA = OUT_DIR / "presence_absence_matrix_gca.csv"
OUT_PA_MIN2 = OUT_DIR / "presence_absence_matrix_gca_min2.csv"
OUT_FEATURE_DIR = OUT_DIR / "features_gca_min2"


def build_protein_to_genome_map(faa_dir: Path) -> tuple[dict[str, str], list[str]]:
    protein_to_genome: dict[str, str] = {}
    genomes: list[str] = []

    faa_files = sorted(faa_dir.glob("*.faa"))
    if not faa_files:
        raise FileNotFoundError(f"No .faa files found in {faa_dir}")

    for faa in faa_files:
        genome = faa.stem
        genomes.append(genome)

        with faa.open("r") as f:
            for line in f:
                if line.startswith(">"):
                    protein_id = line[1:].strip().split()[0]
                    protein_to_genome[protein_id] = genome

    return protein_to_genome, genomes


def build_genome_feature_sets(
    best_hits_path: Path, protein_to_genome: dict[str, str]
) -> tuple[dict[str, set[str]], list[str], int]:
    genome_features: dict[str, set[str]] = defaultdict(set)
    all_features: set[str] = set()
    missing = 0

    if not best_hits_path.exists():
        raise FileNotFoundError(f"best_hits.tsv not found: {best_hits_path}")

    with best_hits_path.open("r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            query_protein = parts[0]
            family_id = parts[1]

            genome = protein_to_genome.get(query_protein)
            if genome is None:
                missing += 1
                continue

            genome_features[genome].add(family_id)
            all_features.add(family_id)

    return genome_features, sorted(all_features), missing


def write_presence_absence_matrix(
    out_path: Path,
    genomes: list[str],
    features: list[str],
    genome_features: dict[str, set[str]],
) -> None:
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Genome"] + features)

        for genome in sorted(genomes):
            featset = genome_features.get(genome, set())
            row = [genome] + [1 if feat in featset else 0 for feat in features]
            writer.writerow(row)


def filter_features_by_prevalence(
    genome_features: dict[str, set[str]], min_prevalence: int
) -> list[str]:
    counts: dict[str, int] = defaultdict(int)

    for featset in genome_features.values():
        for feat in featset:
            counts[feat] += 1

    kept = sorted([feat for feat, count in counts.items() if count >= min_prevalence])
    return kept


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_FEATURE_DIR.mkdir(parents=True, exist_ok=True)

    print("Building protein -> GCA genome mapping...")
    protein_to_genome, genomes = build_protein_to_genome_map(FAA_DIR)
    print(f"Mapped proteins to genomes: {len(protein_to_genome):,}")
    print(f"Genomes found: {len(genomes):,}")

    print("Reading best hits and building genome feature sets...")
    genome_features, all_features, missing = build_genome_feature_sets(
        BEST_HITS, protein_to_genome
    )
    print(f"Unique raw features: {len(all_features):,}")
    print(f"Queries missing genome mapping: {missing:,}")

    print("Writing full GCA-based presence/absence matrix...")
    write_presence_absence_matrix(
        OUT_PA, genomes, all_features, genome_features
    )
    print(f"Wrote: {OUT_PA}")

    print("Applying prevalence filter: features present in >= 2 genomes...")
    kept_features = filter_features_by_prevalence(genome_features, min_prevalence=2)
    print(f"Kept features (>=2 genomes): {len(kept_features):,}")

    print("Writing filtered GCA-based presence/absence matrix...")
    write_presence_absence_matrix(
        OUT_PA_MIN2, genomes, kept_features, genome_features
    )
    print(f"Wrote: {OUT_PA_MIN2}")

    print("Running GenoPHI host feature assignment on filtered GCA matrix...")
    run_feature_assignment(
        input_file=str(OUT_PA_MIN2),
        output_dir=str(OUT_FEATURE_DIR),
        source="strain",
        select="none",
        select_column="strain",
        max_ram=12,
        threads=12,
    )
    print("Host feature assignment complete.")
    print(f"Feature outputs saved in: {OUT_FEATURE_DIR}")


if __name__ == "__main__":
    main()
