# BAPS × GenoPHI Generalisation Experiments

## Project Overview

This project investigates the generalisation performance of GenoPHI on
novel Escherichia coli host–phage interactions derived from the BAPS
lytic phage dataset.

Objectives:

1.  Extract host–phage interactions from the BAPS dataset.
2.  Construct a balanced positive/negative dataset.
3.  Map GenoPHI training strain names to NCBI genome assemblies.
4.  Download host genomes from NCBI.
5.  Quantify genomic distance between GenoPHI training strains and novel
    BAPS strains.
6.  Evaluate generalisation using MASH distances.

------------------------------------------------------------------------

## Data Sources

BAPS Dataset: - current_BAPS_anns.tsv

GenoPHI Training Data: -
data/interation_matrices/ecoli_interaction_matrix_subset.csv

NCBI Genomes: Downloaded using: datasets download genome accession

------------------------------------------------------------------------

## Pipeline Overview

### 1. Extract BAPS Positive Host–Phage Pairs

Script: scripts/scriptA_build_baps_ecoli_pos_pairs.py

Input: current_BAPS_anns.tsv

Output: - outputs/baps_ecoli_pos_pairs.tsv -
outputs/baps_ecoli_accessions.txt

Purpose: Filter BAPS annotations for E. coli hosts and extract host
assembly accession (GCA\_\*) and associated phage contig. All extracted
pairs are labelled as positive (interaction = 1).

------------------------------------------------------------------------

### 2. Construct Balanced Positive/Negative Dataset

Script: scripts/scriptB_build_posneg_dataset.py

Input: outputs/baps_ecoli_pos_pairs.tsv

Output (example):
outputs/baps_ecoli_posneg_hosts200_pos10_neg1_seed42.csv

Method: - Randomly sample 200 E. coli host genomes - For each host: -
Sample up to 10 positive interactions - Generate 1 negative per positive
by pairing with a phage not annotated for that host - Fixed random seed
ensures reproducibility

Resulting dataset is balanced (e.g., 914 positives / 914 negatives).

------------------------------------------------------------------------

### 3. Map GenoPHI Strains to NCBI Assemblies

Script: scripts/scriptC_map_genophi_strains_to_ncbi_accessions.py

Input: data/interation_matrices/ecoli_interaction_matrix_subset.csv

Output: - outputs/genophi_ecoli_strain_to_accession.tsv -
outputs/genophi_ecoli_accessions.txt

Purpose: Map named GenoPHI training strains (e.g., ECOR48, BL21) to NCBI
assembly accessions (GCA/GCF). This enables downloading the correct
genome assemblies for distance analysis.

------------------------------------------------------------------------

## Genome Download (NCBI)

Genomes were downloaded in batches to avoid archive errors.

Example: datasets download genome accession –inputfile
genophi_batch_00.txt

Extract FASTA files: find genophi_batch\_\* -name “\*genomic.fna” \>
genophi_ecoli_157_genomes.txt

------------------------------------------------------------------------

## Genomic Distance Analysis (MASH)

Sketch genomes: mash sketch -o mash/genophi_ecoli_downloaded
genophi_ecoli_157_genomes.txt mash sketch -o mash/baps_ecoli_200
baps_ecoli_200_genomes.txt

Compute pairwise distances: mash dist mash/genophi_ecoli_downloaded.msh
mash/baps_ecoli_200.msh \> outputs/mash_dist_genophiTrain_vs_baps200.tsv

Extract minimum distance per BAPS host:
outputs/baps200_min_dist_to_genophi_train.csv

Example summary: - Mean Mash distance ≈ 0.008 - Approximate ANI ≈
0.992 - Some genomes nearly identical to GenoPHI training strains

------------------------------------------------------------------------

## Key Finding

The sampled BAPS E. coli genomes are highly genetically similar to the
GenoPHI training strains. This indicates that current evaluation
primarily tests near-distribution performance rather than far
generalisation.

------------------------------------------------------------------------

## Repository Structure

baps-genophi-generalisation/ ├── scripts/ ├── README.md └── .gitignore

Large genome downloads and FASTA files are excluded via .gitignore.

------------------------------------------------------------------------

## Reproducibility

Recommended environment:

conda create -n baps python=3.10 pandas mash

------------------------------------------------------------------------

## Future Work

-   Expand training set using larger BAPS interaction matrices
-   Retrain GenoPHI on more phylogenetically diverse host–phage pairs
-   Evaluate model performance versus genomic distance bins
