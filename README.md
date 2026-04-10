# baps-genophi-generalisation

BEng research project repository investigating the external generalisation of GenoPHI for strain-level phage–host prediction and exploratory prioritisation of candidate phages for TXTL rebootability testing.

## Project overview

This repository contains three main analytical strands:

1. External evaluation of the original GenoPHI framework on a BAPS-derived *E. coli* dataset
2. Reduced pilot retraining of GenoPHI on biologically separated host and phage partitions
3. Feature-based prioritisation of collaborator-provided *E. coli* phage genomes for future TXTL rebootability testing

## Repository structure

- `scripts/` — analysis and plotting scripts
- `results/` — processed result tables and split files
- `report_figures/` — figures used in the dissertation
- `aim2_ecoli96/` — compact derived outputs for Aim 2 prioritisation
- `retrain_baseline/` — smaller baseline inputs and metadata
- `references.bib` — bibliography for dissertation
- `environment.yml` — conda environment definition

## Environment setup

```bash
conda env create -f environment.yml
conda activate genophi
