# BAPS × GenoPHI generalisation experiments

This repo contains scripts used to:
1) extract E. coli positive phage-host pairs from BAPS metadata,
2) construct a balanced pos/neg dataset for training/testing,
3) map GenoPHI E. coli strain names to NCBI assembly accessions (GCA/GCF),
and support downstream MASH distance analysis.

## Scripts
- scripts/scriptA_build_baps_ecoli_pos_pairs.py
  - Input: current_BAPS_anns.tsv (BAPS annotations table)
  - Output: outputs/baps_ecoli_pos_pairs.tsv and outputs/baps_ecoli_accessions.txt
  - Purpose: extract E. coli host accessions and associated lytic phage contigs (positive pairs).

- scripts/scriptB_build_posneg_dataset.py
  - Input: outputs/baps_ecoli_pos_pairs.tsv
  - Output: outputs/baps_ecoli_posneg_hosts200_pos10_neg1_seed42.csv (example)
  - Purpose: sample hosts and build balanced positive + negative pairs.

- scripts/scriptC_map_genophi_strains_to_ncbi_accessions.py
  - Input: GenoPHI interaction matrix strain list (e.g. ecoli_interaction_matrix_subset.csv)
  - Output: outputs/genophi_ecoli_strain_to_accession.tsv and outputs/genophi_ecoli_accessions.txt
  - Purpose: map strain names (e.g. ECOR48) to assemblies in NCBI (GCA/GCF).

## Notes
- Large NCBI genome downloads and raw FASTA files are not committed (see .gitignore).
