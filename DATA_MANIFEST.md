# Data manifest

Large raw datasets and bulky intermediate files are not tracked in this repository.

## Required external inputs

### BAPS-derived evaluation data
Required for:
- external GenoPHI evaluation
- construction of evaluation pairs
- comparison against labelled interactions

Recommended location:
- `data_raw/baps/`

### Original GenoPHI input data
Required for:
- feature construction
- retraining workflow
- mapping to the original protein-family space

Recommended location:
- `data_raw/genophi_original/`

### Reduced pilot subset inputs
Required for:
- reduced pilot retraining
- biologically separated feature construction

Recommended location:
- `data_raw/reduced_pilot/`

### Collaborator-provided E. coli phage genomes
Required for:
- Aim 2 TXTL candidate prioritisation

Recommended location:
- `data_raw/aim2_ecoli96/`

## Not tracked in Git

The following are intentionally excluded:
- raw FASTA / FNA / FAA files
- zip archives
- MMseqs2 databases and temporary files
- prodigal output directories
- large model artefacts
- large feature-construction folders
