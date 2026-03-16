# GenoPHI Generalisation Analysis on BAPS Dataset

This repository contains the code and analysis used to evaluate the generalisation performance of the GenoPHI phage–host interaction prediction framework on an external dataset derived from BAPS (Bacterial Assembly-Associated Phage Sequences).

## Project Overview

Predicting phage–host interactions is important for applications such as phage therapy and microbiome engineering. GenoPHI is a machine learning framework that predicts strain-level phage–host interactions using genome-derived features.

This project evaluates whether GenoPHI generalises to an independent dataset of phage–host interactions derived from BAPS.

The key goals were:

- Apply trained GenoPHI models to an external dataset
- Evaluate predictive performance on unseen phage genomes
- Investigate causes of reduced performance
- Analyse protein family overlap between training and evaluation datasets

## Dataset

The evaluation dataset is derived from the **BAPS phage–host interaction dataset**.

A **randomly sampled subset** of the BAPS dataset was used to create an independent evaluation test set.

Dataset characteristics:

- 200 bacterial strains
- 144 phages
- 28,800 strain–phage prediction combinations
- 428 labelled interaction pairs

## Methods

The evaluation pipeline followed these steps:

1. Assign proteins from BAPS strains and phages to GenoPHI protein families using **MMseqs2**
2. Convert assignments to **binary presence/absence feature tables**
3. Generate predictions using trained **CatBoost models**
4. Aggregate prediction confidence scores across training runs
5. Evaluate predictions using standard metrics:

- ROC-AUC
- PR-AUC
- Accuracy
- Precision
- Recall
- MCC

## Results Summary

Model performance on the BAPS dataset:

| Metric | Value |
|------|------|
| ROC-AUC | 0.44 |
| Accuracy | 0.16 |
| Precision | 0 |
| Recall | 0 |

The model predicted **all interactions as negative** at the default threshold.

Protein family overlap analysis shows that many BAPS phages contain proteins not present in the GenoPHI training dataset, limiting feature representation and reducing prediction accuracy.

## Repository Structure
# GenoPHI Generalisation Analysis on BAPS Dataset

This repository contains the code and analysis used to evaluate the generalisation performance of the GenoPHI phage–host interaction prediction framework on an external dataset derived from BAPS (Bacterial Assembly-Associated Phage Sequences).

## Project Overview

Predicting phage–host interactions is important for applications such as phage therapy and microbiome engineering. GenoPHI is a machine learning framework that predicts strain-level phage–host interactions using genome-derived features.

This project evaluates whether GenoPHI generalises to an independent dataset of phage–host interactions derived from BAPS.

The key goals were:

- Apply trained GenoPHI models to an external dataset
- Evaluate predictive performance on unseen phage genomes
- Investigate causes of reduced performance
- Analyse protein family overlap between training and evaluation datasets

## Dataset

The evaluation dataset is derived from the **BAPS phage–host interaction dataset**.

A **randomly sampled subset** of the BAPS dataset was used to create an independent evaluation test set.

Dataset characteristics:

- 200 bacterial strains
- 144 phages
- 28,800 strain–phage prediction combinations
- 428 labelled interaction pairs

## Methods

The evaluation pipeline followed these steps:

1. Assign proteins from BAPS strains and phages to GenoPHI protein families using **MMseqs2**
2. Convert assignments to **binary presence/absence feature tables**
3. Generate predictions using trained **CatBoost models**
4. Aggregate prediction confidence scores across training runs
5. Evaluate predictions using standard metrics:

- ROC-AUC
- PR-AUC
- Accuracy
- Precision
- Recall
- MCC

## Results Summary

Model performance on the BAPS dataset:

| Metric | Value |
|------|------|
| ROC-AUC | 0.44 |
| Accuracy | 0.16 |
| Precision | 0 |
| Recall | 0 |

The model predicted **all interactions as negative** at the default threshold.

Protein family overlap analysis shows that many BAPS phages contain proteins not present in the GenoPHI training dataset, limiting feature representation and reducing prediction accuracy.

## Repository Structure
scripts/
Prediction and analysis scripts

figures/
Generated plots and visualisations

results/
Processed prediction outputs and evaluation tables

report/
Dissertation sections and write-ups

README.md
.gitignore


## Key Findings

- GenoPHI performs well on its training datasets but **fails to generalise to the BAPS dataset**
- Many BAPS phages lack proteins present in the GenoPHI training feature space
- Dataset distribution shift likely contributes to reduced predictive performance

## Author

Irvin Kshirsagar  
Biomedical Engineering
