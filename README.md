# kp-pmf-mfm-ml

Computational pipeline accompanying the manuscript:

**Sharma A., Sharma A., Rani A., Kaur S., Kuksal K., Nile S.H. (2026).** *Position-specific O-methylation drives selectivity of Kaempferia parviflora polymethoxyflavones against the colorectal-cancer mitochondrial survival network.*

> **Note:** This repository was originally named for the "Mitochondrial Fragmentation-to-Metastasis (MFM) axis." The cascade is referred to as the **Mitochondrial Survival Network (MSN)** in the published manuscript. The two names refer to the same nine-target panel.

---

## Overview

This repository contains all custom analysis code used in the study, organised into four modules corresponding to the main pipeline stages:

| Folder | Purpose |
|---|---|
| `01_target_selection/` | Three-layer transcriptomic integration (GSE80320, GSE32323, TCGA-COAD/READ), composite Z-score scoring, druggability filtering. R-based, using limma, DESeq2, TCGAbiolinks. |
| `02_data_preparation/` | LC–MS feature curation, compound library construction, descriptor calculation (RDKit), positional methoxylation annotation. Python. |
| `03_ml_pipeline/` | XGBoost regression with nested CV, Random Forest classification, Gaussian Process Regression, Lasso with bootstrap CIs, SHAP, conformal prediction, UMAP/PCA embeddings, Tanimoto similarity, applicability domain. Python. |
| `04_mds/` | Molecular dynamics analysis: replicate aggregation, RMSD calculation, PBC-aware minimum-distance, PLIP per-frame interaction fingerprinting, figure generation. Python and shell. |

---

## Requirements

Python ≥ 3.12 with:
- scikit-learn 1.4
- XGBoost 2.1
- SHAP 0.44
- UMAP-learn 0.5
- RDKit 2024.03
- MDAnalysis (latest)
- matplotlib, numpy, pandas

R ≥ 4.3 with:
- limma 3.54
- DESeq2 1.42
- TCGAbiolinks 2.30

External tools:
- Schrödinger Maestro / Glide XP / Prime MM-GBSA (commercial, version 2024)
- GROMACS 2024.5
- PLIP (Protein-Ligand Interaction Profiler)

---

## Data availability

- **Raw transcriptomic data:** GSE80320 and GSE32323 from Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/); TCGA-COAD and TCGA-READ from Genomic Data Commons (https://portal.gdc.cancer.gov/).
- **Processed datasets, model outputs, and analytical tables:** see Supplementary Tables S1–S11 of the manuscript.
- **Molecular dynamics trajectories** (3 × 100 ns × 6 systems): not included due to size; available from the corresponding author on reasonable request.

---

## Citation

If you use this code, please cite the published manuscript (DOI to be added upon acceptance) and this Zenodo archive (DOI to be added upon Zenodo deposit).

---

## License

MIT License — see [LICENSE](LICENSE) for full terms.

---

## Contact

Corresponding author: Dr. Shivraj Hariram Nile, National Agri-Food Biotechnology Institute (NABI), Mohali, India. Email: shivraj.nile@nabi.res.in
