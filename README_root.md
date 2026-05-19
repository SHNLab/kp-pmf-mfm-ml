# kp-pmf-mfm-ml

**Polymethoxyflavones from *Kaempferia parviflora* modulate the Mitochondrial Fragmentation-to-Metastasis (MFM) Axis in colorectal cancer: a multi-omic computational study**

> Aman Sharma · SHN Lab · NABI Mohali / Panjab University · 2026  
> Manuscript under review

---

## Overview

This repository contains the complete computational pipeline for a multi-omic
study identifying and validating polymethoxyflavone (PMF) leads from black
ginger (*Kaempferia parviflora*) against a 9-target mitochondrial cascade in
colorectal cancer (CRC).

The pipeline runs in four sequential modules:

```
01_target_selection/    Multi-omic target prioritisation (R)
02_data_preparation/    Feature matrix construction (R)
03_ml_pipeline/         ML validation — XGBoost + GPR (Python)
04_molecular_docking/   100 ns MD simulation + MM-GBSA (GROMACS, HPC)
```

---

## Repository structure

```
kp-pmf-mfm-ml/
│
├── README.md                        ← you are here
├── .gitignore
│
├── 01_target_selection/
│   ├── target_selection_v4.R        ← multi-omic target prioritisation
│   └── README.md
│
├── 02_data_preparation/
│   ├── AI_PREP_FILE_v5.R            ← feature matrix + Pareto cascade scoring
│   └── README.md
│
├── 03_ml_pipeline/
│   ├── v6.py                        ← XGBoost + GPR + 8-analysis validation stack
│   └── README.md
│
├── 04_molecular_docking/
│   ├── 01_setup.sh                  ← system preparation
│   ├── 02_equilibrate.sh            ← EM + NVT + NPT
│   ├── 03_production.sh             ← 100 ns production MD
│   ├── 04_analyze.sh                ← trajectory analysis
│   └── README.md
│
└── data/                            ← gitignored (see below)
```

---

## The MFM axis — 9 targets across 5 nodes

| Node | Target | CRC role | Pharmacological mode |
|---|---|---|---|
| Node 1: Proteostasis | HSP90AB1 | Upregulated | INHIBIT |
| Node 2: Apoptotic gate | BCL2 | Upregulated | INHIBIT |
| Node 2: Apoptotic gate | VDAC1 | Upregulated | INHIBIT |
| Node 3: Dynamics/Mitophagy | MFN2 | Downregulated | ACTIVATE |
| Node 3: Dynamics/Mitophagy | PINK1 | Upregulated | INHIBIT |
| Node 3: Dynamics/Mitophagy | DRP1 | Hyperactive | INHIBIT |
| Node 4: Lipid coupling | CPT1A | Context-dependent | MODULATE |
| Node 4: Lipid coupling | HK2 | Upregulated | INHIBIT |
| Node 5: Deacetylase axis | SIRT3 | Downregulated | ACTIVATE |

---

## Key results

| Analysis | Result |
|---|---|
| Compound library | 37 compounds (33 PMFs + 4 gingeroids) across 9 targets |
| Top lead | PMF-01 (5,7,3',4'-tetramethoxyflavone) |
| XGBoost LOCO-CV | r = 0.776 |
| GPR LOCO-CV | r = 0.831 |
| Y-randomization ROC-AUC | 0.97 ± 0.02 |
| Conformal coverage | >= 95% (target met across all 9 targets) |
| MD simulation | 100 ns x 6 systems (BCL2, HSP90, SIRT3 + respective controls) |

---

## How to reproduce

### Module 1 — Target selection (R)

```r
# Install packages on first run (auto-handled by script)
source("01_target_selection/target_selection_v4.R")
# Outputs written to output/ and output/figures/
```

Requires internet access for GEO and TCGA data fetch. GEPIA fallback CSVs
(`coad.csv`, `read.csv`) can be placed in the working directory if offline.

### Module 2 — Data preparation (R)

```r
# Place mmgbsa_final.xlsx and COMPOUNDS_MF.xlsx in working directory
source("02_data_preparation/AI_PREP_FILE_v5.R")
# Outputs: 7 CSVs including FINAL_MASTER_AI_DATASET.csv
```

### Module 3 — ML pipeline (Python)

```bash
# Place output CSVs from Module 2 in working directory
pip install numpy pandas scikit-learn xgboost matplotlib seaborn scipy openpyxl
python 03_ml_pipeline/v6.py
# Outputs: 10 CSVs + 9 figures in figures_v5/
# Runtime: 20-40 min on 16GB CPU machine
```

### Module 4 — Molecular dynamics (HPC only)

```bash
# Requires Param Smriti HPC access or equivalent SLURM cluster with GPU nodes
# See 04_molecular_docking/README.md for full setup and SLURM commands
cd 04_molecular_docking/
./01_setup.sh ~/inputs/bcl2_p12.pdb LIG bcl2_p12
# followed by 02, 03, 04 via sbatch
```

---

## Data availability

Raw docking and compound data (`mmgbsa_final.xlsx`, `COMPOUNDS_MF.xlsx`) are
not included in this repository as they contain unpublished experimental
results. They will be deposited as supplementary data upon manuscript
acceptance.

GEO datasets (GSE80320, GSE32323) and TCGA data are publicly available and
fetched automatically by `target_selection_v4.R`.

---

## Software versions

| Tool | Version |
|---|---|
| R | >= 4.3.0 |
| Python | >= 3.10 |
| GROMACS | 2024.5 (CUDA build) |
| XGBoost | any recent |
| scikit-learn | any recent |
| MDAnalysis | >= 2.10 |
| gmx_MMPBSA | 1.6.4 |

Full R session info is written to `output/session_info.txt` on each run.

---

## Citation

If you use this pipeline or any part of it, please cite:

> Sharma A. et al. *Polymethoxyflavones from Kaempferia parviflora modulate
> the Mitochondrial Fragmentation-to-Metastasis axis in colorectal cancer:
> a multi-omic computational study.* (Under review, 2026)

Individual tool citations are listed in each module README.

---

## License

MIT License — see `LICENSE` file.

---

## Contact

**Aman Sharma**  
NABI Mohali / Panjab University  
SHN Lab  
GitHub: [SHNLab](https://github.com/SHNLab)
