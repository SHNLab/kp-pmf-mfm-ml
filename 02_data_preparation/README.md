# 02 — Data Preparation

**Pipeline:** Black-ginger PMF × MPMA Cascade — Master Feature Matrix for ML  
**Script:** `AI_PREP_FILE.R`  
**Language:** R (≥ 4.3.0)

---

## What this module does

Takes raw docking (MM-GBSA) and compound property (SwissADME) data and
produces a clean, fully-documented feature matrix ready for the ML pipeline
in `03_ml_pipeline/`.

Key operations:

- Recomputes methoxy group counts directly from IUPAC names — corrects stored
  `Meth_Count` values where they diverged from the IUPAC-derived count
- Classifies all 37 compounds into chemical classes (PMF, Meth_Flav,
  Standard_Flav, Gingeroid) using a rule-based engine
- Flags drug-likeness under strict and lenient criteria (Lipinski, Veber,
  Egan, PAINS, Brenk, TPSA)
- Merges MM-GBSA binding energies with all compound-level features
- Assigns compounds to lead tiers via Pareto non-dominated sorting across
  three objectives: inhibit-pool MM-GBSA, activate-pool MM-GBSA, and
  number of strong binders (< −50 kcal/mol)
- Computes per-target Mann-Whitney statistics and Cliff's delta (PMF vs
  non-PMF binders)
- Generates a selectivity index and promiscuity classification per compound

---

## Inputs

| File | Description |
|---|---|
| `mmgbsa_final.xlsx` | MM-GBSA binding energies — 12 targets × ~38 compounds |
| `COMPOUNDS_MF.xlsx` | 37 compounds with SwissADME ADMET output + positional substitution data |

Place both files in the working directory before running.

---

## Chemical classification logic

| Class | Rule |
|---|---|
| Gingeroid | Compound code in `{p33, p34, p35, p36, p37}` |
| PMF | IUPAC-derived methoxy count ≥ 2 |
| Meth_Flav | IUPAC-derived methoxy count = 1 |
| Standard_Flav | IUPAC-derived methoxy count = 0 |

Classification is based entirely on IUPAC name parsing — not on stored fields —
to prevent propagation of upstream annotation errors.

---

## MPMA target groups

| Group | Targets | Pharmacological goal |
|---|---|---|
| Inhibit | BCL2, HSP90, VDAC1, DRP1, PINK1, HK2 | Suppress oncogenic activity |
| Activate | MFN2, CPT1A, SIRT3 | Restore tumour-suppressive function |
| Extended | KEAP1, SIRT1, MMP9 | Secondary panel (excluded from primary ML) |

---

## Pareto non-dominated sorting

Each compound is evaluated on three simultaneous objectives:

1. Inhibit-pool mean MM-GBSA (lower = better)
2. Activate-pool mean MM-GBSA (lower = better)
3. Number of strong binders across all targets (> −50 kcal/mol threshold; higher = better)

A compound is **Pareto-optimal** if no other compound is strictly better on
all three objectives simultaneously.

Lead tiers are assigned as:

| Tier | Criteria |
|---|---|
| Tier 1 Gold | Pareto-optimal + top 10% single-target binder + drug-like strict |
| Tier 2 Silver | Pareto-optimal + top 10% single-target binder |
| Tier 3 Bronze | Pareto-optimal only |
| Tier 4 Potent | Top 10% single-target binder only |
| Tier 5 Other | Neither |

---

## ML feature schema (key columns)

| Role | Columns | Notes |
|---|---|---|
| Primary outcome (Y) | `MMGBSA dG Bind` | Regression target |
| Secondary outcome | `XP_GScore` | Validation target |
| Structural features | `M3–M8, M3p–M5p` | Positional methoxy flags |
| Structural features | `OH3, OH5, OH7, OH3p, OH4p` | Positional hydroxyl flags |
| Physicochemical | `LogP, TPSA, MR, Fraction Csp3` | SwissADME-derived |
| Binary classifier Y | `Is_PMF` | PMF vs non-PMF classification |
| Cascade classifier Y | `Pareto_Optimal` | Multi-objective lead selection |
| Leakage — exclude | `IUPAC_Methoxy, Chem_Class` | Directly defines PMF label |
| Groups | `p_code, Target_Protein` | Used for GroupKFold CV |

Full schema is saved to `ML_FEATURE_SCHEMA.csv`.

---

## Drug-likeness filters

| Filter | Threshold |
|---|---|
| Lipinski | 0 violations |
| Veber | 0 violations |
| PAINS | 0 alerts |
| TPSA | ≤ 140 Å² |
| Strict (all four) | `Drug_Like_Strict = 1` |
| Lenient (Lipinski + TPSA only) | `Drug_Like_Lenient = 1` |

---

## Outputs

| File | Rows | Description |
|---|---|---|
| `FINAL_MASTER_AI_DATASET.csv` | ~451 | Full feature matrix (37 compounds × 12 targets + controls) |
| `ML_FEATURE_SCHEMA.csv` | — | Column role documentation for ML pipeline |
| `COMPOUND_CASCADE_RANKING.csv` | 37 | Pareto + cascade scores, lead tier, promiscuity |
| `DRUGLIKE_SUBSET.csv` | subset | Compounds passing strict drug-likeness filter |
| `TARGET_ANALYSIS_SUMMARY.csv` | 12 | Per-target Mann-Whitney p, Cliff's delta, best compound |
| `CLASSIFICATION_AUDIT.csv` | 37 | Stored vs IUPAC-derived methoxy correction log |
| `SELECTIVITY_LONG.csv` | ~451 | Per compound × target selectivity index |

---

## Required packages

**CRAN only**
```
readxl, dplyr, tidyr, stringr, purrr, tibble
```

All packages are standard tidyverse — no Bioconductor dependencies.

---

## How to run

```r
# Place mmgbsa_final.xlsx and COMPOUNDS_MF.xlsx in working directory
source("AI_PREP_FILE_v5.R")
```

All 7 output CSVs are written to the working directory.
The script ends with a prompt to run the ML pipeline on `FINAL_MASTER_AI_DATASET.csv`.

---

## Unit tests

The methoxy counter is validated at startup against five known structures:

| IUPAC name | Expected count |
|---|---|
| 5-hydroxy-7-methoxychromen-4-one | 1 |
| 5,7-dimethoxychromen-4-one | 2 |
| 2-(3,4-dimethoxyphenyl)-5,7-dimethoxychromen-4-one | 4 |
| 5,7-dimethoxy-2-(3,4,5-trimethoxyphenyl)chromen-4-one | 5 |
| NA | 0 |

Script halts immediately if any test fails — prevents silent misclassification.

---

## Notes on leakage prevention

`IUPAC_Methoxy` and `Chem_Class` are flagged as **EXCLUDE** in the ML schema.
Both are derived from IUPAC names and directly define the `Is_PMF` label.
Including them as features would constitute target leakage. The positional
methoxy flags (`M3–M5p`) are the structurally informative equivalents and
are the correct inputs for the ML model.
