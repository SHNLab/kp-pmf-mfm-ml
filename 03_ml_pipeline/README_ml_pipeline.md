# 03 — ML Pipeline

**Script:** `v6.py` (internal version: `run_ml_pipeline_v5`)  
**Language:** Python (≥ 3.10)  
**Target journal:** Computers in Biology and Medicine (CBM)  
**Runtime:** approximately 20–40 minutes on a 16GB machine

---

## What this module does

End-to-end ML validation for predicting MM-GBSA binding energies of
*Kaempferia parviflora* PMFs against the 9-target MPMA cascade in CRC.

The pipeline runs 8 analyses and generates 9 publication-ready figures
in a single execution. No intermediate steps required.

---

## Target panel

9-target MPMA cascade — final locked panel:

| Tier | Targets | Basis |
|---|---|---|
| Tier 1 (above Q3 + druggability ≥ 7.5) | BCL2, HSP90, PINK1, HK2, CPT1A, MFN2 | Evidence + druggability |
| Tier 2 (cascade-obligate) | DRP1, VDAC1, SIRT3 | Below Q3 but post-translationally regulated; retained for biological completeness |

Removed from earlier versions: MMP9 (not mitochondrial), SIRT1 (nuclear-predominant), KEAP1 (below Q3, not cascade-obligate).

---

## Feature sets

| Set | Columns | Count |
|---|---|---|
| Positional | `M3, M5, M6, M7, M8, M3p, M4p, M5p, OH3, OH5, OH7, OH3p, OH4p` | 13 |
| Structural | `Is_Sugar, Sugar_OH_Count, Alkyl_Chain_Len, Has_5OH_IntramolecularHB` | 4 |
| Physicochemical | `LogP, TPSA, #Heavy atoms, #Rotatable bonds, #H-bond acceptors, #H-bond donors, Fraction Csp3, MR, Synthetic Accessibility` | 9 |
| Full (default) | All three sets combined | 26 |

---

## Analysis stack (8 analyses)

### 1. Nested Cross-Validation
- Outer loop: Leave-One-Compound-Out (LOCO) — each compound left out once
- Inner loop: 5-fold grid search on remaining training data
- Prevents hyperparameter overfitting to test compounds
- Reports honest R² and MAE per target

Hyperparameter grid searched:

| Parameter | Values |
|---|---|
| `n_estimators` | 150, 300, 500 |
| `max_depth` | 2, 3, 4 |
| `learning_rate` | 0.03, 0.05, 0.10 |
| `subsample` | 0.7, 0.85 |

### 2. Conformal Prediction Intervals
- Distribution-free 95% coverage guarantee (no normality assumption)
- Train/calibration split: 70/30 within each LOCO fold
- Calibration residuals determine interval width per compound
- Reports empirical coverage vs 95% target per target

### 3. Y-Randomization
- 500 permutations per target
- Real R² compared against null distribution
- p-value = fraction of null R² ≥ real R²
- Significant if p < 0.05 — confirms model learns signal, not noise

### 4. Gaussian Process Regression (GPR)
- Kernel: `ConstantKernel × Matérn(ν=2.5) + WhiteKernel`
- LOCO-CV predictions with calibrated uncertainty (σ)
- 95% CI: prediction ± 1.96σ
- Reports coverage and flags high-uncertainty compounds

### 5. Applicability Domain (Williams Plot)
- KNN distance in standardised feature space (k=5)
- Threshold: mean + 3 × SD of KNN distances
- Compounds beyond threshold flagged as out-of-domain

### 6. Baseline Benchmarks
Three reduced models run per target to isolate signal source:
- Physchem-only (no structural features)
- Position-only (methoxy/hydroxyl flags only)
- Full features (default)

Delta metrics quantify the contribution of positional vs physicochemical features — the key biological interpretation (methoxylation pattern governs binding).

### 7. Learning Curve
- Training fractions: 30%, 40%, 50%, 60%, 70%, 80%, 90%
- 10 replicates per fraction (random splits)
- Reports R² mean ± SD vs training size

### 8. Active Learning Simulation
- GPR uncertainty sampling vs random acquisition
- 20 simulation replicates
- Reported on CPT1A and VDAC1 (PMF-favoured targets)
- Shows efficiency gain of uncertainty-guided compound selection

---

## Inputs

| File | Source | Notes |
|---|---|---|
| `FINAL_MASTER_AI_DATASET.csv` | `02_data_preparation/` | Primary input |
| `FINAL_MASTER_AI_DATASET_corrected.csv` | Optional override | Preferred if positional flags were re-derived; script auto-detects |
| `COMPOUND_CASCADE_RANKING.csv` | `02_data_preparation/` | Chemotype mapping |
| `CLASSIFICATION_AUDIT.csv` | `02_data_preparation/` | Fallback for chemotype assignment |

---

## Outputs

**CSV / XLSX:**

| File | Description |
|---|---|
| `NESTED_CV_RESULTS.csv` | Per-target nested R² and MAE |
| `HYPERPARAMETER_SCAN.csv` | Modal best hyperparameters per target |
| `CONFORMAL_PREDICTIONS.csv` | Per-compound prediction intervals + coverage flag |
| `Y_RANDOMIZATION_RESULTS.csv` | Real R², null distribution, p-value per target |
| `GPR_UNCERTAINTY_PREDICTIONS.csv` | GPR predictions, σ, CI bounds, within-CI flag |
| `APPLICABILITY_DOMAIN.csv` | KNN distances, threshold, in/out-domain flag |
| `BASELINE_BENCHMARKS.csv` | R² for physchem-only, position-only, full features |
| `LEARNING_CURVE.csv` | R² mean ± SD at each training fraction |
| `ACTIVE_LEARNING_SIMULATION.csv` | Uncertainty vs random acquisition curves |
| `ML_VALIDATION_COMPLETE.xlsx` | All above tables in one workbook |

**Figures** (saved to `figures_v5/` as `.png` + `.svg` at 600 DPI):

| File | Description |
|---|---|
| `FigV5A_y_randomization` | Real R² vs null distribution — 9 targets |
| `FigV5B_applicability_domain` | Williams plot — KNN distance vs leverage |
| `FigV5C_gpr_uncertainty` | GPR predictions with 95% CI by chemotype |
| `FigV5D_learning_curve` | R² vs training fraction |
| `FigV5E_active_learning` | Uncertainty vs random acquisition — CPT1A + VDAC1 |
| `FigV5F_baseline_benchmarks` | Physchem vs position vs full feature R² |
| `FigV5G_conformal_coverage` | Empirical vs target coverage + interval width |
| `FigV5H_nested_cv_comparison` | Nested CV vs single-fold XGBoost R² |
| `FigV5_hero_validation` | 4-panel summary (Y-rand, nested CV, baselines, AL) |

All figures use Times New Roman 12pt, 600 DPI — publication-ready for CBM submission.

---

## Required packages

```
numpy
pandas
scikit-learn
xgboost
matplotlib
seaborn
scipy
openpyxl      # for XLSX output
```

Install:
```bash
pip install numpy pandas scikit-learn xgboost matplotlib seaborn scipy openpyxl
```

---

## How to run

```bash
# Place all input CSVs in the working directory
python v6.py
```

All outputs are written to the working directory.  
Figures are written to `figures_v5/`.

---

## Key results (locked from validated run)

| Metric | Value |
|---|---|
| XGBoost LOCO-CV (best target) | r = 0.776 |
| GPR LOCO-CV (best target) | r = 0.831 |
| Y-randomization ROC-AUC | 0.97 ± 0.02 |
| Conformal mean coverage | ≥ 95% (target met) |

---

## Design decisions

- **LOCO-CV over random split:** each compound is the test set once, preventing data leakage between structurally similar compounds docked to the same target
- **Nested CV:** hyperparameters tuned only on training data per fold — no test-set leakage from grid search
- **Conformal intervals:** distribution-free coverage guarantee; no assumption of Gaussian residuals — appropriate for a small, chemically diverse dataset
- **Baseline benchmarks:** position-only vs physchem-only comparison directly supports the methoxylation-position-governs-binding hypothesis in the manuscript
- **`tree_method='hist'`:** fast histogram-based XGBoost — required for reasonable runtime on CPU; equivalent results to exact method on this dataset size
