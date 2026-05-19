# 01 — Target Selection

**Pipeline:** Multi-omic Target Prioritisation for the Mitochondrial Fragmentation-to-Metastasis (MFM) Axis in Colorectal Cancer  
**Script:** `target_selection_v4.R`  
**Language:** R (≥ 4.3.0)

---

## What this module does

Identifies and ranks drug-targetable genes within the MFM axis using evidence integrated across three independent transcriptomic layers:

| Layer | Dataset | Design |
|---|---|---|
| Mechanistic | GSE80320 | mtDNA-deficient vs wild-type CRC cells; 2×2 factorial (oxygen × mtDNA status) |
| Functional | GSE32323 | Paired tumour/normal CRC tissue; patient-blocked limma contrast |
| Clinical | TCGA-COAD + TCGA-READ | Primary tumour vs solid tissue normal; DESeq2 full-count pipeline |

Per-layer log₂ fold-changes are Z-scored, summed into a composite score, and ranked genome-wide. The final panel is filtered at Q3 (≥75th percentile composite score) and druggability ≥ 7.5/10.

---

## Target panel — MPMA (9 targets across 5 nodes)

| Gene | Node | CRC Role | Pharmacological Mode |
|---|---|---|---|
| HSP90AB1 | Node 1: Proteostasis | Upregulated | INHIBIT — ATP pocket |
| BCL2 | Node 2: Apoptotic gate | Upregulated | INHIBIT — BH3 groove |
| VDAC1 | Node 2: Apoptotic gate | Upregulated | INHIBIT — HK2 interface |
| MFN2 | Node 3: Dynamics/Mitophagy | Downregulated | ACTIVATE — GTPase pocket |
| PINK1 | Node 3: Dynamics/Mitophagy | Upregulated | INHIBIT — kinase ATP pocket |
| DRP1 | Node 3: Dynamics/Mitophagy | Hyperactive | INHIBIT — GTPase pocket |
| CPT1A | Node 4: Lipid coupling (FAO) | Context-dependent | MODULATE — CoA-binding pocket |
| HK2 | Node 4: Lipid coupling (FAO) | Upregulated | INHIBIT — catalytic glucose pocket |
| SIRT3 | Node 5: Deacetylase axis | Downregulated | ACTIVATE — allosteric NAD+ site |

---

## Pipeline phases

```
Phase 0  — Package management + namespace conflict resolution
Phase 1  — Target panel definition + gene alias engine
Phase 2  — Mechanistic layer (GSE80320) — limma DE
Phase 2b — Path A vs Path B sensitivity (2×2 factorial validation)
Phase 3  — Functional layer (GSE32323) — paired limma DE
Phase 4  — Clinical layer (TCGA-COAD + TCGA-READ) — DESeq2
Phase 5  — Multi-omic integration + imputation sensitivity analysis
Phase 5b — Legacy vs current config diagnostic
Phase 6  — GSEA validation (2-layer and 3-layer)
Phase 7  — Druggability filter
Phase 8  — Visualisations (5 publication-ready figures)
Phase 9  — Session info
```

---

## Imputation sensitivity

Missing layer data (gene absent from platform) is handled under four strategies and compared:

| Strategy | Description |
|---|---|
| zero | Missing → 0 (original baseline) |
| mean | Missing → layer mean \|logFC\| (adopted primary) |
| complete | Inner join — genes present in all 3 layers only |
| min2of3 | Genes present in ≥2 layers; mean-impute the missing layer |

A target is **imputation-stable** if it passes Q3 in ≥3 of 4 strategies. All 9 MPMA targets satisfy this criterion in the validated pipeline.

---

## Required packages

**CRAN**
```
dplyr, ggplot2, ggrepel, tidyr, stringr, scales, patchwork, tibble, readr, conflicted
```

**Bioconductor**
```
GEOquery, limma, TCGAbiolinks, DESeq2, fgsea, msigdbr, SummarizedExperiment
```

Packages are auto-installed on first run via the `install_if_missing()` helper at the top of the script.

---

## Inputs

| File | Source | Notes |
|---|---|---|
| GSE80320 | GEO (auto-fetched) | mtDNA-deficient CRC cells |
| GSE32323 | GEO (auto-fetched) | Paired CRC tumour/normal |
| TCGA-COAD | GDC API (auto-fetched) | Requires stable internet; fallback to `coad.csv` |
| TCGA-READ | GDC API (auto-fetched) | Requires stable internet; fallback to `read.csv` |

If network access fails for TCGA, place `coad.csv` and `read.csv` (GEPIA-format DE tables) in the working directory. The script detects and uses these automatically.

---

## Outputs

All outputs are written to `output/` and `output/figures/`:

| File | Description |
|---|---|
| `target_summary_MPMA.csv` | Per-target integrated scores, ranks, percentiles |
| `final_MPMA_panel.csv` | Q3 + druggability filter flags |
| `master_ranking_full.csv` | Full genome-wide composite ranking |
| `coverage_audit.csv` | Per-target layer presence audit |
| `imputation_sensitivity_MPMA.csv` | Q3-pass stability across 4 strategies |
| `pathA_vs_pathB_logFC.csv` | GSE80320 design sensitivity comparison |
| `legacy_vs_current_diagnostic.csv` | Tier change audit across pipeline versions |
| `TCGA_COAD_DE_full.csv` | Full TCGA-COAD DESeq2 results |
| `TCGA_READ_DE_full.csv` | Full TCGA-READ DESeq2 results |
| `gsea_3layer.csv` | GSEA results — 3-layer ranking |
| `gsea_2layer.csv` | GSEA results — 2-layer ranking (sparsity-robust) |
| `session_info.txt` | R session info for reproducibility |
| `figures/Fig_Galaxy_v4.png` | Rank vs composite Z-score |
| `figures/Fig_Druggability_v4.png` | Evidence × druggability dual-filter scatter |
| `figures/Fig_GSEA_v4.png` | 2-layer vs 3-layer GSEA enrichment |
| `figures/Fig_LayerPercentile_v4.png` | Per-layer percentile by target |
| `figures/Fig_MPMA_Cascade.png` | 5-node MPMA cascade schematic |
| `figures/Fig_Landscape_MPMA.png` | Full evidence × druggability landscape |
| `figures/Fig2_MultiOmics_Heatmap.png` | Multi-omics heatmap — per-layer logFC |

---

## How to run

```r
# From the repo root or the 01_target_selection/ directory:
source("target_selection_v4.R")
```

Set a stable working directory before running. All outputs are written relative to `./output/`.

---

## Key design decisions

- **Path B adopted** over Path A for GSE80320: the clean 2×2 factorial design (samples 5,6,7,8 only) separates mtDNA and oxygen effects without pooling biologically distinct rho0n and rho0x sub-variants. All 9 MPMA targets are sign-concordant across both paths.
- **Mean imputation adopted** over zero-imputation as primary strategy: prevents inflation of Z-scores for platform-absent genes. Zero-imputation baseline retained for sensitivity reporting.
- **TCGA full-count DESeq2** preferred over GEPIA: direct access to raw counts avoids pre-filtering ambiguities. GEPIA fallback retained for offline use.

---

## Citation

If you use this pipeline, please cite the associated manuscript:

> Sharma A. et al. *Polymethoxyflavones from Kaempferia parviflora modulate the Mitochondrial Fragmentation-to-Metastasis axis in colorectal cancer: a multi-omic computational study.* (Under review, 2026)
