# ============================================================================
# AI_PREP_FILE_v5.R
# ----------------------------------------------------------------------------
# Black-ginger PMF x MPMA cascade dataset — master prep for ML pipeline
#
# INPUTS:
#   - mmgbsa_final.xlsx      (12 targets x ~38 rows each)
#   - COMPOUNDS_MF.xlsx      (37 compounds + SwissADME ADMET output)
#
# OUTPUTS (7 CSVs):
#   - FINAL_MASTER_AI_DATASET.csv    (full feature matrix, ~451 rows)
#   - ML_FEATURE_SCHEMA.csv          (column role documentation)
#   - COMPOUND_CASCADE_RANKING.csv   (Pareto + cascade scores, 37 rows)
#   - DRUGLIKE_SUBSET.csv            (compounds passing bioavailability)
#   - TARGET_ANALYSIS_SUMMARY.csv    (per-target Mann-Whitney, effect size)
#   - CLASSIFICATION_AUDIT.csv       (Meth_Count correction log)
#   - SELECTIVITY_LONG.csv           (compound x target selectivity index)
#
# Designed for Q1 in silico manuscript (Molecules / IJMS tier)
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ============================================================================
# CONFIG
# ============================================================================
MMGBSA_FILE    <- "mmgbsa_final.xlsx"
COMPOUNDS_FILE <- "COMPOUNDS_MF.xlsx"
OUT_DIR        <- "."

GINGEROIDS          <- c("p33", "p34", "p35", "p36", "p37")
PMF_METH_THRESHOLD  <- 2

INHIBIT_TARGETS  <- c("BCL2", "HSP90", "VDAC1", "DRP1", "PINK1", "HK2")
ACTIVATE_TARGETS <- c("MFN2", "CPT1A", "SIRT3")
EXTENDED_TARGETS <- c("KEAP1", "SIRT1", "MMP9")

STRONG_BINDER_MMGBSA <- -50.0
TPSA_MAX             <- 140.0

# ============================================================================
# 1. CORRECTED METHOXY COUNTER from IUPAC name
# ============================================================================
count_methoxy <- function(name) {
  if (is.na(name) || is.null(name)) return(0L)
  s <- tolower(as.character(name))
  
  prefix_map <- list(
    hexamethoxy  = 6, pentamethoxy = 5, tetramethoxy = 4,
    trimethoxy   = 3, dimethoxy    = 2, monomethoxy  = 1
  )
  
  total <- 0L
  for (prefix in names(prefix_map)) {
    if (grepl(prefix, s, fixed = TRUE)) {
      matches <- length(gregexpr(prefix, s, fixed = TRUE)[[1]])
      n <- prefix_map[[prefix]]
      total <- total + n * matches
      s <- gsub(prefix, strrep("|", n), s, fixed = TRUE)
    }
  }
  if (grepl("methoxy", s, fixed = TRUE)) {
    bare <- length(gregexpr("methoxy", s, fixed = TRUE)[[1]])
    total <- total + bare
  }
  as.integer(total)
}

stopifnot(count_methoxy("5-hydroxy-7-methoxychromen-4-one") == 1)
stopifnot(count_methoxy("5,7-dimethoxychromen-4-one") == 2)
stopifnot(count_methoxy("2-(3,4-dimethoxyphenyl)-5,7-dimethoxychromen-4-one") == 4)
stopifnot(count_methoxy("5,7-dimethoxy-2-(3,4,5-trimethoxyphenyl)chromen-4-one") == 5)
stopifnot(count_methoxy(NA) == 0)
cat("count_methoxy unit tests passed\n")

count_hydroxy_iupac <- function(name) {
  if (is.na(name) || is.null(name)) return(0L)
  s <- tolower(as.character(name))
  prefix_map <- list(
    hexahydroxy  = 6, pentahydroxy = 5, tetrahydroxy = 4,
    trihydroxy   = 3, dihydroxy    = 2
  )
  total <- 0L
  for (prefix in names(prefix_map)) {
    if (grepl(prefix, s, fixed = TRUE)) {
      matches <- length(gregexpr(prefix, s, fixed = TRUE)[[1]])
      n <- prefix_map[[prefix]]
      total <- total + n * matches
      s <- gsub(prefix, strrep("|", n), s, fixed = TRUE)
    }
  }
  if (grepl("hydroxy", s, fixed = TRUE)) {
    bare <- length(gregexpr("hydroxy", s, fixed = TRUE)[[1]])
    total <- total + bare
  }
  as.integer(total)
}

# ============================================================================
# 2. LOAD DATA
# ============================================================================
cat("\n=== LOADING FILES ===\n")
mmgbsa    <- read_xlsx(MMGBSA_FILE)
compounds <- read_xlsx(COMPOUNDS_FILE)
cat(sprintf("mmgbsa_final.xlsx: %d rows x %d cols\n", nrow(mmgbsa), ncol(mmgbsa)))
cat(sprintf("COMPOUNDS_MF.xlsx: %d rows x %d cols\n", nrow(compounds), ncol(compounds)))

if ("glide gscore" %in% colnames(mmgbsa)) {
  mmgbsa$XP_GScore <- mmgbsa$`glide gscore`
} else if ("XP GScore" %in% colnames(mmgbsa)) {
  mmgbsa$XP_GScore <- mmgbsa$`XP GScore`
} else {
  mmgbsa$XP_GScore <- mmgbsa$`docking score`
}

# ============================================================================
# 3. RECOMPUTE METHOXY + CLASSIFY
# ============================================================================
cat("\n=== AUDITING STORED Meth_Count vs IUPAC-DERIVED ===\n")
compounds <- compounds %>%
  mutate(
    p_code            = annotation,
    IUPAC_Methoxy     = map_int(Official_Name, count_methoxy),
    IUPAC_Hydroxy     = map_int(Official_Name, count_hydroxy_iupac),
    Stored_Methoxy    = Meth_Count,
    Methoxy_Corrected = IUPAC_Methoxy != Stored_Methoxy
  )

n_corrected <- sum(compounds$Methoxy_Corrected, na.rm = TRUE)
cat(sprintf("Meth_Count corrected: %d/%d compounds\n", n_corrected, nrow(compounds)))

compounds <- compounds %>%
  mutate(
    Chem_Class = case_when(
      p_code %in% GINGEROIDS              ~ "Gingeroid",
      IUPAC_Methoxy >= PMF_METH_THRESHOLD ~ "PMF",
      IUPAC_Methoxy %in% c(1L, 2L)        ~ "Meth_Flav",
      IUPAC_Methoxy == 0L                 ~ "Standard_Flav",
      TRUE                                 ~ "Unknown"
    ),
    Is_PMF       = as.integer(Chem_Class == "PMF"),
    Is_Gingeroid = as.integer(Chem_Class == "Gingeroid"),
    Is_Flavonoid = as.integer(Chem_Class %in% c("PMF", "Meth_Flav", "Standard_Flav"))
  )

cat("\nFinal classification:\n")
print(table(compounds$Chem_Class))

# ============================================================================
# 4. POSITIONAL VALIDATION + 5-OH FLAG
# ============================================================================
positional_meth_cols <- c("M3", "M5", "M6", "M7", "M8", "M3p", "M4p", "M5p")
positional_oh_cols   <- c("OH3", "OH5", "OH7", "OH3p", "OH4p")

compounds <- compounds %>%
  rowwise() %>%
  mutate(
    Positional_Meth_Sum = sum(c_across(all_of(positional_meth_cols)), na.rm = TRUE),
    Positional_OH_Sum   = sum(c_across(all_of(positional_oh_cols)),   na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Meth_Position_Valid = (Positional_Meth_Sum == IUPAC_Methoxy) | (Chem_Class == "Gingeroid"),
    OH_Position_Valid   = abs(Positional_OH_Sum - Core_OH_Count) <= 1,
    Has_5OH_IntramolecularHB = as.integer(OH5 == 1 & !is.na(OH5))
  )

n_pos_mismatch <- sum(!compounds$Meth_Position_Valid, na.rm = TRUE)
cat(sprintf("\nPositional methoxy mismatches: %d\n", n_pos_mismatch))

# ============================================================================
# 5. DRUG-LIKENESS FLAGS (SwissADME)
# ============================================================================
get_col <- function(df, pattern) {
  m <- grep(pattern, colnames(df), value = TRUE)
  if (length(m) == 0) return(rep(NA_real_, nrow(df)))
  df[[m[1]]]
}

compounds$Lipinski_Violations <- get_col(compounds, "^Lipinski")
compounds$Veber_Violations    <- get_col(compounds, "^Veber")
compounds$Egan_Violations     <- get_col(compounds, "^Egan")
compounds$Ghose_Violations    <- get_col(compounds, "^Ghose")
compounds$Muegge_Violations   <- get_col(compounds, "^Muegge")
compounds$PAINS_Alerts        <- get_col(compounds, "^PAINS")
compounds$Brenk_Alerts        <- get_col(compounds, "^Brenk")

compounds <- compounds %>%
  mutate(
    Pass_Lipinski = as.integer(Lipinski_Violations == 0),
    Pass_Veber    = as.integer(Veber_Violations == 0),
    Pass_Egan     = as.integer(Egan_Violations == 0),
    Pass_PAINS    = as.integer(PAINS_Alerts == 0),
    Pass_Brenk    = as.integer(Brenk_Alerts == 0),
    Pass_TPSA     = as.integer(TPSA <= TPSA_MAX),
    Drug_Like_Strict  = as.integer(Pass_Lipinski & Pass_Veber & Pass_PAINS & Pass_TPSA),
    Drug_Like_Lenient = as.integer(Pass_Lipinski & Pass_TPSA),
    Total_RuleOf5_Violations = Lipinski_Violations + Veber_Violations +
      Egan_Violations + PAINS_Alerts + Brenk_Alerts
  )

cat(sprintf("\nDrug-like strict:  %d / %d\n",
            sum(compounds$Drug_Like_Strict),  nrow(compounds)))
cat(sprintf("Drug-like lenient: %d / %d\n",
            sum(compounds$Drug_Like_Lenient), nrow(compounds)))

# ============================================================================
# 6. MERGE
# ============================================================================
cat("\n=== MERGING MMGBSA + COMPOUND FEATURES ===\n")

mmgbsa <- mmgbsa %>%
  mutate(
    p_code = ifelse(
      grepl("_pc$", `Entry Name`, ignore.case = TRUE),
      NA_character_,
      str_extract(`Entry Name`, "^p\\d+")
    ),
    Is_Control = as.integer(grepl("_pc", `Entry Name`, ignore.case = TRUE))
  )

compound_features <- compounds %>%
  select(p_code, Compound, Official_Name, Pub_CID, smiles,
         LogP, TPSA, Total_OH_HBD, Core_OH_Count, Sugar_OH_Count,
         Is_Sugar, Alkyl_Chain_Len,
         M3, M5, M6, M7, M8, M3p, M4p, M5p,
         OH3, OH5, OH7, OH3p, OH4p,
         `#Heavy atoms`, `#Aromatic heavy atoms`, `Fraction Csp3`,
         `#Rotatable bonds`, `#H-bond acceptors`, `#H-bond donors`,
         MR, `Consensus Log P`, `Bioavailability Score`,
         `Synthetic Accessibility`,
         `GI absorption`, `BBB permeant`, `Pgp substrate`,
         `CYP1A2 inhibitor`, `CYP2C19 inhibitor`, `CYP2C9 inhibitor`,
         `CYP2D6 inhibitor`, `CYP3A4 inhibitor`,
         IUPAC_Methoxy, IUPAC_Hydroxy, Stored_Methoxy, Methoxy_Corrected,
         Chem_Class, Is_PMF, Is_Gingeroid, Is_Flavonoid,
         Has_5OH_IntramolecularHB,
         Lipinski_Violations, Veber_Violations, PAINS_Alerts, Brenk_Alerts,
         Drug_Like_Strict, Drug_Like_Lenient, Total_RuleOf5_Violations)

master <- mmgbsa %>%
  left_join(compound_features, by = "p_code") %>%
  mutate(
    Chem_Class    = ifelse(Is_Control == 1, "Positive_Control", Chem_Class),
    Is_PMF        = ifelse(Is_Control == 1, 0L, Is_PMF),
    Is_Gingeroid  = ifelse(Is_Control == 1, 0L, Is_Gingeroid),
    Is_Flavonoid  = ifelse(Is_Control == 1, 0L, Is_Flavonoid),
    IUPAC_Methoxy = ifelse(Is_Control == 1, 0L, IUPAC_Methoxy)
  )

cat(sprintf("Master dataset: %d rows, %d cols\n", nrow(master), ncol(master)))

# ============================================================================
# 7. CASCADE SCORING
# ============================================================================
cat("\n=== CASCADE SCORING ===\n")

cascade <- master %>%
  filter(Is_Control == 0) %>%
  group_by(p_code, Compound, Official_Name, Chem_Class,
           IUPAC_Methoxy, Is_PMF, Drug_Like_Strict, Drug_Like_Lenient,
           TPSA, LogP, `Bioavailability Score`) %>%
  summarise(
    N_Targets_Docked     = n(),
    Inhibit_Pool_MMGBSA  = mean(`MMGBSA dG Bind`[Target_Protein %in% INHIBIT_TARGETS],  na.rm = TRUE),
    Activate_Pool_MMGBSA = mean(`MMGBSA dG Bind`[Target_Protein %in% ACTIVATE_TARGETS], na.rm = TRUE),
    Extended_Pool_MMGBSA = mean(`MMGBSA dG Bind`[Target_Protein %in% EXTENDED_TARGETS], na.rm = TRUE),
    All_MPMA_MMGBSA      = mean(`MMGBSA dG Bind`, na.rm = TRUE),
    Best_Single_MMGBSA   = min(`MMGBSA dG Bind`, na.rm = TRUE),
    Best_Single_Target   = Target_Protein[which.min(`MMGBSA dG Bind`)],
    N_Strong_Binders_50  = sum(`MMGBSA dG Bind` < STRONG_BINDER_MMGBSA, na.rm = TRUE),
    N_Inhibit_Strong     = sum(`MMGBSA dG Bind` < STRONG_BINDER_MMGBSA &
                                 Target_Protein %in% INHIBIT_TARGETS, na.rm = TRUE),
    N_Activate_Strong    = sum(`MMGBSA dG Bind` < STRONG_BINDER_MMGBSA &
                                 Target_Protein %in% ACTIVATE_TARGETS, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================================
# 8. PARETO NON-DOMINATED SORTING
# ============================================================================
pareto_check <- function(objectives) {
  n <- nrow(objectives)
  is_pareto <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!is_pareto[i]) next
    for (j in seq_len(n)) {
      if (i == j) next
      if (all(objectives[j, ] <= objectives[i, ]) &&
          any(objectives[j, ] <  objectives[i, ])) {
        is_pareto[i] <- FALSE; break
      }
    }
  }
  is_pareto
}

pareto_obj <- cascade %>%
  select(Inhibit_Pool_MMGBSA, Activate_Pool_MMGBSA, N_Strong_Binders_50) %>%
  mutate(N_Strong_Binders_50 = -N_Strong_Binders_50) %>%
  as.data.frame()

cascade$Pareto_Optimal <- pareto_check(pareto_obj)

best_threshold <- quantile(cascade$Best_Single_MMGBSA, 0.10, na.rm = TRUE)
cascade <- cascade %>%
  mutate(
    In_Top10pct_AnyTarget = Best_Single_MMGBSA <= best_threshold,
    Lead_Tier = case_when(
      Pareto_Optimal & In_Top10pct_AnyTarget & Drug_Like_Strict == 1 ~ "Tier1_Gold",
      Pareto_Optimal & In_Top10pct_AnyTarget                         ~ "Tier2_Silver",
      Pareto_Optimal                                                 ~ "Tier3_Bronze",
      In_Top10pct_AnyTarget                                          ~ "Tier4_Potent",
      TRUE                                                           ~ "Tier5_Other"
    )
  ) %>%
  arrange(Lead_Tier, Inhibit_Pool_MMGBSA)

cat("\nLead tier distribution:\n")
print(table(cascade$Lead_Tier))
cat("\nTop leads:\n")
print(cascade %>%
        filter(Lead_Tier %in% c("Tier1_Gold", "Tier2_Silver", "Tier3_Bronze")) %>%
        select(p_code, Compound, Chem_Class, IUPAC_Methoxy,
               Inhibit_Pool_MMGBSA, Activate_Pool_MMGBSA,
               N_Strong_Binders_50, Best_Single_Target,
               Drug_Like_Strict, Lead_Tier))

# ============================================================================
# 9. SELECTIVITY + PROMISCUITY
# ============================================================================
selectivity <- master %>%
  filter(Is_Control == 0) %>%
  group_by(p_code) %>%
  mutate(
    Mean_OffTarget_MMGBSA = (sum(`MMGBSA dG Bind`, na.rm = TRUE) - `MMGBSA dG Bind`) /
      (n() - 1),
    Selectivity_Index = `MMGBSA dG Bind` / Mean_OffTarget_MMGBSA
  ) %>%
  ungroup() %>%
  select(p_code, Target_Protein, `MMGBSA dG Bind`,
         Mean_OffTarget_MMGBSA, Selectivity_Index)

top3 <- master %>%
  filter(Is_Control == 0) %>%
  group_by(Target_Protein) %>%
  slice_min(order_by = `MMGBSA dG Bind`, n = 3) %>%
  ungroup()

promiscuity <- top3 %>%
  count(p_code, name = "N_Targets_InTop3") %>%
  mutate(
    Promiscuity_Class = case_when(
      N_Targets_InTop3 >= 5 ~ "Highly_Promiscuous",
      N_Targets_InTop3 >= 3 ~ "Moderately_Promiscuous",
      N_Targets_InTop3 >= 1 ~ "Target_Selective",
      TRUE                  ~ "Non_Binder"
    )
  )

cascade <- cascade %>%
  left_join(promiscuity, by = "p_code") %>%
  mutate(
    N_Targets_InTop3  = ifelse(is.na(N_Targets_InTop3),  0L,          N_Targets_InTop3),
    Promiscuity_Class = ifelse(is.na(Promiscuity_Class), "Non_Binder", Promiscuity_Class)
  )

# ============================================================================
# 10. PER-TARGET STATISTICS (Mann-Whitney + Cliff's delta)
# ============================================================================
cat("\n=== PER-TARGET STATISTICS ===\n")

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  nlt <- sum(outer(x, y, "<"))
  ngt <- sum(outer(x, y, ">"))
  (ngt - nlt) / (length(x) * length(y))
}

target_stats <- master %>%
  filter(Is_Control == 0, !is.na(Chem_Class)) %>%
  group_by(Target_Protein) %>%
  summarise(
    N_Total       = n(),
    N_PMF         = sum(Is_PMF == 1, na.rm = TRUE),
    N_NonPMF      = sum(Is_PMF == 0, na.rm = TRUE),
    Median_PMF    = median(`MMGBSA dG Bind`[Is_PMF == 1], na.rm = TRUE),
    Median_NonPMF = median(`MMGBSA dG Bind`[Is_PMF == 0], na.rm = TRUE),
    MannWhitney_p = tryCatch(
      wilcox.test(`MMGBSA dG Bind`[Is_PMF == 1],
                  `MMGBSA dG Bind`[Is_PMF == 0])$p.value,
      error = function(e) NA_real_
    ),
    Cliffs_Delta  = cliffs_delta(`MMGBSA dG Bind`[Is_PMF == 1],
                                 `MMGBSA dG Bind`[Is_PMF == 0]),
    Best_Compound = p_code[which.min(`MMGBSA dG Bind`)],
    Best_MMGBSA   = min(`MMGBSA dG Bind`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    PMF_Significant = as.integer(MannWhitney_p < 0.05),
    Effect_Size_Cat = case_when(
      abs(Cliffs_Delta) >= 0.474 ~ "Large",
      abs(Cliffs_Delta) >= 0.33  ~ "Medium",
      abs(Cliffs_Delta) >= 0.147 ~ "Small",
      TRUE                        ~ "Negligible"
    ),
    Therapeutic_Mode = case_when(
      Target_Protein %in% INHIBIT_TARGETS  ~ "Inhibit",
      Target_Protein %in% ACTIVATE_TARGETS ~ "Activate",
      Target_Protein %in% EXTENDED_TARGETS ~ "Extended",
      TRUE                                  ~ "Unknown"
    )
  )

print(target_stats)

# ============================================================================
# 11. CLASSIFICATION AUDIT
# ============================================================================
classification_audit <- compounds %>%
  select(p_code, Compound, Official_Name,
         Stored_Methoxy, IUPAC_Methoxy, Methoxy_Corrected,
         Chem_Class, Is_PMF, Is_Gingeroid,
         Drug_Like_Strict, Drug_Like_Lenient) %>%
  arrange(p_code)

# ============================================================================
# 12. ML FEATURE SCHEMA
# ============================================================================
feature_schema <- tribble(
  ~Column,                        ~Role,      ~Use_In_ML,   ~Notes,
  "Target_Protein",               "ID",       "GROUP",      "GroupKFold groups",
  "Title",                        "ID",       "EXCLUDE",    "PubChem CID",
  "Entry Name",                   "ID",       "EXCLUDE",    "Pose ID",
  "p_code",                       "ID",       "GROUP",      "Compound group for LOCO-CV",
  "Compound",                     "ID",       "EXCLUDE",    "Common name",
  "Official_Name",                "ID",       "EXCLUDE",    "IUPAC name",
  "Pub_CID",                      "ID",       "EXCLUDE",    "PubChem CID",
  "smiles",                       "ID",       "EXCLUDE",    "SMILES string",
  "MMGBSA dG Bind",               "OUTCOME",  "Y_PRIMARY",  "Primary regression target",
  "XP_GScore",                    "OUTCOME",  "Y_SECONDARY","Secondary validation target",
  "M3",  "FEATURE", "USE", "Methoxy at 3 (C-ring)",
  "M5",  "FEATURE", "USE", "Methoxy at 5 (A-ring)",
  "M6",  "FEATURE", "USE", "Methoxy at 6 (A-ring)",
  "M7",  "FEATURE", "USE", "Methoxy at 7 (A-ring)",
  "M8",  "FEATURE", "USE", "Methoxy at 8 (A-ring)",
  "M3p", "FEATURE", "USE", "Methoxy at 3' (B-ring)",
  "M4p", "FEATURE", "USE", "Methoxy at 4' (B-ring)",
  "M5p", "FEATURE", "USE", "Methoxy at 5' (B-ring)",
  "OH3",  "FEATURE", "USE", "Hydroxyl at 3",
  "OH5",  "FEATURE", "USE", "Hydroxyl at 5 (5-OH/C4=O intramol HB)",
  "OH7",  "FEATURE", "USE", "Hydroxyl at 7",
  "OH3p", "FEATURE", "USE", "Hydroxyl at 3'",
  "OH4p", "FEATURE", "USE", "Hydroxyl at 4'",
  "Is_Sugar",                 "FEATURE", "USE", "Glycosylation flag",
  "Sugar_OH_Count",           "FEATURE", "USE", "Sugar hydroxyl count",
  "Alkyl_Chain_Len",          "FEATURE", "USE", "Gingeroid alkyl chain",
  "Has_5OH_IntramolecularHB", "FEATURE", "USE", "5-OH/C4=O flag",
  "LogP",                     "FEATURE", "USE", "Consensus LogP",
  "TPSA",                     "FEATURE", "USE", "Topological polar surface area",
  "#Heavy atoms",             "FEATURE", "USE", "MW proxy",
  "#Rotatable bonds",         "FEATURE", "USE", "Flexibility",
  "#H-bond acceptors",        "FEATURE", "USE", "HBA",
  "#H-bond donors",           "FEATURE", "USE", "HBD",
  "Fraction Csp3",            "FEATURE", "USE", "sp3 fraction",
  "MR",                       "FEATURE", "USE", "Molar refractivity",
  "Synthetic Accessibility",  "FEATURE", "USE", "SA score",
  "Total_RuleOf5_Violations", "FEATURE", "USE", "Summed drug-likeness violations",
  "Chem_Class",        "METADATA", "EXCLUDE",   "Circular with IUPAC_Methoxy",
  "Is_PMF",            "METADATA", "CLASS_Y",   "Binary classification target",
  "Is_Gingeroid",      "METADATA", "EXCLUDE",   "Subclass",
  "Is_Flavonoid",      "METADATA", "EXCLUDE",   "Subclass",
  "Is_Control",        "METADATA", "FILTER",    "Filter out PC rows",
  "IUPAC_Methoxy",     "METADATA", "EXCLUDE",   "Defines PMF — LEAKAGE",
  "IUPAC_Hydroxy",     "METADATA", "EXCLUDE",   "Implicit in positional OH",
  "Drug_Like_Strict",  "METADATA", "SUBSET",    "Tier-2 ML filter",
  "Drug_Like_Lenient", "METADATA", "SUBSET",    "Tier-2 ML filter (permissive)",
  "Inhibit_Pool_MMGBSA",  "DERIVED", "EXCLUDE",   "Leakage",
  "Activate_Pool_MMGBSA", "DERIVED", "EXCLUDE",   "Leakage",
  "Pareto_Optimal",       "DERIVED", "CASCADE_Y", "Cascade classifier target",
  "Promiscuity_Class",    "DERIVED", "EXCLUDE",   "Post-hoc label"
)

# ============================================================================
# 13. DRUG-LIKE SUBSET
# ============================================================================
drug_like_subset <- master %>%
  filter(Is_Control == 0, Drug_Like_Strict == 1)

cat(sprintf("\nDrug-like subset: %d rows (%d compounds x %d targets)\n",
            nrow(drug_like_subset),
            length(unique(drug_like_subset$p_code)),
            length(unique(drug_like_subset$Target_Protein))))

excluded <- compounds %>%
  filter(Drug_Like_Strict == 0) %>%
  select(p_code, Compound, Chem_Class, TPSA, Lipinski_Violations, Veber_Violations,
         PAINS_Alerts, Brenk_Alerts)
cat("\nCompounds excluded from drug-like subset:\n")
print(excluded)

# ============================================================================
# 14. WRITE OUTPUTS
# ============================================================================
cat("\n=== WRITING OUTPUT FILES ===\n")

write.csv(master,               file.path(OUT_DIR, "FINAL_MASTER_AI_DATASET.csv"),  row.names = FALSE)
write.csv(feature_schema,       file.path(OUT_DIR, "ML_FEATURE_SCHEMA.csv"),        row.names = FALSE)
write.csv(cascade,              file.path(OUT_DIR, "COMPOUND_CASCADE_RANKING.csv"), row.names = FALSE)
write.csv(drug_like_subset,     file.path(OUT_DIR, "DRUGLIKE_SUBSET.csv"),          row.names = FALSE)
write.csv(target_stats,         file.path(OUT_DIR, "TARGET_ANALYSIS_SUMMARY.csv"),  row.names = FALSE)
write.csv(classification_audit, file.path(OUT_DIR, "CLASSIFICATION_AUDIT.csv"),     row.names = FALSE)
write.csv(selectivity,          file.path(OUT_DIR, "SELECTIVITY_LONG.csv"),         row.names = FALSE)

cat("\nAll 7 output files written.\n")
cat(sprintf("Master dataset:   %d rows, %d cols\n", nrow(master), ncol(master)))
cat(sprintf("Pareto-optimal:   %d compounds\n", sum(cascade$Pareto_Optimal)))
cat(sprintf("Drug-like strict: %d compounds\n", sum(compounds$Drug_Like_Strict)))
cat("\n>>> Next: run run_ml_pipeline.py on FINAL_MASTER_AI_DATASET.csv <<<\n")

summary(drug_like_subset)

