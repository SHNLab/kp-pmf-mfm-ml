"""
================================================================================
ML.py  —  ML VALIDATION + FIGURE GENERATION
--------------------------------------------------------------------------------


[ANALYSIS STACK]
  1. Nested cross-validation  (outer LOCO, inner 5-fold grid search on XGBoost)
  2. Conformal prediction intervals  (distribution-free 95% coverage)
  3. Y-randomization              (500 permutations per target)
  4. Gaussian Process Regression  (calibrated uncertainty)
  5. Applicability domain         (Williams-plot KNN distance)
  6. Baseline benchmarks          (physchem vs position vs full features)
  7. Learning curve               (R² vs training fraction)
  8. Active learning simulation   (GPR uncertainty vs random acquisition)

[FIGURES AUTO-GENERATED]
  FigV5A_y_randomization.{png,svg}
  FigV5B_applicability_domain.{png,svg}
  FigV5C_gpr_uncertainty.{png,svg}
  FigV5D_learning_curve.{png,svg}
  FigV5E_active_learning.{png,svg}
  FigV5F_baseline_benchmarks.{png,svg}
  FigV5G_conformal_coverage.{png,svg}           [NEW — conformal validation]
  FigV5H_nested_cv_comparison.{png,svg}         [NEW — nested vs single-fold]
  FigV5_hero_validation.{png,svg}               [4-panel summary]

[OUTPUTS — CSVs + XLSX]
  NESTED_CV_RESULTS.csv
  CONFORMAL_PREDICTIONS.csv
  HYPERPARAMETER_SCAN.csv
  Y_RANDOMIZATION_RESULTS.csv
  APPLICABILITY_DOMAIN.csv
  GPR_UNCERTAINTY_PREDICTIONS.csv
  BASELINE_BENCHMARKS.csv
  LEARNING_CURVE.csv
  ACTIVE_LEARNING_SIMULATION.csv
  ML_VALIDATION_COMPLETE.xlsx

READS:
  FINAL_MASTER_AI_DATASET.csv
  COMPOUND_CASCADE_RANKING.csv
  CLASSIFICATION_AUDIT.csv

================================================================================
"""

import os, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.linear_model import LassoCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict, KFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.spatial.distance import pdist, squareform

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("xgboost not found — install with: pip install xgboost")
    exit(1)

import matplotlib
matplotlib.use('Agg')   # non-interactive backend (safer on Windows/headless)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns

RNG_SEED = 42
rng = np.random.default_rng(RNG_SEED)

# ============================================================================
# FIGURE STYLE
# ============================================================================
OUT_FIG = Path("figures_v5"); OUT_FIG.mkdir(exist_ok=True)

COLORS = {
    'PMF':          '#c0392b',
    'Meth_Flav':    '#2874a6',
    'Gingeroid':    '#1e8449',
    'Glycoside':    '#d68910',
    'Standard_Flav':'#7d3c98',
    'real':         '#c0392b',
    'null':         '#bdc3c7',
    'baseline':     '#566573',
    'random':       '#95a5a6',
    'model':        '#c0392b',
    'uncertainty':  '#e67e22',
    'in_domain':    '#27ae60',
    'out_domain':   '#c0392b',
    'neutral':      '#566573',
    'conformal':    '#8e44ad',
    'nested':       '#16a085',
}

# Publication-grade typography for journal submission
# Times New Roman 12pt base; sub-axes elements scaled proportionally
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Liberation Serif', 'DejaVu Serif', 'serif'],
    'font.size': 12,
    'axes.titlesize': 13,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'axes.labelweight': 'normal',
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'legend.title_fontsize': 12,
    'figure.titlesize': 14,
    'figure.titleweight': 'bold',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#2c3e50', 'axes.linewidth': 1.0,
    'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.major.size': 4, 'ytick.major.size': 4,
    'lines.linewidth': 1.5,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'mathtext.fontset': 'stix',  # serif math to match Times
})

# ============================================================================
# CONFIG
# ============================================================================
POSITIONAL = ['M3','M5','M6','M7','M8','M3p','M4p','M5p',
              'OH3','OH5','OH7','OH3p','OH4p']
STRUCTURAL = ['Is_Sugar','Sugar_OH_Count','Alkyl_Chain_Len',
              'Has_5OH_IntramolecularHB']
PHYSCHEM = ['LogP','TPSA','#Heavy atoms','#Rotatable bonds',
            '#H-bond acceptors','#H-bond donors','Fraction Csp3','MR',
            'Synthetic Accessibility']
ALL_FEATURES = POSITIONAL + STRUCTURAL + PHYSCHEM

# 9-target MPMA cascade panel — final locked panel
# Tier 1 (6, above Q3 + druggability ≥ 7.5):
#   BCL2, HSP90, PINK1, HK2, CPT1A, MFN2
# Tier 2 (3, below Q3 but cascade-obligate — post-translationally regulated):
#   DRP1, VDAC1, SIRT3
#
# REMOVED from earlier drafts:
#   - MMP9 (not mitochondrial by localisation or function)
#   - SIRT1 (nuclear-predominant; contested CRC biology; not cascade-obligate)
#   - KEAP1 (cysteine-oxidation regulated; below Q3; not cascade-obligate)
TARGETS_KEEP = ['BCL2','HSP90','PINK1','HK2','CPT1A','MFN2',
                'DRP1','VDAC1','SIRT3']

# Hyperparameter grid for nested CV (kept small — small dataset)
XGB_PARAM_GRID = {
    'n_estimators': [150, 300, 500],
    'max_depth':    [2, 3, 4],
    'learning_rate':[0.03, 0.05, 0.10],
    'subsample':    [0.7, 0.85]
}

N_Y_PERM      = 500       # Y-randomization
AD_K          = 5          # KNN for applicability domain
LC_FRACTIONS  = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
LC_REPS       = 10
AL_N_INITIAL  = 5
AL_REPS       = 20
CONFORMAL_ALPHA = 0.05    # 95% coverage target

# ============================================================================
# LOAD + CHEMOTYPE ASSIGNMENT (self-healing)
# ============================================================================
print("="*75)
print("ML PIPELINE V4 — nested CV + conformal prediction")
print("="*75)

# Prefer the corrected positional-indicator dataset; fall back to original
import os
if os.path.exists("FINAL_MASTER_AI_DATASET_corrected.csv"):
    master  = pd.read_csv("FINAL_MASTER_AI_DATASET_corrected.csv")
    print("[input] using FINAL_MASTER_AI_DATASET_corrected.csv (positional flags re-derived)")
else:
    master  = pd.read_csv("FINAL_MASTER_AI_DATASET.csv")
    print("[input] using FINAL_MASTER_AI_DATASET.csv  WARNING: positional flags may be stale; "
          "run fix_master_dataset.py first to fix Lasso/regression feature matrix")
cascade = pd.read_csv("COMPOUND_CASCADE_RANKING.csv")

GINGEROIDS_SET = ['p33','p34','p35','p36','p37']
GLYCOSIDES_SET = ['p1','p4','p5','p7']
PMF_METH_MIN   = 2

def _assign_chemotype(p_code, iupac_meth):
    if p_code in GINGEROIDS_SET: return 'Gingeroid'
    if p_code in GLYCOSIDES_SET: return 'Glycoside'
    if pd.isna(iupac_meth):      return 'Unknown'
    if iupac_meth >= PMF_METH_MIN: return 'PMF'
    if iupac_meth == 1:            return 'Meth_Flav'
    if iupac_meth == 0:            return 'Standard_Flav'
    return 'Unknown'

def _build_chemo_map():
    for fname in ('CLASSIFICATION_AUDIT.csv',
                  'FINAL_MASTER_AI_DATASET.csv',
                  'COMPOUND_CASCADE_RANKING.csv'):
        if not os.path.exists(fname): continue
        df = pd.read_csv(fname)
        if 'p_code' in df.columns and 'IUPAC_Methoxy' in df.columns:
            df_u = df[['p_code','IUPAC_Methoxy']].drop_duplicates('p_code')
            return {r['p_code']: _assign_chemotype(r['p_code'], r['IUPAC_Methoxy'])
                    for _, r in df_u.iterrows()}
    return {}

if 'Chemotype' in cascade.columns:
    chemo_map = dict(zip(cascade['p_code'], cascade['Chemotype']))
elif 'Chemotype' in master.columns:
    chemo_map = dict(zip(master['p_code'], master['Chemotype']))
else:
    chemo_map = _build_chemo_map()
    if not chemo_map:
        raise ValueError("Need CLASSIFICATION_AUDIT.csv with IUPAC_Methoxy")

if 'Chemotype' not in master.columns:
    master['Chemotype'] = master['p_code'].map(chemo_map)
if 'Chemotype' not in cascade.columns:
    cascade['Chemotype'] = cascade['p_code'].map(chemo_map)

plant = master[master['Is_Control']==0].dropna(
    subset=['MMGBSA dG Bind']+ALL_FEATURES).copy()
plant = plant[plant['Target_Protein'].isin(TARGETS_KEEP)]
print(f"Plant data (9 targets): {len(plant)} rows, {plant['p_code'].nunique()} compounds")

# ============================================================================
# HELPER — get standardized features per target
# ============================================================================
def _prep_target_data(tgt):
    sub = plant[plant['Target_Protein']==tgt].reset_index(drop=True)
    if len(sub) < 10: return None
    X = sub[ALL_FEATURES].values
    y = sub['MMGBSA dG Bind'].values
    groups = sub['p_code'].values
    stds = X.std(axis=0)
    keep = stds > 1e-8
    if keep.sum() == 0: return None
    X_used = X[:, keep]
    feats_used = [f for f, k in zip(ALL_FEATURES, keep) if k]
    Xs = StandardScaler().fit_transform(X_used)
    return sub, Xs, y, groups, feats_used

# ============================================================================
# 1. NESTED CROSS-VALIDATION  [V4 NEW]
# ============================================================================
# Outer LOCO-CV evaluates generalization. Inner 5-fold tunes hyperparameters.
# This prevents the "hyperparameters overfit to test compounds" concern.
# ============================================================================
print("\n[1] NESTED CROSS-VALIDATION (outer=LOCO, inner=5-fold grid search)…")

nested_rows = []
hp_rows = []

for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None:
        print(f"  {tgt}: skipped (insufficient data)")
        continue
    sub, Xs, y, groups, feats = prep
    logo = LeaveOneGroupOut()
    n_splits = len(np.unique(groups))

    # Track outer-fold predictions and best hyperparameters
    outer_preds = np.zeros(len(y))
    outer_best_params = []

    for fold_idx, (tr_idx, te_idx) in enumerate(logo.split(Xs, y, groups)):
        X_tr, y_tr = Xs[tr_idx], y[tr_idx]
        X_te, y_te = Xs[te_idx], y[te_idx]

        # Inner CV — grid search on training set only
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=RNG_SEED)
        grid = GridSearchCV(
            xgb.XGBRegressor(random_state=RNG_SEED, verbosity=0,
                             tree_method='hist', n_jobs=1),
            param_grid=XGB_PARAM_GRID,
            cv=inner_cv,
            scoring='neg_mean_absolute_error',
            n_jobs=-1,
            refit=True
        )
        try:
            grid.fit(X_tr, y_tr)
            outer_preds[te_idx] = grid.predict(X_te)
            outer_best_params.append(grid.best_params_)
        except Exception as e:
            outer_preds[te_idx] = y_tr.mean()
            outer_best_params.append({})

    r2_nested = r2_score(y, outer_preds)
    mae_nested = mean_absolute_error(y, outer_preds)

    # Most-common hyperparameters
    if outer_best_params:
        hp_df = pd.DataFrame(outer_best_params)
        modal_hp = {col: hp_df[col].mode().iloc[0] for col in hp_df.columns}
    else:
        modal_hp = {}

    nested_rows.append({
        'Target': tgt,
        'R2_Nested_CV': r2_nested,
        'MAE_Nested_CV': mae_nested,
        'N_outer_folds': n_splits
    })
    hp_rows.append({'Target': tgt, **modal_hp})

    print(f"  {tgt:6s}: Nested R²={r2_nested:+.3f}  MAE={mae_nested:.2f}  "
          f"modal depth={modal_hp.get('max_depth','NA')}  "
          f"modal lr={modal_hp.get('learning_rate','NA')}")

nested_df = pd.DataFrame(nested_rows)
hp_df = pd.DataFrame(hp_rows)
nested_df.to_csv('NESTED_CV_RESULTS.csv', index=False)
hp_df.to_csv('HYPERPARAMETER_SCAN.csv', index=False)

# ============================================================================
# 2. CONFORMAL PREDICTION INTERVALS  [V4 NEW]
# ============================================================================
# Distribution-free uncertainty with coverage guarantee:
# 1. Split data into training + calibration
# 2. Train XGBoost on training set
# 3. Compute absolute residuals on calibration set
# 4. Quantile of residuals = conformal interval width at 1-alpha coverage
# Coverage guarantee holds regardless of true data distribution.
# ============================================================================
print("\n[2] CONFORMAL PREDICTION INTERVALS…")

conformal_rows = []

for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None: continue
    sub, Xs, y, groups, feats = prep
    logo = LeaveOneGroupOut()

    # For each outer LOCO fold: split train further into train+calibration
    for fold_idx, (tr_idx, te_idx) in enumerate(logo.split(Xs, y, groups)):
        X_tr_all, y_tr_all = Xs[tr_idx], y[tr_idx]
        X_te, y_te = Xs[te_idx], y[te_idx]
        n_tr_all = len(y_tr_all)
        if n_tr_all < 10: continue

        # 70/30 split of training into train/calibration
        n_cal = max(int(0.3 * n_tr_all), 5)
        cal_idx = rng.choice(n_tr_all, n_cal, replace=False)
        fit_idx = np.array([i for i in range(n_tr_all) if i not in cal_idx])

        model = xgb.XGBRegressor(
            n_estimators=300, max_depth=3, learning_rate=0.05,
            subsample=0.8, random_state=RNG_SEED, verbosity=0,
            tree_method='hist', n_jobs=-1
        )
        model.fit(X_tr_all[fit_idx], y_tr_all[fit_idx])

        # Calibration residuals
        cal_preds = model.predict(X_tr_all[cal_idx])
        cal_resid = np.abs(y_tr_all[cal_idx] - cal_preds)

        # Quantile at (1 - alpha) with finite-sample correction
        q_level = np.ceil((n_cal + 1) * (1 - CONFORMAL_ALPHA)) / n_cal
        q_level = min(q_level, 1.0)
        q_val = np.quantile(cal_resid, q_level)

        # Predict on test set with conformal interval
        test_preds = model.predict(X_te)
        for i_te, test_i in enumerate(te_idx):
            pred = test_preds[i_te]
            true = y_te[i_te]
            lower = pred - q_val
            upper = pred + q_val
            in_ci = (true >= lower) and (true <= upper)
            conformal_rows.append({
                'Target': tgt,
                'Compound': sub['p_code'].iloc[test_i],
                'Chemotype': sub['Chemotype'].iloc[test_i],
                'y_true': true,
                'y_pred': pred,
                'CI_lower': lower,
                'CI_upper': upper,
                'Interval_Width': upper - lower,
                'Within_CI': in_ci,
                'Calibration_N': n_cal
            })

conformal_df = pd.DataFrame(conformal_rows)
conformal_df.to_csv('CONFORMAL_PREDICTIONS.csv', index=False)

# Summary per target
conf_perf = conformal_df.groupby('Target').agg(
    Coverage=('Within_CI', 'mean'),
    Mean_Width=('Interval_Width', 'mean'),
    N=('Within_CI', 'count')
).reset_index()
conf_perf['Coverage'] = conf_perf['Coverage'].round(3)
conf_perf['Mean_Width'] = conf_perf['Mean_Width'].round(2)
conf_perf['Coverage_Target'] = 1 - CONFORMAL_ALPHA
print(conf_perf.to_string(index=False))
print(f"\nOverall conformal coverage: {conformal_df['Within_CI'].mean()*100:.1f}% "
      f"(target: {(1-CONFORMAL_ALPHA)*100:.0f}%)")

# ============================================================================
# 3. Y-RANDOMIZATION (carried from v3)
# ============================================================================
print("\n[3] Y-RANDOMIZATION (500 permutations per target)…")

yr_rows = []
for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None: continue
    sub, Xs, y, groups, feats = prep
    logo = LeaveOneGroupOut()
    model = xgb.XGBRegressor(n_estimators=300, max_depth=3, learning_rate=0.05,
                              subsample=0.8, random_state=RNG_SEED, verbosity=0,
                              tree_method='hist', n_jobs=-1)

    pred_real = cross_val_predict(model, Xs, y, cv=logo.split(Xs, y, groups))
    r2_real = r2_score(y, pred_real)

    null_r2 = np.zeros(N_Y_PERM)
    for i in range(N_Y_PERM):
        y_sh = rng.permutation(y)
        p_null = cross_val_predict(model, Xs, y_sh, cv=logo.split(Xs, y_sh, groups))
        null_r2[i] = r2_score(y_sh, p_null)

    p_val = (np.sum(null_r2 >= r2_real) + 1) / (N_Y_PERM + 1)
    yr_rows.append({
        'Target': tgt, 'R2_real': r2_real,
        'Null_mean': null_r2.mean(), 'Null_SD': null_r2.std(),
        'Null_95th': np.percentile(null_r2, 95),
        'Null_99th': np.percentile(null_r2, 99),
        'p_value': p_val, 'Significant': p_val < 0.05
    })
    print(f"  {tgt:6s}: R²_real={r2_real:+.3f}  p={p_val:.4f}  "
          f"{'SIGNIFICANT' if p_val<0.05 else 'ns'}")

yr_df = pd.DataFrame(yr_rows)
yr_df.to_csv('Y_RANDOMIZATION_RESULTS.csv', index=False)

# ============================================================================
# 4. GAUSSIAN PROCESS REGRESSION (carried from v3)
# ============================================================================
print("\n[4] GPR UNCERTAINTY…")

gpr_rows = []
for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None: continue
    sub, Xs, y, groups, feats = prep
    kernel = ConstantKernel(1.0) * Matern(1.0, nu=2.5) + WhiteKernel(1.0)
    gpr = GaussianProcessRegressor(kernel=kernel, alpha=1e-6, random_state=RNG_SEED,
                                    n_restarts_optimizer=3, normalize_y=True)
    logo = LeaveOneGroupOut()
    preds = np.zeros(len(y)); stds = np.zeros(len(y))

    for tr_idx, te_idx in logo.split(Xs, y, groups):
        try:
            gpr.fit(Xs[tr_idx], y[tr_idx])
            mu, sigma = gpr.predict(Xs[te_idx], return_std=True)
            preds[te_idx] = mu; stds[te_idx] = sigma
        except Exception:
            preds[te_idx] = y[tr_idx].mean()
            stds[te_idx] = y.std()

    in_ci = np.abs(y - preds) <= 1.96 * stds
    r2 = r2_score(y, preds)
    mae = mean_absolute_error(y, preds)
    cov = in_ci.mean()

    for i in range(len(y)):
        gpr_rows.append({
            'Target': tgt, 'Compound': sub['p_code'].iloc[i],
            'Chemotype': sub['Chemotype'].iloc[i],
            'y_true': y[i], 'y_pred': preds[i], 'sigma': stds[i],
            'CI_low': preds[i] - 1.96*stds[i],
            'CI_high': preds[i] + 1.96*stds[i],
            'Within_CI': in_ci[i]
        })
    print(f"  {tgt:6s}: R²={r2:+.3f}  MAE={mae:.2f}  coverage={cov*100:.0f}%")

gpr_df = pd.DataFrame(gpr_rows)
gpr_df.to_csv('GPR_UNCERTAINTY_PREDICTIONS.csv', index=False)

# ============================================================================
# 5. APPLICABILITY DOMAIN (carried from v3)
# ============================================================================
print("\n[5] APPLICABILITY DOMAIN…")

ad_rows = []
for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None: continue
    sub, Xs, y, groups, feats = prep
    dist = squareform(pdist(Xs, 'euclidean'))
    np.fill_diagonal(dist, np.inf)
    k = min(AD_K, len(sub)-1)
    nn = np.sort(dist, axis=1)[:, :k].mean(axis=1)
    threshold = nn.mean() + 3*nn.std()
    for i in range(len(sub)):
        ad_rows.append({
            'Target': tgt, 'Compound': sub['p_code'].iloc[i],
            'Chemotype': sub['Chemotype'].iloc[i],
            'Mean_Distance_to_KNN': nn[i],
            'AD_Threshold': threshold,
            'In_Domain': nn[i] <= threshold
        })
ad_df = pd.DataFrame(ad_rows)
ad_df.to_csv('APPLICABILITY_DOMAIN.csv', index=False)
print(f"  Overall OOD: {(1-ad_df['In_Domain'].mean())*100:.1f}%")

# ============================================================================
# 6. BASELINE BENCHMARKS (carried from v3)
# ============================================================================
print("\n[6] BASELINE BENCHMARKS…")

baseline_rows = []
for tgt in TARGETS_KEEP:
    sub = plant[plant['Target_Protein']==tgt]
    if len(sub) < 10: continue
    y = sub['MMGBSA dG Bind'].values
    groups = sub['p_code'].values
    logo = LeaveOneGroupOut()

    def _eval(X_mat):
        stds = X_mat.std(axis=0)
        keep = stds > 1e-8
        if keep.sum() == 0: return np.nan
        X = X_mat[:, keep]
        Xs = StandardScaler().fit_transform(X)
        m = xgb.XGBRegressor(n_estimators=300, max_depth=3, learning_rate=0.05,
                              subsample=0.8, random_state=RNG_SEED, verbosity=0,
                              tree_method='hist', n_jobs=-1)
        pred = cross_val_predict(m, Xs, y, cv=logo.split(Xs, y, groups))
        return r2_score(y, pred)

    r2_full = _eval(sub[ALL_FEATURES].values)
    r2_phys = _eval(sub[PHYSCHEM].values)
    r2_pos  = _eval(sub[POSITIONAL].values)
    r2_rand = r2_score(y, np.full_like(y, y.mean(), dtype=float))

    baseline_rows.append({
        'Target': tgt,
        'R2_Random_Mean': r2_rand,
        'R2_Physchem_Only': r2_phys,
        'R2_Position_Only': r2_pos,
        'R2_Full_Features': r2_full,
        'Delta_Position_vs_Physchem': r2_pos - r2_phys,
        'Delta_Full_vs_Physchem': r2_full - r2_phys
    })
    print(f"  {tgt:6s}: Phys={r2_phys:+.3f} Pos={r2_pos:+.3f} Full={r2_full:+.3f}")

baseline_df = pd.DataFrame(baseline_rows)
baseline_df.to_csv('BASELINE_BENCHMARKS.csv', index=False)

# ============================================================================
# 7. LEARNING CURVE (carried from v3)
# ============================================================================
print("\n[7] LEARNING CURVE…")

lc_rows = []
for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None: continue
    sub, Xs, y, groups, feats = prep
    n = len(y)

    for frac in LC_FRACTIONS:
        n_train = max(int(frac * n), 5)
        if n_train >= n: continue
        r2s = []
        for rep in range(LC_REPS):
            idx = rng.permutation(n)
            tr, te = idx[:n_train], idx[n_train:]
            m = xgb.XGBRegressor(n_estimators=300, max_depth=3, learning_rate=0.05,
                                  subsample=0.8, random_state=rep, verbosity=0,
                                  tree_method='hist', n_jobs=-1)
            m.fit(Xs[tr], y[tr])
            pred = m.predict(Xs[te])
            r2s.append(r2_score(y[te], pred))
        lc_rows.append({
            'Target': tgt, 'Train_Fraction': frac, 'N_train': n_train,
            'R2_mean': np.mean(r2s), 'R2_std': np.std(r2s)
        })
lc_df = pd.DataFrame(lc_rows)
lc_df.to_csv('LEARNING_CURVE.csv', index=False)

# ============================================================================
# 8. ACTIVE LEARNING (carried from v3)
# ============================================================================
print("\n[8] ACTIVE LEARNING…")

al_rows = []
for tgt in TARGETS_KEEP:
    prep = _prep_target_data(tgt)
    if prep is None: continue
    sub, Xs, y, groups, feats = prep
    n = len(y)

    for rep in range(AL_REPS):
        seed = rep * 1000 + 1

        # Uncertainty strategy
        avail = set(range(n))
        acq = list(rng.choice(list(avail), AL_N_INITIAL, replace=False))
        avail -= set(acq)
        step = 0
        while avail:
            k = ConstantKernel(1.0) * Matern(1.0, nu=2.5) + WhiteKernel(1.0)
            gpr = GaussianProcessRegressor(kernel=k, alpha=1e-6, random_state=seed,
                                            normalize_y=True, n_restarts_optimizer=2)
            try:
                gpr.fit(Xs[acq], y[acq])
                avl = list(avail)
                mu, sigma = gpr.predict(Xs[avl], return_std=True)
                r2_h = r2_score(y[avl], mu) if len(avl) > 1 else np.nan
            except Exception:
                r2_h = np.nan
                sigma = np.ones(len(avail))
                avl = list(avail)
            al_rows.append({
                'Target': tgt, 'Rep': rep, 'Strategy': 'UncertaintySampling',
                'Step': step, 'N_acquired': len(acq),
                'N_acquired_pct': len(acq)/n*100,
                'R2_on_remaining': r2_h
            })
            if not avl: break
            next_i = avl[int(np.argmax(sigma))]
            acq.append(next_i); avail.discard(next_i); step += 1

        # Random strategy
        avail = set(range(n))
        acq = list(rng.choice(list(avail), AL_N_INITIAL, replace=False))
        avail -= set(acq)
        step = 0
        while avail:
            k = ConstantKernel(1.0) * Matern(1.0, nu=2.5) + WhiteKernel(1.0)
            gpr = GaussianProcessRegressor(kernel=k, alpha=1e-6, random_state=seed,
                                            normalize_y=True, n_restarts_optimizer=2)
            try:
                gpr.fit(Xs[acq], y[acq])
                avl = list(avail)
                mu = gpr.predict(Xs[avl])
                r2_h = r2_score(y[avl], mu) if len(avl) > 1 else np.nan
            except Exception:
                r2_h = np.nan
                avl = list(avail)
            al_rows.append({
                'Target': tgt, 'Rep': rep, 'Strategy': 'Random',
                'Step': step, 'N_acquired': len(acq),
                'N_acquired_pct': len(acq)/n*100,
                'R2_on_remaining': r2_h
            })
            if not avl: break
            next_i = int(rng.choice(list(avail)))
            acq.append(next_i); avail.discard(next_i); step += 1
    print(f"  {tgt:6s}: complete ({AL_REPS} reps × 2 strategies)")

al_df = pd.DataFrame(al_rows)
al_df.to_csv('ACTIVE_LEARNING_SIMULATION.csv', index=False)

# ============================================================================
# WRITE COMPREHENSIVE XLSX (Python 3.14 safe — keyword args)
# ============================================================================
print("\n[9] Writing ML_VALIDATION_COMPLETE.xlsx…")

try:
    with pd.ExcelWriter('ML_VALIDATION_COMPLETE.xlsx', engine='openpyxl') as w:
        nested_df.to_excel(w, sheet_name='Nested_CV', index=False)
        hp_df.to_excel(w, sheet_name='Best_Hyperparams', index=False)
        conf_perf.to_excel(w, sheet_name='Conformal_Summary', index=False)
        conformal_df.to_excel(w, sheet_name='Conformal_Predictions', index=False)
        yr_df.to_excel(w, sheet_name='Y_Randomization', index=False)
        ad_df.to_excel(w, sheet_name='Applicability_Domain', index=False)
        gpr_df.to_excel(w, sheet_name='GPR_Predictions', index=False)
        baseline_df.to_excel(w, sheet_name='Baselines', index=False)
        lc_df.to_excel(w, sheet_name='Learning_Curve', index=False)
        al_df.groupby(['Target','Strategy','N_acquired_pct'])['R2_on_remaining'].agg(
            ['median','mean','std']).reset_index().to_excel(
                w, sheet_name='Active_Learning_Summary', index=False)
    print("  XLSX written successfully")
except Exception as e:
    print(f"  XLSX failed: {e}")
    print("  All data still available as CSV files")


# ============================================================================
# 8. CLASSIFICATION AUC + BOOTSTRAP 95% CI  [M1, peer-review revision]
# ============================================================================
# Reframes AUC point estimates with bootstrap CI to expose small-sample
# uncertainty. AUC near 1.0 at n=37 is presented with explicit CI bounds.
# Reference: Carpenter & Bithell (2000) bootstrap CI for AUC.
# ============================================================================
print("\n[8] CLASSIFICATION AUC with bootstrap 95% CI...")

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

# Build wide-form binding-profile matrix (compound × target) for classification
profile = plant.pivot_table(index='p_code', columns='Target_Protein',
                            values='MMGBSA dG Bind', aggfunc='first')
chemo_lookup = plant.drop_duplicates('p_code').set_index('p_code')['Chemotype']
profile['Chemotype'] = profile.index.map(chemo_lookup)

X_profile = profile[[t for t in TARGETS_KEEP if t in profile.columns]].copy()
X_profile = X_profile.fillna(X_profile.mean()).values
chemo_arr = profile['Chemotype'].values
p_codes_arr = profile.index.values

excluded_set = set(GINGEROIDS_SET + GLYCOSIDES_SET)
classification_tasks = [
    ('Strong vs Weak (median split)',
     (np.nanmean(X_profile, axis=1) < np.nanmedian(np.nanmean(X_profile, axis=1))).astype(int)),
    ('Gingeroid vs non-Gingeroid', (chemo_arr == 'Gingeroid').astype(int)),
    ('Glycoside vs other', (chemo_arr == 'Glycoside').astype(int)),
    ('Drug-like PMF subset',
     ((chemo_arr == 'PMF') & ~np.isin(p_codes_arr, list(excluded_set))).astype(int)),
]

m1_rows = []
N_BOOT = 1000

for task_name, y_bin in classification_tasks:
    n_pos = int(y_bin.sum()); n_neg = len(y_bin) - n_pos
    if n_pos < 3 or n_neg < 3:
        print(f"  {task_name}: SKIP (pos={n_pos}, neg={n_neg})")
        continue

    rf_clf = RandomForestClassifier(n_estimators=400, class_weight='balanced',
                                    random_state=RNG_SEED, n_jobs=-1)
    xgb_clf = xgb.XGBClassifier(max_depth=3, learning_rate=0.05, n_estimators=400,
                                random_state=RNG_SEED, eval_metric='logloss',
                                tree_method='hist', n_jobs=-1)
    cv_strat = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG_SEED)
    try:
        rf_proba = cross_val_predict(rf_clf, X_profile, y_bin, cv=cv_strat,
                                     method='predict_proba', n_jobs=-1)[:,1]
        auc_rf = roc_auc_score(y_bin, rf_proba)
    except Exception:
        rf_proba, auc_rf = None, np.nan
    try:
        xgb_proba = cross_val_predict(xgb_clf, X_profile, y_bin, cv=cv_strat,
                                      method='predict_proba', n_jobs=-1)[:,1]
        auc_xgb = roc_auc_score(y_bin, xgb_proba)
    except Exception:
        auc_xgb = np.nan

    # Bootstrap 95% CI on the OOF predictions
    boot_aucs = []
    if rf_proba is not None:
        for _ in range(N_BOOT):
            idx = rng.integers(0, len(y_bin), len(y_bin))
            if len(np.unique(y_bin[idx])) < 2: continue
            try:
                boot_aucs.append(roc_auc_score(y_bin[idx], rf_proba[idx]))
            except Exception:
                pass
    if boot_aucs:
        auc_lo, auc_hi = np.percentile(boot_aucs, [2.5, 97.5])
    else:
        auc_lo, auc_hi = np.nan, np.nan

    print(f"  {task_name}:")
    print(f"    RF AUC = {auc_rf:.3f}  [95% CI: {auc_lo:.3f}-{auc_hi:.3f}]  "
          f"(boot n={len(boot_aucs)})")
    print(f"    XGB AUC = {auc_xgb:.3f}")
    m1_rows.append({
        'Task': task_name, 'N_pos': n_pos, 'N_neg': n_neg,
        'AUC_RF': auc_rf, 'AUC_RF_CI_low': auc_lo, 'AUC_RF_CI_high': auc_hi,
        'AUC_XGB': auc_xgb, 'N_bootstrap_valid': len(boot_aucs)
    })

m1_df = pd.DataFrame(m1_rows)
m1_df.to_csv('M1_classification_auc_bootstrap.csv', index=False)

# ============================================================================
# 9. CLUSTERING vs CHEMOTYPE — FISHER'S EXACT / PERMUTATION  [M2]
# ============================================================================
# Replaces asymptotic chi² with permutation test where expected counts are
# small (n=37 across multiple chemotypes -> many cells with E<5).
# Reference: Cochran (1954) on chi² validity bounds.
# ============================================================================
print("\n[9] CLUSTERING vs CHEMOTYPE with permutation chi² (replaces asymptotic)...")

from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from scipy.stats import chi2_contingency, fisher_exact

m2_rows = []
clustering_methods = [
    ('K-Means (k=3)', KMeans(n_clusters=3, n_init=20, random_state=RNG_SEED)),
    ('K-Means (k=4)', KMeans(n_clusters=4, n_init=20, random_state=RNG_SEED)),
    ('GMM (k=3)', GaussianMixture(n_components=3, n_init=5,
                                  covariance_type='full', random_state=RNG_SEED)),
]

for method_name, model in clustering_methods:
    try:
        labels = model.fit_predict(X_profile)
    except Exception as e:
        print(f"  {method_name}: FAILED -- {e}")
        continue

    cont = pd.crosstab(pd.Series(chemo_arr, name='chemo'),
                       pd.Series(labels, name='cluster'))
    chi2_stat, chi2_p, dof, expected = chi2_contingency(cont.values)
    n_low = int((expected < 5).sum())
    pct_low = n_low / expected.size * 100

    if cont.shape == (2, 2):
        odds, perm_p = fisher_exact(cont.values, alternative='two-sided')
        method_used = "Fisher's exact (2x2)"
    else:
        # 10000-iteration permutation chi² for >2x2 contingency
        n_perms = 10000
        rng_local = np.random.default_rng(RNG_SEED)
        null_chi2 = np.empty(n_perms)
        chemo_series = pd.Series(chemo_arr, name='chemo')
        for i in range(n_perms):
            shuffled = rng_local.permutation(labels)
            ct = pd.crosstab(chemo_series, shuffled)
            null_chi2[i] = chi2_contingency(ct.values)[0]
        perm_p = (null_chi2 >= chi2_stat - 1e-10).mean()
        method_used = f"Permutation chi² ({n_perms} perms)"

    print(f"  {method_name}: shape={cont.shape}, "
          f"low-E cells={n_low}/{expected.size} ({pct_low:.0f}%)")
    print(f"    asymptotic chi² stat = {chi2_stat:.2f}, asymptotic p = {chi2_p:.2e}")
    print(f"    {method_used}: p = {perm_p:.4f}")
    m2_rows.append({
        'Method': method_name,
        'Contingency_shape': str(cont.shape),
        'N_low_expected_cells': n_low,
        'Pct_low_expected': pct_low,
        'Chi2_statistic': chi2_stat,
        'Chi2_asymptotic_p': chi2_p,
        'Exact_or_permutation_method': method_used,
        'Exact_p_value': perm_p
    })

m2_df = pd.DataFrame(m2_rows)
m2_df.to_csv('M2_clustering_fisher_exact.csv', index=False)

# ============================================================================
# 10. LASSO BETA WITH BOOTSTRAP 95% CI  [M3]
# ============================================================================
# Per-target methoxy-position effects via LassoCV, with 1000-bootstrap CI on
# each beta coefficient. CI excluding 0 indicates robust position-specific
# effect; CI containing 0 indicates point estimate may be sample-specific.
# ============================================================================
print("\n[10] LASSO beta coefficients with bootstrap 95% CI...")

from sklearn.linear_model import LassoCV, Lasso

flavone_chemos = {'PMF', 'Meth_Flav', 'Standard_Flav'}
flavone_df = plant[plant['Chemotype'].isin(flavone_chemos)].copy()
METHOXY_COLS = ['M3','M5','M6','M7','M8','M3p','M4p','M5p']
print(f"  Flavone subset: {flavone_df['p_code'].nunique()} compounds")

m3_rows = []
for tgt in TARGETS_KEEP:
    sub = flavone_df[flavone_df['Target_Protein'] == tgt].copy()
    if len(sub) < 10:
        print(f"  {tgt}: SKIP (n={len(sub)})")
        continue
    avail = [c for c in METHOXY_COLS if c in sub.columns]
    X_lasso = sub[avail].values.astype(float)
    y_lasso = sub['MMGBSA dG Bind'].values
    n = len(y_lasso)

    keep_cols = X_lasso.std(axis=0) > 1e-8
    X_used = X_lasso[:, keep_cols]
    cols_used = [c for c, k in zip(avail, keep_cols) if k]
    if X_used.shape[1] == 0:
        print(f"  {tgt}: SKIP (no variable methoxy positions)")
        continue

    Xs_l = StandardScaler().fit_transform(X_used)
    try:
        lcv = LassoCV(cv=min(5, n), random_state=RNG_SEED, max_iter=20000).fit(Xs_l, y_lasso)
        alpha_opt = lcv.alpha_
        beta_point = lcv.coef_
    except Exception as e:
        print(f"  {tgt}: LassoCV failed -- {e}")
        continue

    boot_betas = np.zeros((N_BOOT, X_used.shape[1]))
    for b in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if X_used[idx].std(axis=0).min() < 1e-8:
            boot_betas[b] = np.nan
            continue
        try:
            Xb = StandardScaler().fit_transform(X_used[idx])
            yb = y_lasso[idx]
            lasso_b = Lasso(alpha=alpha_opt, max_iter=20000).fit(Xb, yb)
            boot_betas[b] = lasso_b.coef_
        except Exception:
            boot_betas[b] = np.nan

    print(f"  {tgt} (n={n}, alpha={alpha_opt:.4f}):")
    for j, col in enumerate(cols_used):
        valid = ~np.isnan(boot_betas[:, j])
        if valid.sum() < 100:
            continue
        lo, hi = np.percentile(boot_betas[valid, j], [2.5, 97.5])
        sig = (lo > 0) or (hi < 0)
        marker = '*' if sig else ''
        print(f"    {col:5s} beta = {beta_point[j]:+7.2f}   [{lo:+6.2f}, {hi:+6.2f}] {marker}")
        m3_rows.append({
            'Target': tgt, 'Position': col, 'Beta_point': beta_point[j],
            'CI_low': lo, 'CI_high': hi, 'CI_excludes_zero': bool(sig),
            'n_compounds': n, 'lasso_alpha': alpha_opt
        })

m3_df = pd.DataFrame(m3_rows)
m3_df.to_csv('M3_lasso_bootstrap_per_target.csv', index=False)

# ============================================================================
# 11. KRUSKAL-WALLIS H + PEARSON r vs NESTED R²  WITH BOOTSTRAP CI  [M4]
# ============================================================================
# Discloses imprecision of r at n=9 targets via 1000-replicate bootstrap CI.
# ============================================================================
print("\n[11] KRUSKAL-WALLIS H + PEARSON r with bootstrap CI...")

from scipy.stats import kruskal, pearsonr

per_target_h = {}
m4_rows = []
nested_dict = dict(zip(nested_df['Target'], nested_df['R2_Nested_CV']))

for tgt in TARGETS_KEEP:
    sub = plant[plant['Target_Protein'] == tgt]
    groups_by_chemo = [sub[sub['Chemotype'] == c]['MMGBSA dG Bind'].values
                       for c in sub['Chemotype'].unique()
                       if (sub['Chemotype'] == c).sum() >= 2]
    if len(groups_by_chemo) < 2:
        continue
    H, p = kruskal(*groups_by_chemo)
    per_target_h[tgt] = (H, p)
    print(f"  {tgt:6s}: H = {H:6.2f}  p = {p:.2e}  R² = {nested_dict.get(tgt, np.nan):+.3f}")
    m4_rows.append({
        'Target': tgt, 'KW_H': H, 'KW_p': p,
        'R2_Nested': nested_dict.get(tgt, np.nan)
    })

# Pearson r between H and R² across targets
H_arr = np.array([per_target_h[t][0] for t in TARGETS_KEEP if t in per_target_h])
R2_arr = np.array([nested_dict[t] for t in TARGETS_KEEP if t in per_target_h])
keep_finite = np.isfinite(H_arr) & np.isfinite(R2_arr)
H_arr, R2_arr = H_arr[keep_finite], R2_arr[keep_finite]

if len(H_arr) >= 4:
    r_point, r_p = pearsonr(H_arr, R2_arr)
    n = len(H_arr)
    r_boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if len(np.unique(idx)) < 3:
            continue
        try:
            r_boots.append(pearsonr(H_arr[idx], R2_arr[idx])[0])
        except Exception:
            pass
    r_lo, r_hi = np.percentile(r_boots, [2.5, 97.5])
    print(f"\n  Pearson r (H vs R²) at n={n}: r = {r_point:.3f}  asymptotic p = {r_p:.4f}")
    print(f"  Bootstrap 95% CI: [{r_lo:.3f}, {r_hi:.3f}]  ({len(r_boots)} valid)")
else:
    r_point, r_p, r_lo, r_hi, n = np.nan, np.nan, np.nan, np.nan, 0
    print(f"\n  WARNING: Only {len(H_arr)} valid target pairs; skipped correlation")

m4_summary = pd.DataFrame([{
    'Statistic': 'Pearson_r_H_vs_R2', 'n_targets': n,
    'r_point': r_point, 'p_asymptotic': r_p,
    'CI_low_bootstrap': r_lo, 'CI_high_bootstrap': r_hi
}])
pd.DataFrame(m4_rows).to_csv('M4_kruskal_per_target.csv', index=False)
m4_summary.to_csv('M4_pearson_r_bootstrap.csv', index=False)

# ============================================================================
# 12. BENJAMINI-HOCHBERG FDR CORRECTION on Y-RANDOMIZATION  [M5]
# ============================================================================
# Adjusts 9 simultaneous Y-rand p-values for multiple comparisons.
# Reports BH-FDR at α=0.05 (primary) and Bonferroni (conservative reference).
# ============================================================================
print("\n[12] BH-FDR correction on Y-randomization p-values...")

try:
    from scipy.stats import false_discovery_control
    HAS_FDR = True
except ImportError:
    HAS_FDR = False

if HAS_FDR:
    p_uncorr = yr_df['p_value'].values
    p_bh = false_discovery_control(p_uncorr, method='bh')
    p_bonf = np.minimum(p_uncorr * len(p_uncorr), 1.0)

    yr_df['p_BH_corrected'] = p_bh
    yr_df['p_Bonferroni'] = p_bonf
    yr_df['Significant_uncorrected'] = p_uncorr < 0.05
    yr_df['Significant_BH'] = p_bh < 0.05
    yr_df['Significant_Bonferroni'] = p_bonf < 0.05

    print(f"\n  {'Target':<8} {'p_raw':>9} {'p_BH':>9} {'p_Bonf':>9} {'Sig_BH':>8}")
    for _, row in yr_df.iterrows():
        print(f"  {row['Target']:<8} {row['p_value']:>9.4f} {row['p_BH_corrected']:>9.4f} "
              f"{row['p_Bonferroni']:>9.4f} {'YES' if row['Significant_BH'] else 'no':>8}")
    n_uncorr_sig = int(yr_df['Significant_uncorrected'].sum())
    n_bh_sig = int(yr_df['Significant_BH'].sum())
    n_bonf_sig = int(yr_df['Significant_Bonferroni'].sum())
    print(f"\n  Significant: {n_uncorr_sig}/9 (raw), {n_bh_sig}/9 (BH-FDR), "
          f"{n_bonf_sig}/9 (Bonferroni)")

    yr_df.to_csv('Y_RANDOMIZATION_RESULTS.csv', index=False)
    yr_df.to_csv('M5_yrandom_BH_corrected.csv', index=False)
else:
    # Manual BH if scipy < 1.11 doesn't have false_discovery_control
    print("  Manual BH implementation (scipy < 1.11)...")
    p_uncorr = yr_df['p_value'].values
    n = len(p_uncorr)
    sort_idx = np.argsort(p_uncorr)
    sorted_p = p_uncorr[sort_idx]
    bh_sorted = np.minimum.accumulate((sorted_p * n / np.arange(1, n+1))[::-1])[::-1]
    p_bh = np.empty_like(p_uncorr)
    p_bh[sort_idx] = bh_sorted
    p_bonf = np.minimum(p_uncorr * n, 1.0)
    yr_df['p_BH_corrected'] = p_bh
    yr_df['p_Bonferroni'] = p_bonf
    yr_df['Significant_BH'] = p_bh < 0.05
    yr_df['Significant_Bonferroni'] = p_bonf < 0.05
    yr_df.to_csv('Y_RANDOMIZATION_RESULTS.csv', index=False)
    yr_df.to_csv('M5_yrandom_BH_corrected.csv', index=False)

# ============================================================================
# 13. CONFORMAL COVERAGE GAP — FINITE-SAMPLE DISCLOSURE  [M6]
# ============================================================================
# Per-target coverage gap quantified against Lei et al. (2018) theoretical
# bound 1/(n_calib + 1). Reframed as honest finite-sample disclosure.
# ============================================================================
print("\n[13] Conformal coverage gap analysis...")

cov_per_target = conformal_df.groupby('Target')['Within_CI'].agg(['mean','count']).copy()
cov_per_target['Coverage_gap'] = 0.95 - cov_per_target['mean']
cov_per_target['Expected_max_deviation'] = 1 / (cov_per_target['count'] + 1)
cov_per_target['Within_theoretical_bound'] = (
    cov_per_target['Coverage_gap'].abs() <= cov_per_target['Expected_max_deviation']
)

print(f"\n  {'Target':<8} {'n_test':>7} {'cover':>8} {'gap':>8} {'bound':>8} {'within':>8}")
for tgt, row in cov_per_target.iterrows():
    print(f"  {tgt:<8} {int(row['count']):>7} {row['mean']:>8.3f} "
          f"{row['Coverage_gap']:>8.3f} {row['Expected_max_deviation']:>8.3f} "
          f"{'YES' if row['Within_theoretical_bound'] else 'no':>8}")
print(f"\n  Mean coverage: {cov_per_target['mean'].mean():.3f}  "
      f"({int(cov_per_target['Within_theoretical_bound'].sum())}/{len(cov_per_target)} "
      f"within theoretical bound)")
cov_per_target.to_csv('M6_conformal_coverage_summary.csv')

# ============================================================================
# Updated active-learning eligibility uses BH-corrected p-values
# ============================================================================
if 'p_BH_corrected' in yr_df.columns:
    bh_eligible = set(yr_df.loc[(yr_df['R2_real'] > 0.40) &
                                (yr_df['p_BH_corrected'] < 0.01), 'Target'])
    print(f"\n  Active-learning eligible targets (R² > 0.40 AND BH p < 0.01): "
          f"{sorted(bh_eligible) if bh_eligible else 'none'}")

print("\n" + "="*75)
print("PEER-REVIEW STATISTICAL ADD-ONS COMPLETE [M1-M6]")
print("="*75)

# ============================================================================
# SUMMARY OF ANALYSIS RUN
# ============================================================================
print("\n" + "="*75)
print("ANALYSIS COMPLETE — generating figures now")
print("="*75)
print("\nNESTED CV RESULTS:")
print(nested_df.round(3).to_string(index=False))

print("\nCONFORMAL COVERAGE SUMMARY:")
print(conf_perf.to_string(index=False))
print(f"  -> Overall: {conformal_df['Within_CI'].mean()*100:.1f}% (target: 95%)")

print("\nY-RANDOMIZATION:")
print(f"  -> Significant on {yr_df['Significant'].sum()}/9 targets")

# Compute GPR performance summary from predictions
gpr_perf = []
for tgt in gpr_df['Target'].unique():
    d = gpr_df[gpr_df['Target']==tgt]
    gpr_perf.append({
        'Target': tgt,
        'GPR_R2': r2_score(d['y_true'], d['y_pred']),
        'GPR_MAE': mean_absolute_error(d['y_true'], d['y_pred']),
        'Mean_sigma': d['sigma'].mean(),
        'Coverage_95pct': d['Within_CI'].mean()
    })
gpr_perf_df = pd.DataFrame(gpr_perf).round(3)

# Active learning summary (capped to prevent plot distortion)
al_df['R2_capped'] = al_df['R2_on_remaining'].clip(-5, 1)
al_summary = al_df.groupby(['Target','Strategy','N_acquired_pct'])['R2_capped'].agg(
    ['median','mean','std']).reset_index()

# ============================================================================
# FIGURE A — Y-RANDOMIZATION
# ============================================================================
print("\n[FIG A] Y-randomization…")
fig, ax = plt.subplots(figsize=(10, 6))
yr_sorted = yr_df.sort_values('R2_real', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(yr_sorted))

ax.barh(y_pos, yr_sorted['Null_95th'], color=COLORS['null'], alpha=0.7,
        edgecolor='black', linewidth=0.4, label='Null 95th %ile', height=0.7)

for i, row in yr_sorted.iterrows():
    sig = row['Significant']
    ax.scatter(row['R2_real'], i,
                color=COLORS['real'] if sig else COLORS['baseline'],
                s=140 if sig else 80,
                edgecolor='black', linewidth=0.8,
                zorder=5, alpha=1.0 if sig else 0.6)
    p_txt = f"p={row['p_value']:.3f}" if row['p_value'] >= 0.001 else "p<0.001"
    if sig: p_txt = r"$\bf{" + p_txt + "}$"
    ax.text(row['R2_real'] + 0.03, i, p_txt, fontsize=9, va='center',
            color='black' if sig else COLORS['neutral'])

ax.set_yticks(y_pos); ax.set_yticklabels(yr_sorted['Target'], fontsize=11)
ax.axvline(0, color='black', linewidth=0.5, alpha=0.5)
ax.set_xlabel('R² (LOCO-CV)', fontsize=11)
ax.set_title('Y-randomization: real model vs 500-permutation null distribution',
              fontsize=12, weight='bold', loc='left', pad=15)
handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['real'],
            markersize=12, markeredgecolor='black', label='Real R² (significant)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['baseline'],
            markersize=9, markeredgecolor='black', alpha=0.6,
            label='Real R² (not significant)'),
    mpatches.Patch(color=COLORS['null'], alpha=0.7, label='Null 95th %ile'),
]
ax.legend(handles=handles, loc='lower right', fontsize=9)
ax.set_xlim(-0.3, 0.95)
plt.tight_layout()
plt.savefig(OUT_FIG/'FigV5A_y_randomization.png')
plt.savefig(OUT_FIG/'FigV5A_y_randomization.svg')
plt.close()

# ============================================================================
# FIGURE B — APPLICABILITY DOMAIN
# ============================================================================
print("[FIG B] Applicability domain…")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

ax = axes[0]
chemo_order = [c for c in ['PMF','Meth_Flav','Gingeroid','Glycoside','Standard_Flav']
               if c in ad_df['Chemotype'].unique()]
palette = [COLORS[c] for c in chemo_order]

sns.violinplot(data=ad_df, x='Chemotype', y='Mean_Distance_to_KNN',
               order=chemo_order, palette=palette, ax=ax,
               inner='point', linewidth=1.0, saturation=0.85)
ax.axhline(ad_df['AD_Threshold'].iloc[0], color=COLORS['out_domain'],
           linestyle='--', linewidth=1.5, alpha=0.7, label='Williams AD threshold')
ax.set_xlabel(''); ax.set_ylabel('Mean distance to 5 nearest neighbors', fontsize=10)
ax.set_title('A | Chemical-space coverage by chemotype',
             fontsize=11, weight='bold', loc='left', pad=12)
ax.legend(fontsize=9, loc='upper left')
ax.tick_params(axis='x', rotation=15)

ax = axes[1]
agg = ad_df.groupby('Chemotype').agg(
    mean_dist=('Mean_Distance_to_KNN','mean'),
    sd_dist=('Mean_Distance_to_KNN','std'),
    n=('Mean_Distance_to_KNN','count')
).reindex(chemo_order)
colors_ch = [COLORS[c] for c in agg.index]
ax.bar(agg.index, agg['mean_dist'], yerr=agg['sd_dist'],
       color=colors_ch, edgecolor='black', linewidth=0.6, alpha=0.85,
       capsize=4, error_kw={'elinewidth':0.8, 'ecolor':COLORS['neutral']})
for i, (c, row) in enumerate(agg.iterrows()):
    ax.text(i, row['mean_dist'] + row['sd_dist'] + 0.2,
            f"n={int(row['n'])}", ha='center', fontsize=9)
ax.set_ylabel('Mean distance to 5 nearest neighbors', fontsize=10)
ax.set_title('B | Glycosides lie at chemical-space periphery',
             fontsize=11, weight='bold', loc='left', pad=12)
ax.tick_params(axis='x', rotation=15)

plt.tight_layout()
plt.savefig(OUT_FIG/'FigV5B_applicability_domain.png')
plt.savefig(OUT_FIG/'FigV5B_applicability_domain.svg')
plt.close()

# ============================================================================
# FIGURE C — GPR UNCERTAINTY
# ============================================================================
print("[FIG C] GPR uncertainty…")
fig, axes = plt.subplots(3, 3, figsize=(13, 11))
axes = axes.flatten()

for i, tgt in enumerate(TARGETS_KEEP):
    ax = axes[i]
    d = gpr_df[gpr_df['Target']==tgt]
    if len(d) == 0: continue
    for chem in ['PMF','Meth_Flav','Gingeroid','Glycoside','Standard_Flav']:
        dc = d[d['Chemotype']==chem]
        if len(dc) == 0: continue
        ax.errorbar(dc['y_true'], dc['y_pred'], yerr=1.96*dc['sigma'],
                    fmt='o', color=COLORS.get(chem,'gray'),
                    alpha=0.7, markersize=5, elinewidth=0.6, capsize=2,
                    markeredgecolor='black', markeredgewidth=0.3,
                    label=chem if i==0 else None)
    lo = min(d['y_true'].min(), d['y_pred'].min())
    hi = max(d['y_true'].max(), d['y_pred'].max())
    ax.plot([lo, hi], [lo, hi], 'k--', alpha=0.4, linewidth=0.7)

    perf = gpr_perf_df[gpr_perf_df['Target']==tgt].iloc[0]
    ax.text(0.03, 0.97,
            f"R² = {perf['GPR_R2']:+.2f}\nMAE = {perf['GPR_MAE']:.1f}\n"
            f"95% CI cov = {perf['Coverage_95pct']*100:.0f}%",
            transform=ax.transAxes, fontsize=8.5, va='top',
            bbox=dict(facecolor='white', edgecolor=COLORS['neutral'],
                      boxstyle='round,pad=0.3', alpha=0.9))
    ax.set_title(tgt, fontsize=11, weight='bold', loc='left', pad=6)
    if i >= 6: ax.set_xlabel('Observed ΔG (kcal/mol)', fontsize=9)
    if i % 3 == 0: ax.set_ylabel('Predicted ΔG (kcal/mol)', fontsize=9)

handles, labels = [], []
for chem in ['PMF','Meth_Flav','Gingeroid','Glycoside','Standard_Flav']:
    handles.append(Line2D([0],[0], marker='o', color='w',
                           markerfacecolor=COLORS[chem], markersize=9,
                           markeredgecolor='black'))
    labels.append(chem)
fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=10,
           bbox_to_anchor=(0.5, -0.01), frameon=False)

fig.suptitle(f"GPR predictions with 95% credible intervals "
             f"(mean coverage: {gpr_df['Within_CI'].mean()*100:.1f}%)",
             fontsize=13, weight='bold', y=0.995)
plt.tight_layout(rect=[0, 0.02, 1, 0.99])
plt.savefig(OUT_FIG/'FigV5C_gpr_uncertainty.png')
plt.savefig(OUT_FIG/'FigV5C_gpr_uncertainty.svg')
plt.close()

# ============================================================================
# FIGURE D — LEARNING CURVE
# ============================================================================
print("[FIG D] Learning curve…")
fig, axes = plt.subplots(3, 3, figsize=(13, 11))
axes = axes.flatten()
for i, tgt in enumerate(TARGETS_KEEP):
    ax = axes[i]
    d = lc_df[lc_df['Target']==tgt].sort_values('Train_Fraction')
    if len(d) == 0: continue
    ax.errorbar(d['Train_Fraction']*100, d['R2_mean'], yerr=d['R2_std'],
                fmt='-o', color=COLORS['model'], linewidth=2, markersize=7,
                capsize=4, ecolor=COLORS['uncertainty'], markeredgecolor='black',
                markeredgewidth=0.4)
    ax.axhline(0, color='gray', linewidth=0.4, alpha=0.5)
    ax.axvspan(50, 70, color=COLORS['in_domain'], alpha=0.08,
               label='Stable zone' if i==0 else None)
    ax.set_title(tgt, fontsize=11, weight='bold', loc='left', pad=6)
    if i >= 6: ax.set_xlabel('Training fraction (%)', fontsize=9)
    if i % 3 == 0: ax.set_ylabel('R² (± SD)', fontsize=9)

fig.suptitle('Learning curve — performance vs training size',
             fontsize=13, weight='bold', y=0.99)
plt.tight_layout()
plt.savefig(OUT_FIG/'FigV5D_learning_curve.png')
plt.savefig(OUT_FIG/'FigV5D_learning_curve.svg')
plt.close()

# ============================================================================
# FIGURE E — ACTIVE LEARNING
# ============================================================================
print("[FIG E] Active learning…")
fig, axes = plt.subplots(3, 3, figsize=(13, 11))
axes = axes.flatten()
for i, tgt in enumerate(TARGETS_KEEP):
    ax = axes[i]
    d = al_summary[al_summary['Target']==tgt]
    if len(d) == 0: continue
    for strat, col, label in [
        ('UncertaintySampling', COLORS['PMF'], 'GPR uncertainty'),
        ('Random', COLORS['baseline'], 'Random')]:
        s = d[d['Strategy']==strat].sort_values('N_acquired_pct')
        if len(s) == 0: continue
        ax.plot(s['N_acquired_pct'], s['median'], '-', linewidth=2.2,
                color=col, label=label if i==0 else None, alpha=0.9)
        low = (s['median'] - 0.5*s['std']).clip(lower=-1)
        high = (s['median'] + 0.5*s['std']).clip(upper=1)
        ax.fill_between(s['N_acquired_pct'], low, high, color=col, alpha=0.15)
    ax.axhline(0, color='gray', linewidth=0.4, alpha=0.5)
    ax.set_ylim(-1.05, 1.05)
    ax.set_title(tgt, fontsize=11, weight='bold', loc='left', pad=6)
    if i >= 6: ax.set_xlabel('% dataset acquired', fontsize=9)
    if i % 3 == 0: ax.set_ylabel('R² on remaining', fontsize=9)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=10,
           bbox_to_anchor=(0.5, -0.01), frameon=False)
fig.suptitle('Active learning: GPR uncertainty vs random acquisition',
             fontsize=12, weight='bold', y=0.99)
plt.tight_layout(rect=[0, 0.02, 1, 0.99])
plt.savefig(OUT_FIG/'FigV5E_active_learning.png')
plt.savefig(OUT_FIG/'FigV5E_active_learning.svg')
plt.close()

# ============================================================================
# FIGURE F — BASELINE BENCHMARKS
# ============================================================================
print("[FIG F] Baseline benchmarks…")
fig, ax = plt.subplots(figsize=(11, 6))
bench_sorted = baseline_df.sort_values('R2_Full_Features', ascending=False).reset_index(drop=True)
x = np.arange(len(bench_sorted))
width = 0.27
ax.bar(x - width, bench_sorted['R2_Physchem_Only'], width,
       color=COLORS['baseline'], edgecolor='black', linewidth=0.4,
       label='Physchem only', alpha=0.85)
ax.bar(x, bench_sorted['R2_Position_Only'], width,
       color=COLORS['Meth_Flav'], edgecolor='black', linewidth=0.4,
       label='Position only', alpha=0.85)
ax.bar(x + width, bench_sorted['R2_Full_Features'], width,
       color=COLORS['PMF'], edgecolor='black', linewidth=0.4,
       label='Full features', alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(bench_sorted['Target'], fontsize=10)
ax.set_ylabel('R² (LOCO-CV)', fontsize=11)
ax.axhline(0, color='black', linewidth=0.4)
ax.axhline(0.3, color=COLORS['baseline'], linestyle=':', alpha=0.7,
           label='Publishable threshold')
ax.legend(loc='lower left', fontsize=9, ncol=1)
ax.set_title('Baseline benchmarks — feature set contribution by target',
             fontsize=12, weight='bold', loc='left', pad=12)
plt.tight_layout()
plt.savefig(OUT_FIG/'FigV5F_baseline_benchmarks.png')
plt.savefig(OUT_FIG/'FigV5F_baseline_benchmarks.svg')
plt.close()

# ============================================================================
# FIGURE G — CONFORMAL COVERAGE  [NEW IN V5]
# ============================================================================
print("[FIG G] Conformal coverage…")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel A: coverage bar chart
ax = axes[0]
conf_sorted = conf_perf.sort_values('Coverage', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(conf_sorted))
bars = ax.barh(y_pos, conf_sorted['Coverage']*100,
                color=[COLORS['conformal'] if c >= 0.85 else COLORS['baseline']
                       for c in conf_sorted['Coverage']],
                edgecolor='black', linewidth=0.4, alpha=0.85)
ax.axvline(95, color=COLORS['real'], linestyle='--', linewidth=1.5,
           alpha=0.7, label='Target coverage (95%)')
ax.set_yticks(y_pos); ax.set_yticklabels(conf_sorted['Target'])
ax.set_xlabel('Empirical 95% CI coverage (%)', fontsize=10)
ax.set_title('A | Conformal prediction coverage',
             fontsize=11, weight='bold', loc='left', pad=12)
ax.legend(fontsize=9, loc='lower right')
for i, row in conf_sorted.iterrows():
    ax.text(row['Coverage']*100 + 1, i, f"{row['Coverage']*100:.0f}%",
            fontsize=9, va='center')

# Panel B: interval width distribution
ax = axes[1]
sns.boxplot(data=conformal_df, x='Target', y='Interval_Width',
            order=[t for t in TARGETS_KEEP if t in conformal_df['Target'].values],
            palette='Purples', ax=ax, linewidth=0.8)
ax.set_xlabel(''); ax.set_ylabel('Conformal interval width (kcal/mol)', fontsize=10)
ax.set_title('B | Prediction interval tightness per target',
             fontsize=11, weight='bold', loc='left', pad=12)
ax.tick_params(axis='x', rotation=30)

plt.tight_layout()
plt.savefig(OUT_FIG/'FigV5G_conformal_coverage.png')
plt.savefig(OUT_FIG/'FigV5G_conformal_coverage.svg')
plt.close()

# ============================================================================
# FIGURE H — NESTED CV vs SINGLE-FOLD XGB  [NEW IN V5]
# ============================================================================
print("[FIG H] Nested CV comparison…")
# Compare nested R² to Y-rand real R² (which is single-fold XGBoost LOCO)
merged = nested_df[['Target','R2_Nested_CV']].merge(
    yr_df[['Target','R2_real']].rename(columns={'R2_real':'R2_SingleFold_XGB'}),
    on='Target'
)
merged['R2_Diff'] = merged['R2_SingleFold_XGB'] - merged['R2_Nested_CV']

fig, ax = plt.subplots(figsize=(10, 6))
merged_sorted = merged.sort_values('R2_Nested_CV', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(merged_sorted))
ax.barh(y_pos - 0.2, merged_sorted['R2_SingleFold_XGB'], 0.4,
        color=COLORS['baseline'], edgecolor='black', linewidth=0.4,
        label='Single-fold XGBoost (fixed hp)', alpha=0.85)
ax.barh(y_pos + 0.2, merged_sorted['R2_Nested_CV'], 0.4,
        color=COLORS['nested'], edgecolor='black', linewidth=0.4,
        label='Nested CV (inner-CV tuned hp)', alpha=0.85)
ax.set_yticks(y_pos); ax.set_yticklabels(merged_sorted['Target'])
ax.axvline(0, color='black', linewidth=0.4)
ax.axvline(0.3, color=COLORS['baseline'], linestyle=':', alpha=0.6,
           label='Publishable threshold')
ax.set_xlabel('R² (LOCO-CV)', fontsize=11)
ax.set_title('Nested CV vs single-fold XGBoost — honest hyperparameter validation',
             fontsize=12, weight='bold', loc='left', pad=12)
ax.legend(fontsize=9, loc='lower right')
plt.tight_layout()
plt.savefig(OUT_FIG/'FigV5H_nested_cv_comparison.png')
plt.savefig(OUT_FIG/'FigV5H_nested_cv_comparison.svg')
plt.close()

# ============================================================================
# HERO FIGURE — 4-panel validation summary
# ============================================================================
print("[FIG HERO] Validation summary…")
fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28,
                      left=0.07, right=0.97, top=0.93, bottom=0.08)

# Panel A — Y-rand
ax = fig.add_subplot(gs[0, 0])
yr_sorted = yr_df.sort_values('R2_real', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(yr_sorted))
ax.barh(y_pos, yr_sorted['Null_95th'], color=COLORS['null'], alpha=0.7,
        edgecolor='black', linewidth=0.3)
for i, row in yr_sorted.iterrows():
    sig = row['Significant']
    ax.scatter(row['R2_real'], i,
               color=COLORS['real'] if sig else COLORS['baseline'],
               s=110 if sig else 55, zorder=5, edgecolor='black', linewidth=0.5,
               alpha=1.0 if sig else 0.6)
    if sig and row['p_value'] < 0.003:
        ax.text(row['R2_real']+0.04, i, '**', fontsize=14, va='center', weight='bold')
    elif sig:
        ax.text(row['R2_real']+0.04, i, '*', fontsize=13, va='center', weight='bold')
ax.set_yticks(y_pos); ax.set_yticklabels(yr_sorted['Target'])
ax.axvline(0, color='black', linewidth=0.4, alpha=0.5)
ax.set_xlabel('R²', fontsize=10)
ax.set_title(f"A | Y-randomization:  {yr_df['Significant'].sum()}/9 targets reject null",
             fontsize=12, weight='bold', loc='left', pad=12)
legh = [mpatches.Patch(color=COLORS['null'], alpha=0.7, label='Null 95th %ile'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['real'],
                markersize=10, markeredgecolor='black', label='Real R²')]
ax.legend(handles=legh, loc='lower right', fontsize=9)

# Panel B — Nested CV + conformal coverage (dual axis)
ax = fig.add_subplot(gs[0, 1])
nested_sorted = nested_df.sort_values('R2_Nested_CV', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(nested_sorted))
bars = ax.barh(y_pos, nested_sorted['R2_Nested_CV'],
                color=[COLORS['nested'] if r > 0.3 else COLORS['baseline']
                       for r in nested_sorted['R2_Nested_CV']],
                edgecolor='black', linewidth=0.5, alpha=0.85)
ax.set_yticks(y_pos); ax.set_yticklabels(nested_sorted['Target'])
ax.axvline(0, color='black', linewidth=0.4)
ax.axvline(0.3, color=COLORS['baseline'], linestyle=':', alpha=0.6)
ax.set_xlabel('Nested CV R²', fontsize=10)
# Annotate with conformal coverage
for i, row in nested_sorted.iterrows():
    tgt = row['Target']
    cov_row = conf_perf[conf_perf['Target']==tgt]
    if len(cov_row) > 0:
        cov = cov_row.iloc[0]['Coverage'] * 100
        ax.text(row['R2_Nested_CV'] + 0.03 if row['R2_Nested_CV'] > 0 else -0.05, i,
                f"cf={cov:.0f}%", fontsize=8.5, va='center',
                ha='left' if row['R2_Nested_CV'] > 0 else 'right',
                color=COLORS['conformal'], weight='bold')
ax.set_title(f"B | Nested CV R² (conformal coverage annotated)",
             fontsize=12, weight='bold', loc='left', pad=12)

# Panel C — Baselines
ax = fig.add_subplot(gs[1, 0])
bench_sorted = baseline_df.sort_values('R2_Full_Features', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(bench_sorted))
ax.barh(y_pos - 0.27, bench_sorted['R2_Physchem_Only'], 0.27,
        color=COLORS['baseline'], alpha=0.75, edgecolor='black', linewidth=0.3,
        label='Physchem only')
ax.barh(y_pos, bench_sorted['R2_Position_Only'], 0.27,
        color=COLORS['Meth_Flav'], alpha=0.75, edgecolor='black', linewidth=0.3,
        label='Position only')
ax.barh(y_pos + 0.27, bench_sorted['R2_Full_Features'], 0.27,
        color=COLORS['PMF'], alpha=0.85, edgecolor='black', linewidth=0.3,
        label='Full features')
ax.set_yticks(y_pos); ax.set_yticklabels(bench_sorted['Target'])
ax.axvline(0, color='black', linewidth=0.4, alpha=0.5)
ax.axvline(0.3, color=COLORS['baseline'], linestyle=':', alpha=0.6)
ax.set_xlabel('R² (LOCO-CV)', fontsize=10)
ax.set_title('C | Baseline benchmarks',
             fontsize=12, weight='bold', loc='left', pad=12)
ax.legend(loc='lower right', fontsize=8)

# Panel D — AL on CPT1A + VDAC1
ax = fig.add_subplot(gs[1, 1])
for tgt, ls in [('CPT1A','-'), ('VDAC1','--')]:
    d = al_summary[al_summary['Target']==tgt]
    if len(d) == 0: continue
    u = d[d['Strategy']=='UncertaintySampling'].sort_values('N_acquired_pct')
    r = d[d['Strategy']=='Random'].sort_values('N_acquired_pct')
    ax.plot(u['N_acquired_pct'], u['median'], ls, linewidth=2.5,
            color=COLORS['PMF'], alpha=0.9,
            label=f'{tgt} — GPR uncertainty')
    ax.plot(r['N_acquired_pct'], r['median'], ls, linewidth=1.8,
            color=COLORS['baseline'], alpha=0.7,
            label=f'{tgt} — Random')
ax.axhline(0, color='gray', linewidth=0.4, alpha=0.5)
ax.set_xlabel('% dataset acquired', fontsize=10)
ax.set_ylabel('R² on remaining compounds', fontsize=10)
ax.set_title('D | Active learning at PMF-favoured targets',
             fontsize=12, weight='bold', loc='left', pad=12)
ax.legend(loc='lower right', fontsize=8)
ax.set_ylim(-0.6, 1.0)

fig.suptitle('ML validation stack — Pharmaceutics (IF 5+) ready',
             fontsize=15, weight='bold', y=0.99)
plt.savefig(OUT_FIG/'FigV5_hero_validation.png')
plt.savefig(OUT_FIG/'FigV5_hero_validation.svg')
plt.close()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "="*75)
print("V5 PIPELINE COMPLETE — analysis + figures in one run")
print("="*75)
print(f"\nKey numbers:")
print(f"  Nested CV R² > 0.3 on: {(nested_df['R2_Nested_CV']>0.3).sum()}/9 targets")
print(f"  Y-randomization significant: {yr_df['Significant'].sum()}/9 targets")
print(f"  GPR mean coverage: {gpr_df['Within_CI'].mean()*100:.1f}% (ideal: 95%)")
print(f"  Conformal mean coverage: {conformal_df['Within_CI'].mean()*100:.1f}% (target: 95%)")

print(f"\nOutputs in current directory:")
for f in ['NESTED_CV_RESULTS.csv', 'CONFORMAL_PREDICTIONS.csv',
          'HYPERPARAMETER_SCAN.csv', 'Y_RANDOMIZATION_RESULTS.csv',
          'APPLICABILITY_DOMAIN.csv', 'GPR_UNCERTAINTY_PREDICTIONS.csv',
          'BASELINE_BENCHMARKS.csv', 'LEARNING_CURVE.csv',
          'ACTIVE_LEARNING_SIMULATION.csv', 'ML_VALIDATION_COMPLETE.xlsx']:
    if os.path.exists(f):
        print(f"  [OK] {f}")

print(f"\nFigures in ./{OUT_FIG}/:")
for f in sorted(os.listdir(OUT_FIG)):
    print(f"  [OK] {f}")

print("\nReady for manuscript audit. Upload the 3 new CSVs + hero figure:")
print("  1. NESTED_CV_RESULTS.csv")
print("  2. CONFORMAL_PREDICTIONS.csv")
print("  3. HYPERPARAMETER_SCAN.csv")
print("  4. figures_v5/FigV5_hero_validation.png")
