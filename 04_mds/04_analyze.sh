#!/bin/bash
# ============================================================
# 04_analyze.sh - Trajectory analysis pipeline
# ============================================================
# Takes a 100 ns production trajectory and produces:
#   - Preprocessed trajectory (PBC + center + fit)
#   - RMSD (backbone, ligand)
#   - RMSF (per residue)
#   - Radius of gyration
#   - H-bond counts (protein-ligand)
#   - MM-GBSA ΔG_bind on last 50 ns
#   - Per-residue contact frequencies (binary residue-level)
#   - Cluster medoid frame (representative structure)
#   - PLIP 2D ligand interaction diagram
#   - Summary report
#
# USAGE: ./04_analyze.sh <sysname>
# EXAMPLE: ./04_analyze.sh test_bcl2pc
#
# REQUIREMENTS:
#   - 03_production.sh must have completed for this sysname
#   - Conda env with: gromacs (CUDA), gmx_MMPBSA, MDAnalysis, plip
#   - Walltime: ~30-60 min on CPU
#
# OUTPUT: All analysis files + summary.txt
# ============================================================

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <sysname>"
    exit 1
fi

SYSNAME="$1"
WORKDIR="/scratch/amans/mds/runs/$SYSNAME"

if [ ! -d "$WORKDIR" ]; then
    echo "ERROR: Directory not found: $WORKDIR"
    exit 1
fi

cd "$WORKDIR"

for f in production.xtc production.tpr topol.top index.ndx; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing $f - run 03_production.sh first"
        exit 1
    fi
done

source /home/amans/mambaforge/etc/profile.d/conda.sh
conda activate md

echo "============================================================"
echo "04_analyze.sh - $SYSNAME"
echo "Started: $(date)"
echo "============================================================"

# Clean up any partial outputs
rm -f production_nojump.xtc production_centered.xtc production_fit.xtc \
      rmsd_*.xvg rmsf_*.xvg gyrate.xvg hbonds_*.xvg \
      contact_frequencies.csv medoid_frame.pdb test_frame.pdb \
      summary.txt
rm -rf plip_output _GMXMMPBSA_* COM_traj_* receptor_traj_* ligand_traj_* \
       FINAL_RESULTS_MMPBSA.dat *.h5

# ============================================================
# STAGE 1: TRAJECTORY PREPROCESSING (PBC + center + fit)
# ============================================================
echo ""
echo "[STAGE 1] Trajectory preprocessing (trjconv 3-step)"

# Step 1: nojump on full System (fix periodic-boundary jumps)
echo "System" | gmx trjconv -s production.tpr -f production.xtc \
    -o production_nojump.xtc -pbc nojump > trjconv1.log 2>&1
[ ! -f production_nojump.xtc ] && { echo "ERROR: trjconv1 failed"; tail -10 trjconv1.log; exit 1; }

# Step 2: center Protein in box, output System
echo -e "Protein\nSystem" | gmx trjconv -s production.tpr -f production_nojump.xtc \
    -o production_centered.xtc -pbc mol -center > trjconv2.log 2>&1
[ ! -f production_centered.xtc ] && { echo "ERROR: trjconv2 failed"; tail -10 trjconv2.log; exit 1; }

# Step 3: fit on Backbone (rotation + translation), output System
echo -e "Backbone\nSystem" | gmx trjconv -s production.tpr -f production_centered.xtc \
    -o production_fit.xtc -fit rot+trans > trjconv3.log 2>&1
[ ! -f production_fit.xtc ] && { echo "ERROR: trjconv3 failed"; tail -10 trjconv3.log; exit 1; }

# Cleanup intermediate files (keep production.xtc + production_fit.xtc)
rm -f production_nojump.xtc production_centered.xtc

echo "  Trajectory preprocessing complete"
ls -la production_fit.xtc


# ============================================================
# STAGE 2: RMSD, RMSF, Rg, H-BONDS (via MDAnalysis)
# ============================================================
# MDAnalysis-based to avoid two known gmx bugs:
#   - gmx rms/rmsf/gyrate: PBC wrap artifacts on some frames
#   - gmx hbond: doesn't auto-detect donors/acceptors for non-standard ligands
echo ""
echo "[STAGE 2] RMSD/RMSF/Rg/H-bonds via MDAnalysis"

python3 << 'PYEOF'
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
import numpy as np

u = mda.Universe("production.tpr", "production_fit.xtc")
n_frames = len(u.trajectory)
print(f"  Loaded {n_frames} frames")

# RMSD - backbone (vs frame 0)
backbone = u.select_atoms("backbone")
ref_pos = backbone.positions.copy()
ligand = u.select_atoms("resname LIG")
ref_lig = ligand.positions.copy()

times, rmsds_bb, rmsds_lig = [], [], []
for ts in u.trajectory:
    times.append(ts.time / 1000.0)  # ps -> ns
    rmsds_bb.append(rms.rmsd(backbone.positions, ref_pos, superposition=False) / 10.0)
    rmsds_lig.append(rms.rmsd(ligand.positions, ref_lig, superposition=False) / 10.0)

with open("rmsd_backbone.xvg", "w") as f:
    f.write("# RMSD backbone vs frame 0 (nm)\n")
    f.write("# time(ns)  rmsd(nm)\n")
    for t, r in zip(times, rmsds_bb):
        f.write(f"{t:.3f}  {r:.4f}\n")

with open("rmsd_ligand.xvg", "w") as f:
    f.write("# RMSD ligand (nm)\n")
    f.write("# time(ns)  rmsd(nm)\n")
    for t, r in zip(times, rmsds_lig):
        f.write(f"{t:.3f}  {r:.4f}\n")

print(f"  RMSD backbone: mean={np.mean(rmsds_bb):.3f}, final={rmsds_bb[-1]:.3f} nm")
print(f"  RMSD ligand:   mean={np.mean(rmsds_lig):.3f}, final={rmsds_lig[-1]:.3f} nm")

# RMSF per Cα
calpha = u.select_atoms("name CA")
rmsf_calc = rms.RMSF(calpha).run()
with open("rmsf_calpha.xvg", "w") as f:
    f.write("# RMSF per residue (nm)\n")
    f.write("# resid  rmsf(nm)\n")
    for resid, rval in zip(calpha.resids, rmsf_calc.results.rmsf):
        f.write(f"{resid}  {rval/10.0:.4f}\n")
print(f"  RMSF: mean={np.mean(rmsf_calc.results.rmsf)/10.0:.3f} nm, max={np.max(rmsf_calc.results.rmsf)/10.0:.3f} nm")

# Radius of gyration
rg_values = []
for ts in u.trajectory:
    rg = u.select_atoms("protein").radius_of_gyration() / 10.0
    rg_values.append(rg)
with open("gyrate.xvg", "w") as f:
    f.write("# Radius of gyration (nm)\n")
    f.write("# time(ns)  Rg(nm)\n")
    for t, r in zip(times, rg_values):
        f.write(f"{t:.3f}  {r:.4f}\n")
print(f"  Rg: mean={np.mean(rg_values):.3f}, final={rg_values[-1]:.3f} nm")

# H-bonds: protein donors -> ligand acceptors
hba1 = HBA(
    universe=u,
    donors_sel="protein and (name N* or name O*) and not (name C* or name CA)",
    acceptors_sel="resname LIG and (name O* or name N*)",
    d_a_cutoff=3.5,
    d_h_a_angle_cutoff=120.0,
    update_selections=False
)
hba1.run()

# H-bonds: ligand donors -> protein acceptors
hba2 = HBA(
    universe=u,
    donors_sel="resname LIG and (name N* or name O*)",
    acceptors_sel="protein and (name O* or name N*)",
    d_a_cutoff=3.5,
    d_h_a_angle_cutoff=120.0,
    update_selections=False
)
hba2.run()

n_hbonds_per_frame = hba1.count_by_time() + hba2.count_by_time()
with open("hbonds_count.xvg", "w") as f:
    f.write("# H-bonds protein-ligand per frame\n")
    f.write("# time(ns)  n_hbonds\n")
    for t, n in zip(times, n_hbonds_per_frame):
        f.write(f"{t:.3f}  {n}\n")
print(f"  H-bonds: mean={np.mean(n_hbonds_per_frame):.2f}, max={int(max(n_hbonds_per_frame))}")
PYEOF

echo "[STAGE 2] Done"

# ============================================================
# STAGE 3: MM-GBSA (last 50 ns, frames 500-1000)
# ============================================================
echo ""
echo "[STAGE 3] MM-GBSA on last 50 ns"

cat > mmgbsa.in << 'MMGBSA'
&general
  startframe=500, endframe=1000, interval=2,
  forcefields="oldff/leaprc.ff14SB",
  PBRadii=4, temperature=310,
/
&gb
  igb=8, saltcon=0.150, surften=0.0072,
/
MMGBSA

# Find Protein and LIG group numbers in index.ndx
PROT_NUM=$(grep -n "^\[ Protein \]$" index.ndx | cut -d: -f1 | head -1)
LIG_NUM=$(grep -n "^\[ LIG \]$" index.ndx | cut -d: -f1 | head -1)
PROT_IDX=$(($(grep -c "^\[ " <(head -n $PROT_NUM index.ndx)) - 1))
LIG_IDX=$(($(grep -c "^\[ " <(head -n $LIG_NUM index.ndx)) - 1))
echo "  Protein group #: $PROT_IDX, LIG group #: $LIG_IDX"

gmx_MMPBSA -O -i mmgbsa.in -cs production.tpr -ci index.ndx \
    -cg $PROT_IDX $LIG_IDX -ct production_fit.xtc -cp topol.top \
    -nogui > mmgbsa_run.log 2>&1

if [ -f FINAL_RESULTS_MMPBSA.dat ]; then
    DG_BIND=$(awk '/Delta \(Complex/,/^$/' FINAL_RESULTS_MMPBSA.dat | grep "ΔTOTAL" | awk '{print $2}')
    echo "  MM-GBSA ΔG_bind: $DG_BIND kcal/mol"
else
    echo "  ⚠️  MM-GBSA failed - check mmgbsa_run.log"
    tail -20 mmgbsa_run.log
fi

# ============================================================
# STAGE 4: PER-RESIDUE CONTACT FREQUENCIES + MEDOID FRAME
# ============================================================
# Binary residue-level counting (a residue is in contact or not)
# Avoids the atom-level multi-counting bug that produced >100% values
echo ""
echo "[STAGE 4] Contact analysis + medoid frame"

python3 << 'PYEOF'
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

u = mda.Universe("production.tpr", "production_fit.xtc")
ligand = u.select_atoms("resname LIG")
protein = u.select_atoms("protein")
n_frames = len(u.trajectory)
residues = list(protein.residues)
print(f"  Frames: {n_frames}, residues: {len(residues)}")

# Binary per-residue contact counting
residue_contact_count = {}
for ts in u.trajectory:
    for res in residues:
        d = distances.distance_array(ligand.positions, res.atoms.positions).min()
        if d < 4.0:
            key = f"{res.resname}{res.resid}"
            residue_contact_count[key] = residue_contact_count.get(key, 0) + 1

# Write CSV (only residues in contact >10% of time)
with open("contact_frequencies.csv", "w") as f:
    f.write("residue,frequency_pct\n")
    n_persistent = 0
    for res, count in sorted(residue_contact_count.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / n_frames
        if pct > 10.0:
            f.write(f"{res},{pct:.2f}\n")
            n_persistent += 1

print(f"  Residues in contact >10% of time: {n_persistent}")
print(f"  Top 5 contacts:")
top5 = sorted(residue_contact_count.items(), key=lambda x: -x[1])[:5]
for res, count in top5:
    pct = 100.0 * count / n_frames
    print(f"    {res}: {pct:.2f}%")

# Cluster medoid frame (representative structure for PLIP)
print("  Computing cluster medoid...")
backbone = u.select_atoms("backbone")
positions_list = []
frame_indices = []
for ts in u.trajectory[::10]:
    positions_list.append(backbone.positions.copy())
    frame_indices.append(ts.frame)

ref_arr = np.array(positions_list)
n_samples = len(ref_arr)
sum_rmsd = np.zeros(n_samples)
for i in range(n_samples):
    for j in range(n_samples):
        if i != j:
            diff = ref_arr[i] - ref_arr[j]
            sum_rmsd[i] += np.sqrt((diff**2).sum() / len(diff))

medoid_idx = int(np.argmin(sum_rmsd))
medoid_frame = frame_indices[medoid_idx]
print(f"  Medoid frame index: {medoid_frame}")

u.trajectory[medoid_frame]
all_atoms = u.select_atoms("protein or resname LIG")
all_atoms.write("medoid_frame.pdb")
print("  Saved medoid_frame.pdb")
PYEOF

echo "[STAGE 4] Done"

# ============================================================
# STAGE 5: PLIP 2D LIGAND INTERACTION DIAGRAM
# ============================================================
echo ""
echo "[STAGE 5] PLIP analysis on medoid frame"

if [ -f medoid_frame.pdb ]; then
    mkdir -p plip_output
    plip -f medoid_frame.pdb -o plip_output --txt --xml > plip_run.log 2>&1
    if ls plip_output/*report.txt &>/dev/null; then
        echo "  PLIP report generated"
        # Quick summary of interactions
        N_HBONDS=$(grep -c "^| " plip_output/*report.txt 2>/dev/null | head -1 || echo 0)
        echo "  PLIP found interactions in: $(ls plip_output/ | head -5)"
    else
        echo "  ⚠️  PLIP did not produce expected output"
        tail -10 plip_run.log
    fi
else
    echo "  ⚠️  No medoid frame to analyze"
fi

# ============================================================
# STAGE 6: SUMMARY REPORT
# ============================================================
echo ""
echo "[STAGE 6] Generating summary"

cat > summary.txt << SUMMARY
============================================================
ANALYSIS SUMMARY: $SYSNAME
Date: $(date)
============================================================

--- TRAJECTORY ---
Production: $(ls -la production.xtc | awk '{print $5}') bytes
Frames: $(awk 'END {print NR-2}' rmsd_backbone.xvg)

--- RMSD ---
Backbone (last 20 ns avg): $(awk 'NF==2 && \$1+0>=80 {sum+=\$2; n+=1} END {if (n>0) printf "%.3f nm", sum/n}' rmsd_backbone.xvg)
Backbone (final): $(awk 'NF==2 {last=\$2} END {printf "%.3f nm", last}' rmsd_backbone.xvg)
Ligand (last 20 ns avg): $(awk 'NF==2 && \$1+0>=80 {sum+=\$2; n+=1} END {if (n>0) printf "%.3f nm", sum/n}' rmsd_ligand.xvg)

--- RMSF ---
Mean Cα RMSF: $(awk 'NF==2 {sum+=\$2; n+=1} END {if (n>0) printf "%.3f nm", sum/n}' rmsf_calpha.xvg)
Max Cα RMSF: $(awk 'NF==2 {if (\$2>m) m=\$2} END {printf "%.3f nm", m}' rmsf_calpha.xvg)

--- RADIUS OF GYRATION ---
Mean: $(awk 'NF==2 {sum+=\$2; n+=1} END {if (n>0) printf "%.3f nm", sum/n}' gyrate.xvg)
Final: $(awk 'NF==2 {last=\$2} END {printf "%.3f nm", last}' gyrate.xvg)

--- H-BONDS ---
Mean per frame (50-100 ns): $(awk 'NF==2 && \$1+0>=50 {sum+=\$2; n+=1} END {if (n>0) printf "%.2f", sum/n}' hbonds_count.xvg)
Max: $(awk 'NF==2 {if (\$2>m) m=\$2} END {printf "%d", m}' hbonds_count.xvg)

--- MM-GBSA ΔG_bind (last 50 ns) ---
$(awk '/Delta \(Complex/,/^$/' FINAL_RESULTS_MMPBSA.dat 2>/dev/null | grep -E "ΔVDWAALS|ΔEEL|ΔEGB|ΔESURF|ΔTOTAL")

--- TOP 10 CONTACT RESIDUES (>10% of frames) ---
$(head -11 contact_frequencies.csv 2>/dev/null)

--- PLIP MEDOID FRAME INTERACTIONS ---
Frame: medoid_frame.pdb
Output: plip_output/

============================================================
SUMMARY

cat summary.txt

echo ""
echo "============================================================"
echo "04_analyze.sh COMPLETE: $SYSNAME"
echo "Finished: $(date)"
echo ""
echo "Files generated:"
echo "  ✓ production_fit.xtc"
echo "  ✓ rmsd_backbone.xvg, rmsd_ligand.xvg"
echo "  ✓ rmsf_calpha.xvg"
echo "  ✓ gyrate.xvg"
echo "  ✓ hbonds_count.xvg"
echo "  ✓ FINAL_RESULTS_MMPBSA.dat"
echo "  ✓ contact_frequencies.csv"
echo "  ✓ medoid_frame.pdb"
echo "  ✓ plip_output/"
echo "  ✓ summary.txt (read this for headline numbers)"
echo "============================================================"