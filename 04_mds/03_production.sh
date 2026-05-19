#!/bin/bash
# ============================================================
# 03_production.sh - 100 ns production MD
# ============================================================
# Takes an NPT-equilibrated system and runs 100 ns production.
# Parrinello-Rahman barostat, no position restraints.
#
# USAGE: ./03_production.sh <sysname>
# EXAMPLE: ./03_production.sh test_bcl2pc
#
# REQUIREMENTS:
#   - 02_equilibrate.sh must have completed for this sysname
#   - GPU partition with CUDA-capable GROMACS
#   - Walltime: ~4-7 hours per system on V100
#
# OUTPUT:
#   production.xtc - 100 ns trajectory (1000 frames at 100 ps)
#   production.gro - Final structure
#   production.cpt - Checkpoint for restart if needed
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

for f in npt.gro npt.cpt topol.top index.ndx; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing $f - run 02_equilibrate.sh first"
        exit 1
    fi
done

source /home/amans/mambaforge/etc/profile.d/conda.sh
conda activate md

echo "============================================================"
echo "03_production.sh - $SYSNAME"
echo "Started: $(date)"
echo "============================================================"

cat > production.mdp << 'EOF'
title                   = Production MD 100 ns at 310 K, 1 atm
integrator              = md
nsteps                  = 50000000
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 50000
nstenergy               = 5000
nstlog                  = 50000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005
coulombtype             = PME
rcoulomb                = 1.0
fourierspacing          = 0.16
pme_order               = 4
vdwtype                 = cutoff
vdw-modifier            = potential-shift
rvdw                    = 1.0
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = Protein_LIG  Water_NA_CL
tau_t                   = 0.1     0.1
ref_t                   = 310     310
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
gen_vel                 = no
EOF

gmx grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx \
           -o production.tpr -maxwarn 2 > grompp_prod.log 2>&1
[ ! -f production.tpr ] && { echo "ERROR: grompp production failed"; tail -20 grompp_prod.log; exit 1; }

echo "[STAGE 1] Production MD starting at $(date)"
gmx mdrun -deffnm production -ntmpi 1 -nb gpu -pme gpu > production_run.log 2>&1
[ ! -f production.gro ] && { echo "ERROR: Production mdrun failed"; tail -20 production_run.log; exit 1; }

PROD_PERF=$(grep "Performance:" production.log | tail -1 | awk '{print $2}')
PROD_SIZE=$(stat -c%s production.xtc)

echo ""
echo "============================================================"
echo "03_production.sh COMPLETE: $SYSNAME"
echo "Finished: $(date)"
echo ""
echo "Output: $WORKDIR"
echo "  ✓ production.xtc ($(($PROD_SIZE / 1024 / 1024)) MB)"
echo "  ✓ production.gro"
echo "  ✓ production.cpt"
echo ""
echo "Performance: $PROD_PERF ns/day"
echo ""
echo "Next: ./04_analyze.sh $SYSNAME"
echo "============================================================"