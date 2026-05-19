#!/bin/bash
# ============================================================
# 02_equilibrate.sh - EM + NVT + NPT equilibration
# ============================================================
# Takes a system prepared by 01_setup.sh and runs:
#   - Energy minimization (steepest descent, Fmax<1000)
#   - NVT equilibration (100 ps, 310K, V-rescale, position restraints)
#   - NPT equilibration (1 ns, 1 atm, C-rescale, position restraints)
#
# USAGE: ./02_equilibrate.sh <sysname>
# EXAMPLE: ./02_equilibrate.sh test_bcl2pc
#
# REQUIREMENTS:
#   - 01_setup.sh must have completed for this sysname
#   - GPU partition with CUDA-capable GROMACS
#
# OUTPUT:
#   em.gro          - Energy-minimized structure
#   nvt.gro, nvt.cpt - NVT-equilibrated (with restart info)
#   npt.gro, npt.cpt - NPT-equilibrated (ready for production)
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
    echo "Run 01_setup.sh first"
    exit 1
fi

cd "$WORKDIR"

# Verify 01_setup.sh outputs exist
for f in complex_ionized.gro topol.top index.ndx; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing $f - run 01_setup.sh first"
        exit 1
    fi
done

source /home/amans/mambaforge/etc/profile.d/conda.sh
conda activate md

echo "============================================================"
echo "02_equilibrate.sh - $SYSNAME"
echo "Started: $(date)"
echo "============================================================"

# ============================================================
# STAGE 1: ENERGY MINIMIZATION
# ============================================================
echo ""
echo "[STAGE 1] Energy minimization (steepest descent)"

cat > em.mdp << 'EOF'
integrator              = steep
emtol                   = 1000.0
emstep                  = 0.01
nsteps                  = 50000
nstlist                 = 10
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
fourierspacing          = 0.16
pme_order               = 4
vdwtype                 = cutoff
vdw-modifier            = potential-shift
rvdw                    = 1.0
DispCorr                = EnerPres
pbc                     = xyz
EOF

gmx grompp -f em.mdp -c complex_ionized.gro -p topol.top -o em.tpr -maxwarn 2 \
    > grompp_em.log 2>&1
[ ! -f em.tpr ] && { echo "ERROR: grompp EM failed"; tail -20 grompp_em.log; exit 1; }

gmx mdrun -deffnm em -nt 4 > em_run.log 2>&1
[ ! -f em.gro ] && { echo "ERROR: EM mdrun failed"; tail -20 em_run.log; exit 1; }

EM_FMAX=$(grep "Maximum force" em.log | tail -1 | awk '{print $4}')
EM_EPOT=$(grep "Potential Energy" em.log | tail -1 | awk '{print $4}')
echo "  EM converged. Fmax=$EM_FMAX, Epot=$EM_EPOT"

# ============================================================
# STAGE 2: NVT EQUILIBRATION (100 ps, 310 K)
# ============================================================
echo ""
echo "[STAGE 2] NVT equilibration (100 ps, 310 K)"

cat > nvt.mdp << 'EOF'
title                   = NVT equilibration 100 ps at 310 K
define                  = -DPOSRES
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 5000
nstenergy               = 500
nstlog                  = 5000
continuation            = no
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
pcoupl                  = no
gen_vel                 = yes
gen_temp                = 310
gen_seed                = -1
EOF

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx \
           -o nvt.tpr -maxwarn 2 > grompp_nvt.log 2>&1
[ ! -f nvt.tpr ] && { echo "ERROR: grompp NVT failed"; tail -20 grompp_nvt.log; exit 1; }

gmx mdrun -deffnm nvt -ntmpi 1 -nb gpu -pme gpu > nvt_run.log 2>&1
[ ! -f nvt.gro ] && { echo "ERROR: NVT mdrun failed"; tail -20 nvt_run.log; exit 1; }

NVT_PERF=$(grep "Performance:" nvt.log | tail -1 | awk '{print $2}')
echo "  NVT done. Performance: $NVT_PERF ns/day"
# ============================================================
# STAGE 3: NPT EQUILIBRATION (1 ns, 1 atm, 310 K)
# ============================================================
echo ""
echo "[STAGE 3] NPT equilibration (1 ns, 1 atm)"

cat > npt.mdp << 'EOF'
title                   = NPT equilibration 1 ns at 310 K, 1 atm
define                  = -DPOSRES
integrator              = md
nsteps                  = 500000
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 5000
nstenergy               = 500
nstlog                  = 5000
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
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
gen_vel                 = no
EOF

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx \
           -o npt.tpr -maxwarn 2 > grompp_npt.log 2>&1
[ ! -f npt.tpr ] && { echo "ERROR: grompp NPT failed"; tail -20 grompp_npt.log; exit 1; }

gmx mdrun -deffnm npt -ntmpi 1 -nb gpu -pme gpu > npt_run.log 2>&1
[ ! -f npt.gro ] && { echo "ERROR: NPT mdrun failed"; tail -20 npt_run.log; exit 1; }

NPT_PERF=$(grep "Performance:" npt.log | tail -1 | awk '{print $2}')
echo "  NPT done. Performance: $NPT_PERF ns/day"

# ============================================================
# FINAL SUMMARY
# ============================================================
echo ""
echo "============================================================"
echo "02_equilibrate.sh COMPLETE: $SYSNAME"
echo "Finished: $(date)"
echo ""
echo "Output files in: $WORKDIR"
echo "  ✓ em.gro (Fmax=$EM_FMAX, Epot=$EM_EPOT)"
echo "  ✓ nvt.gro, nvt.cpt (NVT @ $NVT_PERF ns/day)"
echo "  ✓ npt.gro, npt.cpt (NPT @ $NPT_PERF ns/day)"
echo ""
echo "Next: ./03_production.sh $SYSNAME"
echo "============================================================"

