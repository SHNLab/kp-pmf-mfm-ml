#!/bin/bash
# ============================================================
# 01_setup.sh - System setup for protein-ligand MD
# ============================================================
# Takes a merged protein+ligand PDB, produces a fully prepared
# GROMACS system ready for equilibration (02_equilibrate.sh).
#
# USAGE: ./01_setup.sh <merged.pdb> <ligand_resname> <sysname>
# EXAMPLE: ./01_setup.sh ~/inputs/bcl2_p12.pdb LIG bcl2_p12
#
# Force field stack: AMBER99SB-ILDN (protein) + GAFF2/AM1-BCC (ligand) + TIP3P
# Box: rhombic dodecahedron, 1.2 nm buffer
# Ions: 0.15 M NaCl, neutralized
#
# AUTO-HANDLED: Standard protein-ligand systems
# AUTO-DETECTED + FLAGGED FOR MANUAL HANDLING:
#   - Metalloproteins (Zn, Mg, Fe, Mn, Ni, Cu, Ca)
#   - Cofactors (NAD, ATP, ADP, FAD, GTP, GDP)
#   - Multi-chain proteins
# AUTO-STRIPPED: Crystallization buffers (PEG, EDO, GOL, MES, DMS)
# ============================================================

set -e

# --- Parse arguments ---
if [ $# -ne 3 ]; then
    echo "Usage: $0 <merged.pdb> <ligand_resname> <sysname>"
    echo "Example: $0 ~/inputs/bcl2_p12.pdb LIG bcl2_p12"
    exit 1
fi

INPUT_PDB="$1"
LIG_RESNAME="$2"
SYSNAME="$3"

# --- Validate input ---
if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

# --- Set up workspace ---
WORKDIR="/scratch/amans/mds/runs/$SYSNAME"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Copy input PDB
cp "$INPUT_PDB" merged_input.pdb

# --- Activate conda env ---
source /home/amans/mambaforge/etc/profile.d/conda.sh
conda activate md

echo "============================================================"
echo "01_setup.sh - $SYSNAME"
echo "Input: $INPUT_PDB"
echo "Ligand resname: $LIG_RESNAME"
echo "Started: $(date)"
echo "============================================================"

# --- Validate ligand presence ---
LIG_LINES=$(grep -c "^HETATM.* $LIG_RESNAME " merged_input.pdb || echo 0)
if [ "$LIG_LINES" -eq 0 ]; then
    echo "ERROR: No HETATM lines found with resname '$LIG_RESNAME'"
    echo "Available HETATM resnames in input:"
    grep "^HETATM" merged_input.pdb | awk '{print $4}' | sort -u
    exit 1
fi
echo "[VALIDATE] Found $LIG_LINES ligand atom records (resname $LIG_RESNAME)"

# --- Detect special cases ---
echo ""
echo "[DETECT] Scanning for special cases..."

# Metals
METAL_FOUND=""
for METAL in ZN MG FE MN NI CU CA; do
    if grep -q "^HETATM.* $METAL " merged_input.pdb; then
        COUNT=$(grep -c "^HETATM.* $METAL " merged_input.pdb)
        echo "  ⚠️  METAL DETECTED: $METAL ($COUNT atoms)"
        METAL_FOUND="$METAL_FOUND $METAL"
    fi
done

if [ -n "$METAL_FOUND" ]; then
    echo ""
    echo "============================================================"
    echo "MANUAL INTERVENTION REQUIRED"
    echo "============================================================"
    echo "This system contains metal ions:$METAL_FOUND"
    echo ""
    echo "Metal handling requires manual configuration:"
    echo "  1. Identify metal-coordinating residues (typically Cys, His)"
    echo "  2. Decide handling strategy:"
    echo "     - Non-bonded ion (default, simplest)"
    echo "     - MCPB.py bonded approach (gold standard)"
    echo "  3. For Zn-coordinating Cys: rename to CYM (deprotonated)"
    echo ""
    echo "This automation does not handle metals. Use manual workflow."
    echo "============================================================"
    exit 2
fi

# Cofactors
COFACTOR_FOUND=""
for COFACTOR in NAD NAP NDP NAI NAM ATP ADP FAD FMN GTP GDP COA SAM HEM HEC; do
    if grep -q "^HETATM.* $COFACTOR " merged_input.pdb; then
        COUNT=$(grep -c "^HETATM.* $COFACTOR " merged_input.pdb)
        echo "  ⚠️  COFACTOR DETECTED: $COFACTOR ($COUNT atoms)"
        COFACTOR_FOUND="$COFACTOR_FOUND $COFACTOR"
    fi
done

if [ -n "$COFACTOR_FOUND" ]; then
    echo ""
    echo "============================================================"
    echo "COFACTOR(S) DETECTED:$COFACTOR_FOUND"
    echo "Cofactors require separate parameterization."
    echo "This automation does not handle cofactors."
    echo "============================================================"
    exit 2
fi

# Multi-chain
N_CHAINS=$(awk '/^ATOM/ {print substr($0,22,1)}' merged_input.pdb | sort -u | wc -l)
if [ "$N_CHAINS" -gt 1 ]; then
    echo "  ⚠️  MULTI-CHAIN: $N_CHAINS chains detected"
    echo "     Multi-chain systems may need pdb2gmx -merge interactive setup."
    echo "     Continuing — will fail loudly if multi-chain causes issues."
fi

# Crystallization buffers (auto-strip with notice)
BUFFER_FOUND=""
for BUFFER in PEG EDO GOL MES DMS HOH; do
    if grep -q "^HETATM.* $BUFFER " merged_input.pdb; then
        COUNT=$(grep -c "^HETATM.* $BUFFER " merged_input.pdb)
        echo "  ℹ️  BUFFER: $BUFFER ($COUNT atoms) - will be stripped"
        BUFFER_FOUND="$BUFFER_FOUND $BUFFER"
    fi
done

echo "[DETECT] Done."
echo ""
# ============================================================
# STAGE 1: SPLIT PROTEIN AND LIGAND
# ============================================================
echo "[STAGE 1] Splitting protein and ligand"

# Strip ANISOU lines, alt-locs other than A, and crystallization buffers
PROCESS_PDB="merged_clean.pdb"
grep -v "^ANISOU" merged_input.pdb | \
    awk 'substr($0,17,1)==" " || substr($0,17,1)=="A" {print}' > "$PROCESS_PDB"

# Strip detected buffer molecules
for BUFFER in $BUFFER_FOUND; do
    grep -v "^HETATM.* $BUFFER " "$PROCESS_PDB" > tmp.pdb && mv tmp.pdb "$PROCESS_PDB"
done

# Extract protein (ATOM records only)
grep "^ATOM" "$PROCESS_PDB" > protein.pdb
echo "  Protein atoms: $(wc -l < protein.pdb)"

# Extract ligand (HETATM with our resname)
grep "^HETATM.* $LIG_RESNAME " "$PROCESS_PDB" > ligand.pdb
echo "  Ligand atoms: $(wc -l < ligand.pdb)"

if [ ! -s protein.pdb ] || [ ! -s ligand.pdb ]; then
    echo "ERROR: Empty protein or ligand file after split"
    exit 1
fi

# ============================================================
# STAGE 2: PARAMETERIZE LIGAND (GAFF2 + AM1-BCC)
# ============================================================
echo ""
echo "[STAGE 2] Ligand parameterization"

# Convert PDB to MOL2 with hydrogens added (openbabel handles H placement)
obabel ligand.pdb -O ligand.mol2 -h --partialcharge gasteiger \
    > obabel.log 2>&1
[ ! -f ligand.mol2 ] && { echo "ERROR: obabel failed"; tail obabel.log; exit 1; }
# Force ligand residue name to "LIG" for downstream consistency
# acpype keeps the input resname unless we override it here
sed -i "s/$LIG_RESNAME/LIG/g" ligand.mol2
echo "  Renamed ligand resname '$LIG_RESNAME' → 'LIG' for downstream consistency"

# Determine net charge from openbabel output (rounded to int)
LIG_CHARGE=$(grep -i "total charge" obabel.log 2>/dev/null | tail -1 | awk '{print $NF}')
LIG_CHARGE=${LIG_CHARGE:-0}
LIG_CHARGE_INT=$(printf "%.0f" "$LIG_CHARGE")
echo "  Ligand net charge: $LIG_CHARGE_INT"

# antechamber: assign GAFF2 atom types + AM1-BCC charges
antechamber -i ligand.mol2 -fi mol2 \
            -o ligand_gaff.mol2 -fo mol2 \
            -at gaff2 -c bcc -nc $LIG_CHARGE_INT \
            -pf y -dr no \
            > antechamber.log 2>&1
[ ! -f ligand_gaff.mol2 ] && { echo "ERROR: antechamber failed"; tail -20 antechamber.log; exit 1; }
echo "  GAFF2 atom types + AM1-BCC charges assigned"

# parmchk2: check for missing GAFF2 parameters
parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod -s gaff2 \
    > parmchk2.log 2>&1
if grep -q "ATTN" ligand.frcmod; then
    echo "  ⚠️  Missing GAFF2 parameters detected:"
    grep "ATTN" ligand.frcmod
    echo "  parmchk2 will use estimated parameters - may affect accuracy"
else
    echo "  All GAFF2 parameters found"
fi

# tleap: build AMBER prmtop/inpcrd
cat > tleap.in << TLEAP
source leaprc.gaff2
LIG = loadmol2 ligand_gaff.mol2
loadamberparams ligand.frcmod
saveoff LIG ligand.lib
saveamberparm LIG ligand.prmtop ligand.inpcrd
quit
TLEAP
tleap -f tleap.in > tleap.log 2>&1
[ ! -f ligand.prmtop ] && { echo "ERROR: tleap failed"; tail -20 tleap.log; exit 1; }

# acpype: convert AMBER → GROMACS format
acpype -p ligand.prmtop -x ligand.inpcrd -b LIG > acpype.log 2>&1
[ ! -f LIG.amb2gmx/LIG_GMX.gro ] && { echo "ERROR: acpype failed"; tail -20 acpype.log; exit 1; }

LIG_ATOMS=$(sed -n '2p' LIG.amb2gmx/LIG_GMX.gro | awk '{print $1}')
echo "  Ligand parameterization complete: $LIG_ATOMS atoms"

# ============================================================
# STAGE 3: PROTEIN TOPOLOGY (AMBER99SB-ILDN + TIP3P)
# ============================================================
echo ""
echo "[STAGE 3] Protein topology"

gmx pdb2gmx -f protein.pdb \
            -o protein_processed.gro \
            -p protein.top \
            -i posre_prot.itp \
            -ff amber99sb-ildn \
            -water tip3p \
            -ignh \
            > pdb2gmx.log 2>&1
[ ! -f protein_processed.gro ] && { echo "ERROR: pdb2gmx failed"; tail -20 pdb2gmx.log; exit 1; }

PROT_ATOMS=$(sed -n '2p' protein_processed.gro | awk '{print $1}')
echo "  Protein topology generated: $PROT_ATOMS atoms"

# ============================================================
# STAGE 4: COMBINE PROTEIN + LIGAND, BOX, SOLVATE, IONS
# ============================================================
echo ""
echo "[STAGE 4] System assembly"

# Combine protein + ligand into complex.gro
TOTAL_ATOMS=$((PROT_ATOMS + LIG_ATOMS))
{
    head -1 protein_processed.gro
    echo "$TOTAL_ATOMS"
    sed -n "3,$((PROT_ATOMS + 2))p" protein_processed.gro
    sed -n "3,$((LIG_ATOMS + 2))p" LIG.amb2gmx/LIG_GMX.gro
    tail -1 protein_processed.gro
} > complex.gro
echo "  Combined: $TOTAL_ATOMS atoms"

# Build combined topology
cp protein.top topol.top

python3 << 'PYEOF'
import re
with open('LIG.amb2gmx/LIG_GMX.top') as f:
    content = f.read()
at_match = re.search(r'\[ atomtypes \](.*?)(?=\[ \w+ \])', content, re.DOTALL)
mt_match = re.search(r'(\[ moleculetype \].*?)(?=\[ system \])', content, re.DOTALL)
open('ligand_atomtypes.itp', 'w').write(at_match.group(0) if at_match else '')
open('ligand.itp', 'w').write(mt_match.group(1) if mt_match else '')
PYEOF

python3 << 'PYEOF'
with open('topol.top') as f:
    content = f.read()
ff_pos = content.find('#include')
end = content.find('\n', ff_pos)
inc = '\n; Include ligand atomtypes\n#include "ligand_atomtypes.itp"\n; Include ligand topology\n#include "ligand.itp"\n'
content = content[:end+1] + inc + content[end+1:]
content = content.rstrip() + '\nLIG    1\n'
open('topol.top', 'w').write(content)
PYEOF

# Define box (rhombic dodecahedron, 1.2 nm buffer)
gmx editconf -f complex.gro -o complex_box.gro -c -d 1.2 -bt dodecahedron \
    > editconf.log 2>&1

# Solvate with TIP3P
gmx solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p topol.top \
    > solvate.log 2>&1
SOLV_ATOMS=$(sed -n '2p' complex_solv.gro | awk '{print $1}')
echo "  Solvated: $SOLV_ATOMS atoms"

# Add ions to neutralize + 0.15 M NaCl
cat > ions.mdp << 'IONS'
integrator              = steep
emtol                   = 1000.0
nsteps                  = 50000
nstlist                 = 10
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
pbc                     = xyz
IONS

gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 2 \
    > grompp_ions.log 2>&1
[ ! -f ions.tpr ] && { echo "ERROR: grompp ions failed"; tail -20 grompp_ions.log; exit 1; }

echo "SOL" | gmx genion -s ions.tpr -o complex_ionized.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15 \
    > genion.log 2>&1
[ ! -f complex_ionized.gro ] && { echo "ERROR: genion failed"; tail -20 genion.log; exit 1; }

FINAL_ATOMS=$(sed -n '2p' complex_ionized.gro | awk '{print $1}')
IONS_LINE=$(grep "Will try to add" genion.log)
echo "  Ionized: $FINAL_ATOMS atoms"
echo "  $IONS_LINE"

# ============================================================
# STAGE 5: BUILD INDEX FILE WITH NAMED GROUPS
# ============================================================
echo ""
echo "[STAGE 5] Building index file"

gmx make_ndx -f complex_ionized.gro -o index.ndx << 'NDX' > make_ndx.log 2>&1
"Protein" | "LIG"
"Water" | "NA" | "CL"
q
NDX

if ! grep -q "Protein_LIG" index.ndx; then
    echo "ERROR: Protein_LIG group not created"
    cat make_ndx.log
    exit 1
fi
echo "  index.ndx created with Protein_LIG and Water_NA_CL groups"

# ============================================================
# FINAL VERIFICATION
# ============================================================
echo ""
echo "[VERIFY] Test grompp dry-run"
gmx grompp -f ions.mdp -c complex_ionized.gro -p topol.top -o test.tpr \
    -n index.ndx -maxwarn 2 > grompp_test.log 2>&1
if [ ! -f test.tpr ]; then
    echo "ERROR: Verification grompp failed"
    tail -30 grompp_test.log
    exit 1
fi
rm -f test.tpr

echo ""
echo "============================================================"
echo "01_setup.sh COMPLETE: $SYSNAME"
echo "Finished: $(date)"
echo ""
echo "Output files in: $WORKDIR"
echo "  ✓ complex_ionized.gro ($FINAL_ATOMS atoms)"
echo "  ✓ topol.top"
echo "  ✓ index.ndx (Protein_LIG, Water_NA_CL)"
echo "  ✓ LIG.amb2gmx/ (ligand topology)"
echo ""
echo "Next: ./02_equilibrate.sh $SYSNAME"
echo "============================================================"