# 04 — Molecular Dynamics Simulation

**Scripts:** `01_setup.sh` → `02_equilibrate.sh` → `03_production.sh` → `04_analyze.sh`  
**Platform:** Param Smriti HPC (NABI Mohali) — GPU partition (Tesla V100)  
**Conda environment:** `md`  
**Scheduler:** SLURM  
**Working directory:** `/scratch/amans/mds/runs/<sysname>/`

---

## What this module does

Full protein-ligand MD simulation pipeline validating top PMF leads from
the ML pipeline against key MPMA cascade targets. Each system runs through
setup, equilibration, 100 ns production, and complete trajectory analysis
in four sequential steps.

---

## Pipeline overview

| Script | Phase | Compute | Time per system |
|---|---|---|---|
| `01_setup.sh` | Ligand parameterization, solvation, ionization | CPU | 5-15 min |
| `02_equilibrate.sh` | EM + NVT + NPT equilibration | GPU | 5-10 min |
| `03_production.sh` | 100 ns production MD | GPU | 4-7 hours |
| `04_analyze.sh` | RMSD, RMSF, Rg, H-bonds, MM-GBSA, contacts, PLIP | CPU | 30-60 min |

Total per system: approximately 5-9 hours wallclock on V100.

---

## Systems simulated

| PDB file | System | Role |
|---|---|---|
| `bcl2_p12.pdb` | BCL2 + PMF-12 | Top lead |
| `bcl2pc.pdb` | BCL2 + Venetoclax | Positive control |
| `hspp24.pdb` | HSP90 + PMF-24 | High-efficiency reserve |
| `hsppc.pdb` | HSP90 + Geldanamycin | Positive control |
| `sirt3p9.pdb` | SIRT3 + PMF-09 | Activate target |
| `sirt3pc.pdb` | SIRT3 + positive control | Positive control |

Input PDBs located in `/scratch/amans/mds/inputs/`

---

## Force field stack

| Component | Choice | Reasoning |
|---|---|---|
| Protein | AMBER99SB-ILDN | Standard for GROMACS protein-ligand MD; published precedent in CBM |
| Ligand | GAFF2 | Most validated general force field for drug-like molecules |
| Charges | AM1-BCC (antechamber) | Standard semi-empirical charge model |
| Water | TIP3P | Compatible with AMBER force fields |
| Box | Rhombic dodecahedron, 1.2 nm buffer | Minimizes water count (Lemkul 2024 tutorial) |
| Salt | 0.15 M NaCl, neutralized | Physiological conditions |
| Thermostat | V-rescale (Bussi) | Modern GROMACS recommendation |
| Barostat NPT | C-rescale | Avoids large equilibration oscillations |
| Barostat production | Parrinello-Rahman | Current best practice for production |
| Electrostatics | PME, rcoulomb = 1.0 nm | Standard |
| Constraints | LINCS h-bonds | Enables 2 fs timestep |
| Production length | 100 ns | Standard for binding stability claims in CBM / Eur J Med Chem |

---

## Requirements

### Software (all in conda environment `md`)

| Tool | Version | Purpose |
|---|---|---|
| GROMACS | 2024.5 CUDA build | MD engine |
| AmberTools (antechamber, parmchk2, tleap) | — | Ligand parameterization |
| acpype | — | AMBER to GROMACS format conversion |
| OpenBabel | — | PDB to MOL2 conversion, H placement |
| gmx_MMPBSA | 1.6.4 | Binding free energy |
| MDAnalysis | 2.10+ | Trajectory analysis |
| PLIP | 3.0+ | Protein-ligand interaction profiling |
| Python | 3.11+ | — |

### Conda environment setup (one-time)

```bash
mamba create -n md python=3.11
CONDA_OVERRIDE_CUDA="12.0" mamba install -c conda-forge \
    "gromacs=2024.5=nompi_cuda_h5cb645a_0" -y
mamba install -c conda-forge ambertools openmpi mpi4py gcc gxx -y
pip install Cython>=3.0.0
pip install gmx_MMPBSA==1.6.4 --no-build-isolation
mamba install -c conda-forge MDAnalysis plip -y
```

### Critical - verify GROMACS is CUDA-capable before running

```bash
gmx --version | grep "GPU support"
# Must show: GPU support: CUDA   (NOT OpenCL)
```

---

## Usage

### Standard workflow (one system)

```bash
# Step 1 - setup (login node or CPU job)
./01_setup.sh /path/to/merged.pdb LIG_RESNAME sysname

# Step 2 - equilibration (GPU node)
sbatch --account=YOUR_ACCT --partition=gpu --gres=gpu:1 \
       --time=00:30:00 --cpus-per-task=4 \
       --output=runs/<sysname>/equil_%j.out \
       --wrap="./02_equilibrate.sh <sysname>"

# Step 3 - production (GPU node, 4-7 hours)
sbatch --account=YOUR_ACCT --partition=gpu --gres=gpu:1 \
       --time=08:00:00 --cpus-per-task=4 \
       --output=runs/<sysname>/prod_%j.out \
       --wrap="./03_production.sh <sysname>"

# Step 4 - analysis (CPU, 30-60 min)
sbatch --account=YOUR_ACCT --partition=cpu \
       --time=01:30:00 --cpus-per-task=4 \
       --output=runs/<sysname>/analyze_%j.out \
       --wrap="./04_analyze.sh <sysname>"
```

Note for Param Smriti users: do not include `--mem=` flag. PARAMSMRITI has
`DefMemPerNode=UNLIMITED` and rejects explicit memory requests.

---

## Script-by-script reference

### 01_setup.sh - System preparation

```bash
./01_setup.sh <merged.pdb> <ligand_resname> <sysname>
# Example:
./01_setup.sh ~/inputs/bcl2_p12.pdb LIG bcl2_p12
```

The merged PDB must contain protein as `ATOM` records and ligand as `HETATM`
records with a consistent residue name. The script normalizes the input
resname to `LIG` internally at the antechamber step - any 3-letter input
resname is accepted.

Stages: split protein and ligand, strip buffers, parameterize ligand
(OpenBabel to antechamber to parmchk2 to tleap to acpype), build protein
topology with pdb2gmx, assemble complex, define box, solvate, add ions,
build index file, verify with dry-run grompp.

**Auto-stripped silently with notice:**

| Resname | Description |
|---|---|
| HOH | Crystallographic waters |
| PEG, EDO | Polyethylene glycol / ethylene glycol |
| GOL | Glycerol |
| MES | MES buffer |
| DMS | DMSO |

**Auto-detected and HALTED - manual intervention required:**

| Case | Residues | Required action |
|---|---|---|
| Metalloproteins | Zn, Mg, Fe, Mn, Ni, Cu, Ca | Decide bonded vs non-bonded; consider CYM for Zn-coordinating Cys |
| Cofactors | NAD, ATP, ADP, FAD, FMN, GTP, GDP, COA, SAM, HEM | Separate parameterization via antechamber + acpype |

Multi-chain proteins: detected, processing continues with a warning. May
fail at pdb2gmx - if so, run pdb2gmx interactively with `-merge interactive`.

**Key outputs:** `complex_ionized.gro`, `topol.top`, `index.ndx`, `LIG.amb2gmx/`

---

### 02_equilibrate.sh - EM + NVT + NPT

```bash
./02_equilibrate.sh <sysname>
```

| Stage | Protocol | Duration | Settings |
|---|---|---|---|
| Energy minimization | Steepest descent | Until Fmax < 1000 kJ/mol/nm | max 50,000 steps |
| NVT | V-rescale thermostat | 100 ps | 310 K, position restraints on |
| NPT | C-rescale barostat | 1 ns | 310 K, 1 atm, position restraints on |

**Key outputs:** `em.gro`, `nvt.gro`, `nvt.cpt`, `npt.gro`, `npt.cpt`

---

### 03_production.sh - 100 ns production MD

```bash
./03_production.sh <sysname>
```

100 ns total (50,000,000 steps x 2 fs), Parrinello-Rahman barostat, no
position restraints, trajectory every 100 ps = 1000 frames, GPU-accelerated
(`-ntmpi 1 -nb gpu -pme gpu`).

**Key outputs:** `production.xtc`, `production.gro`, `production.cpt`

---

### 04_analyze.sh - Full trajectory analysis

```bash
./04_analyze.sh <sysname>
```

**Stage 1 - Trajectory preprocessing (3-step trjconv)**
nojump on full system, center Protein in box, fit on Backbone (rot+trans).
Output: `production_fit.xtc`

**Stage 2 - RMSD / RMSF / Rg / H-bonds (MDAnalysis)**

MDAnalysis used over native GROMACS tools to avoid PBC wrap artifacts in
`gmx rms/rmsf` and failure to auto-detect donors/acceptors for non-standard
ligands in `gmx hbond`.

H-bond criteria: donor-acceptor distance <= 3.5 A, D-H-A angle >= 120 degrees.

**Stage 3 - MM-GBSA (last 50 ns, frames 500-1000)**
GB model igb=8 (neck GB, most accurate for protein-ligand), 0.15 M salt,
310 K. Output: `FINAL_RESULTS_MMPBSA.dat`

**Stage 4 - Per-residue contact frequencies + medoid frame**
Binary residue-level counting at 4.0 A cutoff. Avoids atom-level
multi-counting artifact. Reports residues in contact > 10% of simulation
time. Medoid frame = minimum sum-RMSD across all frames, used for PLIP.

**Stage 5 - PLIP interaction diagram**
Runs on medoid frame PDB, outputs `.txt` and `.xml` reports to `plip_output/`.

**Stage 6 - Summary report**
`summary.txt` - read this first for any system.

---

## Output structure per system

```
runs/<sysname>/
├── complex_ionized.gro       # solvated + ionized system
├── topol.top                 # combined topology
├── index.ndx                 # named index groups
├── LIG.amb2gmx/              # ligand topology (GROMACS format)
├── em.gro                    # energy-minimized
├── nvt.gro / nvt.cpt         # NVT equilibrated
├── npt.gro / npt.cpt         # NPT equilibrated
├── production.xtc            # 100 ns raw trajectory
├── production_fit.xtc        # PBC-corrected, centered, backbone-fitted
├── rmsd_backbone.xvg         # backbone RMSD (nm vs ns)
├── rmsd_ligand.xvg           # ligand RMSD (nm vs ns)
├── rmsf_calpha.xvg           # per-residue Ca RMSF (nm)
├── gyrate.xvg                # radius of gyration (nm vs ns)
├── hbonds_count.xvg          # protein-ligand H-bond count per frame
├── FINAL_RESULTS_MMPBSA.dat  # MM-GBSA dG_bind decomposition
├── contact_frequencies.csv   # per-residue contact percentage
├── medoid_frame.pdb          # representative structure for PLIP
├── plip_output/              # PLIP interaction reports (txt + xml)
└── summary.txt               # headline numbers - read this first
```

---

## Common errors and fixes

**gmx shows OpenCL instead of CUDA**
```bash
mamba uninstall gromacs -y
CONDA_OVERRIDE_CUDA="12.0" mamba install -c conda-forge \
    "gromacs=2024.5=nompi_cuda_h5cb645a_0" -y
```

**pdb2gmx fails due to leftover HETATM records**
```bash
grep -v "^HETATM.* WEIRD " input.pdb > cleaned.pdb
```

**trjconv: inconsistent shifts over periodic boundaries**
Known PBC artifact. Do not run `gmx rms` directly on raw trajectory.
`04_analyze.sh` uses MDAnalysis to avoid this.

**gmx hbond: Selection LIG has no donors/acceptors**
GROMACS does not auto-detect donors for non-standard ligands. Already
handled via MDAnalysis HydrogenBondAnalysis in `04_analyze.sh`.

**SLURM rejects --mem= flag**
Param Smriti uses `DefMemPerNode=UNLIMITED`. Omit `--mem=` entirely.

**Long script pastes corrupted over SSH**
Use `nano` and paste in chunks, or scp the file directly from your laptop.

---

## Limitations

This pipeline covers standard protein-ligand binding stability analysis
at the level expected for natural-product MD papers in CBM / J Med Chem /
Eur J Med Chem. It does not support FEP, enhanced sampling, multi-replica
statistics, membrane proteins, covalent inhibitors, or protein-protein
interaction studies.

---

## Citations

- **GROMACS:** Abraham et al. (2015) SoftwareX 1:19-25
- **AMBER99SB-ILDN:** Lindorff-Larsen et al. (2010) Proteins 78:1950-1958
- **GAFF2 / antechamber:** Wang et al. (2004) J Comput Chem 25:1157-1174
- **AM1-BCC charges:** Jakalian et al. (2002) J Comput Chem 23:1623-1641
- **AmberTools:** Case et al. (2024) Amber 2024, UCSF
- **acpype:** Sousa da Silva and Vranken (2012) BMC Res Notes 5:367
- **MDAnalysis:** Michaud-Agrawal et al. (2011) J Comput Chem 32:2319-2327
- **gmx_MMPBSA:** Valdes-Tresanco et al. (2021) J Chem Theory Comput 17:6281-6291
- **PLIP:** Schake et al. (2025) Nucleic Acids Research, gkaf361
- **TIP3P:** Jorgensen et al. (1983) J Chem Phys 79:926-935

---

## Contact

Pipeline developed by Aman Sharma (NABI Mohali / Panjab University) for
research on *Kaempferia parviflora* PMFs targeting the Mitochondrial
Fragmentation-to-Metastasis (MFM) axis in colorectal cancer.

Issues or suggestions: open a GitHub issue on this repository.
