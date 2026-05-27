# ligprep.py — Ligand Preparation

Open-source equivalent of Schrödinger LigPrep/Epik for structure-based drug design workflows.
Uses ChemAxon cxcalc for combined tautomer + protomer enumeration at target pH and CDPKit
CONFORGE for high-quality 3D conformer generation, with optional xTB refinement.

CONFORGE is ~50–100× faster than ETKDGv3 and produces correct ring geometries (chair, not boat/twist-boat) without a post-hoc filter.

---

## Dependencies

| Package | Role | Required |
|---|---|---|
| RDKit | molecule I/O, stereo enumeration, MMFF94s energy evaluation | yes |
| CDPKit | CONFORGE 3D conformer generation | yes |
| cxcalc (ChemAxon JChem Suite) | tautomer + protomer distribution at target pH | yes |
| xtb (v6.x) | GFN2 / GFN-FF geometry refinement | optional (`--xtb`) |
| tqdm | progress bars | optional |
| rich-argparse | coloured help | optional |
| argcomplete | tab completion | optional |

Install in a conda environment:

```bash
conda install -c conda-forge rdkit cdpkit
pip install tqdm rich-argparse argcomplete
# cxcalc: install ChemAxon JChem Suite, ensure cxcalc is on PATH
# xtb: conda install -c conda-forge xtb
```

---

## Pipeline

```
Input SMILES / SDF
  → salt strip (largest fragment, RDKit)
  → cxcalc tautomers -D true -n true -H <pH>
      combined tautomer × protomer distribution at target pH
      states below --min-population filtered out
  → stereo enumeration (RDKit EnumerateStereoisomers, unspecified centers only)
  → CONFORGE 3D conformer generation (CDPKit)
      internal: H addition + MMFF94 minimization + energy window filter + RMSD deduplication
      knowledge-based torsions from CSD/PDB → correct chair geometry by construction
  → dominant/states: lowest-energy conformer per state, optional xTB GFN2 refinement
  → ensemble: optional xTB GFN-FF → RMSD pruning → output
  → output SDF
```

### Why CONFORGE instead of ETKDGv3?

CONFORGE samples torsion angles from a library of experimental CSD/PDB geometries rather than
distance-geometry embedding. This gives correct ring conformations (chair, not boat/twist-boat)
by construction, without a separate filter step. Benchmark: mean RMSD to crystal structures
0.63 Å at 50 conformers, vs 0.76 Å for ETKDGv3 (Seidel et al. 2023), and ~50–100× faster.

---

## Output modes (`--mode`)

| Mode | Tautomers | Conformers | Use case |
|---|---|---|---|
| `dominant` (default) | top state only | 1 lowest-energy conformer per stereo isomer | standard docking |
| `states` | all states ≥ `--min-population` | 1 lowest-energy conformer per state | multi-state docking |
| `ensemble` | all states ≥ `--min-population` | RMSD-pruned conformer set per state | ensemble / pharmacophore |

A **state penalty** is written in all modes:

```
state_penalty = −log10(tautomer_population × stereo_weight)
```

In dominant mode `state_penalty = 0` (single state, population = 1.0). In states/ensemble mode
this penalty can be added to a raw docking score to correct for state probability, analogous to
Glide GlideScore + Epik state penalty.

---

## Usage

```
ligprep.py [-h] -i FILE [-o FILE] [--log FILE]
            [--pH FLOAT] [--min-population FLOAT] [--max-tautomers INT]
            [--mode {dominant,states,ensemble}]
            [--num-confs INT] [--conforge-ewin FLOAT] [--rmsd-threshold FLOAT]
            [--ewin FLOAT] [--max-confs-out INT]
            [--no-stereo] [--max-stereo INT]
            [-j N]
```

### Options

#### General

| Flag | Default | Description |
|---|---|---|
| `-i / --input FILE` | required | Input file: SMILES (`.smi`, `.csv`, `.txt`) or SDF (`.sdf`, `.sd`) |
| `-o / --output FILE` | `<input_stem>_prep.sdf` | Output SDF file |
| `--log FILE` | `<output_stem>.log` | Log file (stdout is tee'd here) |
| `--pH FLOAT` | `7.4` | Target pH for tautomer/protomer distribution |
| `--mode MODE` | `dominant` | Output mode: `dominant`, `states`, or `ensemble` |
| `--no-stereo` | off | Skip stereo enumeration; keep input stereo as-is |
| `--max-stereo INT` | `32` | Max stereoisomers per molecule for unspecified centers |
| `-j / --jobs N` | `cpu_count − 1` | Parallel worker processes for CONFORGE |

SMILES files may use whitespace, comma, or tab as delimiter. Second column is the molecule name.

#### Tautomers / protomers

| Flag | Default | Description |
|---|---|---|
| `--min-population FLOAT` | `0.10` | Minimum fractional population to keep a state (states/ensemble modes) |
| `--max-tautomers INT` | `200` | Maximum tautomers per molecule passed to cxcalc |

#### CONFORGE conformer generation

| Flag | Default | Description |
|---|---|---|
| `--num-confs INT` | 50 / 100 | Max conformers CONFORGE may output per isomer (50 for dominant/states, 100 for ensemble) |
| `--conforge-ewin FLOAT` | `5.0` | CONFORGE energy window kcal/mol for dominant/states mode |
| `--rmsd-threshold FLOAT` | `0.5` | RMSD threshold Å for CONFORGE deduplication |

#### Ensemble mode

| Flag | Default | Description |
|---|---|---|
| `--ewin FLOAT` | `5.0` | CONFORGE energy window kcal/mol for ensemble mode |
| `--max-confs-out INT` | `50` | Max conformers to write per tautomer/stereo state; lowest-energy kept first |

---

## SDF output properties

| Property | Modes | Description |
|---|---|---|
| `_Name` | all | Original molecule name (same across all states of the same molecule) |
| `tautomer_id` | all | 1-based index of tautomer/protomer state |
| `tautomer_population` | all | Fractional population at target pH (0–1) |
| `stereo_id` | all | 1-based stereo isomer index |
| `state_population` | all | `tautomer_population × stereo_weight` |
| `state_penalty` | all | `−log10(state_population)`; add to docking score to correct for state probability |
| `mmff_energy` | all | MMFF94s single-point energy (kcal/mol); arbitrary zero reference, use for relative ranking only |
| `conf_id` | ensemble | 1-based conformer index within the tautomer/stereo state |

---

## Examples

### Standard docking prep (dominant state, best conformer)
```bash
ligprep.py -i compounds.smi -o compounds_prep.sdf
```

### Multi-state docking (all protomers/tautomers ≥ 10%)
```bash
ligprep.py -i compounds.smi -o compounds_prep.sdf --mode states
```

### Ensemble docking (up to 10 conformers per state, 5 kcal/mol energy window)
```bash
ligprep.py -i compounds.smi -o compounds_prep.sdf --mode ensemble --max-confs-out 10
```

### Non-physiological pH (e.g. lysosome)
```bash
ligprep.py -i compounds.smi -o compounds_prep.sdf --pH 5.0
```

### Large library with parallel CONFORGE
```bash
ligprep.py -i library.smi -o library_prep.sdf -j 32
```

---

## Notes

- **cxcalc tautomers with `-D true -H pH`** gives a combined tautomer × protomer Boltzmann
  distribution at target pH in a single call (not sequential tautomer then protonation).
  This is the key advantage over OpenBabel/Dimorphite-DL protonation tools.

- **CONFORGE energy window**: `--conforge-ewin` (dominant/states) and `--ewin` (ensemble) both
  default to 5.0 kcal/mol. For dominant/states only the lowest-energy conformer is kept, so a
  wider window can be used without penalty if challenging geometries are missed; for ensemble the
  window directly controls how many conformers survive.

- **MMFF energies are large positive numbers** by design: MMFF has an arbitrary zero reference
  and is not a thermodynamic energy. Only relative values (within the same molecule) are meaningful.

- **`--max-confs-out` is per tautomer/stereo state**, not per molecule. A molecule with
  two protomeric states will produce up to `2 × max-confs-out` records in the output SDF.

---

## Reference

Seidel T, Baber J, Schindler J, Langer T.  
*CONFORGE: A Highly Versatile and Performant Conformer Generator.*  
J. Chem. Inf. Model. 2023, **63**, 5549–5570.  
DOI: [10.1021/acs.jcim.3c00563](https://doi.org/10.1021/acs.jcim.3c00563)
