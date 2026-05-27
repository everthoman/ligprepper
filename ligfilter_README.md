# ligfilter.py — Ligand Structural Filter

Fast multi-stage ligand filter for SMILES and SDF files.
Removes unwanted compounds based on structural alerts (PAINS, REOS), physicochemical property ranges, drug-likeness rules, and statistical outlier detection.

---

## Dependencies

| Package | Role | Required |
|---|---|---|
| RDKit | molecule I/O, SMARTS matching, descriptor calculation | yes |
| rich-argparse | coloured help | optional |
| argcomplete | tab completion | optional |

```bash
conda install -c conda-forge rdkit
pip install rich-argparse argcomplete
```

---

## Quick start

```bash
# Property filter only
ligfilter.py -i compounds.smi --ro5 --mw 150:500

# PAINS + REOS + Lipinski
ligfilter.py -i compounds.smi --pains --reos --ro5

# Outlier removal in MW/LogP/RB space (default box-plot fences, k=1.5)
ligfilter.py -i compounds.smi --outlier

# Full kitchen-sink
ligfilter.py -i compounds.smi --pains --reos --ro5 --outlier --mw 150:500 --logp -2:5 --tpsa :140
```

---

## Input formats

| Format | Detection | Notes |
|---|---|---|
| `.smi` / `.csv` / `.tsv` | extension | comma- or tab-separated; `#` lines skipped |
| `.sdf` | extension | V2000 / V3000 |

**SMILES column is auto-detected** by scanning the first data rows and picking the column where the most entries parse as valid SMILES. Override with `--smiles-col`.

---

## Preprocessing (all on by default)

Applied in order before any filter runs.

| Step | Flag to disable |
|---|---|
| Salt stripping — keep the largest fragment by heavy-atom count | `--no-strip` |
| Charge neutralization — RDKit Uncharger (quaternary N left intact) | `--no-neutralize` |
| Deduplication — canonical SMILES; keep lexicographically smallest ID | `--no-unique` |

---

## Filters

Filters are applied in the order listed below.  The first failure short-circuits — a rejected molecule is not tested against later filters.

### Structural alerts

| Flag | Description |
|---|---|
| `--pains` | Pan-Assay INterference compoundS (SMARTS from `PAINS.txt`) |
| `--pains-file FILE` | Custom PAINS file — format: `name<TAB>SMARTS` |
| `--reos` | Rapid Elimination Of Swill functional group filter (`REOS.txt`) |
| `--reos-file FILE` | Custom REOS file — format: `SMARTS<TAB>max_count<TAB>description` |
| `--custom` | User-defined SMARTS rules from `custom_filters.txt` (same format as REOS) |
| `--custom-file FILE` | Custom filter file (implies `--custom`) |

### Property ranges

All accept `MIN:MAX`, `MIN:`, `:MAX`, or an exact value.

| Flag | Property |
|---|---|
| `--mw RANGE` | Molecular weight (average) |
| `--logp RANGE` | Wildman-Crippen logP |
| `--hba RANGE` | H-bond acceptors (N + O) |
| `--hbd RANGE` | H-bond donors (NH + OH) |
| `--rb RANGE` | Rotatable bonds |
| `--tpsa RANGE` | Topological polar surface area (Å²) |
| `--qed RANGE` | Quantitative Estimate of Drug-likeness (0–1) |
| `--chiral RANGE` | Stereocentre count |
| `--ha RANGE` | Heavy-atom count |

### Drug-likeness rules

| Flag | Rule |
|---|---|
| `--ro5` | Lipinski Rule of Five — one violation allowed |
| `--ro5-strict` | Lipinski — zero violations allowed |
| `--ro3` | Astex Rule of Three (fragment screening): MW≤300, logP≤3, HBD≤3, HBA≤3, RotBonds≤3 |

### Outlier detection

Tukey IQR fences computed from the **preprocessed molecule set** — bounds reflect the actual input distribution.

| Flag | Default | Description |
|---|---|---|
| `--outlier` | off | Enable outlier removal |
| `--outlier-iqr K` | `1.5` | Fence multiplier: Q1 − k·IQR to Q3 + k·IQR.  Use `3.0` for extreme outliers only |
| `--outlier-props PROPS` | `mw,logp,rb` | Comma-separated properties to check.  Valid: `mw logp hba hbd rb tpsa qed chiral ha` |

---

## Output

| File | Content |
|---|---|
| `<input>_filtered.<ext>` | Passing molecules (same format as input, or use `-o`) |
| `<output>.log` | Mirror of stdout including rejection reasons and property statistics |

SMILES output uses canonical isomeric SMILES.  SDF output auto-generates 2D coordinates.  Extra input columns (scores, IDs) are preserved in the output.

---

## Summary block

At the end of each run the script prints:

- Molecules read / salts stripped / neutralized / duplicates removed
- Passed / rejected counts
- Rejection breakdown by reason (sorted by frequency)
- Property statistics (mean, median, std, min, max) for passing and rejected sets separately

---

## Performance

50K compounds on a 12-core workstation (indicative):

| Filters active | Time |
|---|---|
| Property ranges only | ~15–30 s |
| + PAINS (480 patterns) | ~1–3 min |
| + REOS | ~2–4 min |

Preprocessing (salt strip, neutralize, dedup) is single-threaded.  Filtering is parallelised across all available CPUs (`--jobs N` to override).

---

## Filter file formats

**PAINS** (`PAINS.txt`):
```
name<TAB>SMARTS
azo_A	[N;R0][N;R0]
...
```

**REOS / custom** (`REOS.txt`, `custom_filters.txt`):
```
SMARTS<TAB>max_count<TAB>description
[#6;H0;D3]	0	"quaternary carbon"
...
```
`max_count = 0` — reject if any match found.  
`max_count = N` — reject if more than N matches found.
