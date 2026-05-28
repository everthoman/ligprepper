#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
Ligand Filter — remove molecules with unwanted structural features

Reads molecules from a SMILES (.smi / .csv / .tsv) or SDF file and writes
only the molecules that pass all enabled filters to the output file.
The output format matches the input format unless --output is given with an
explicit extension (.smi or .sdf).

Workflow
--------
1. Read input molecules (SMILES or SDF)
2. Preprocessing (applied in order before filtering; all three are ON by default):
   --strip / --no-strip         strip salts / small fragments, keep the largest fragment
   --neutralize / --no-neutralize  neutralize formal charges (e.g. carboxylate → acid,
                                ammonium → amine) using RDKit MolStandardize Uncharger
   --unique / --no-unique       deduplicate on canonical SMILES; for duplicates keep the
                                entry with the lexicographically smallest identifier
3. Apply each enabled filter in order; the first failure short-circuits
4. Write passing molecules to output; log reason for each rejection

Inputs
------
  SMILES file : first column = SMILES, optional second column = name
                (comma- or tab-separated, lines starting with '#' skipped)
  SDF file    : standard V2000/V3000 SD file

Outputs
-------
  Filtered SMILES / SDF file (same format as input by default)
    SMILES output : canonical isomeric SMILES (tab-separated name)
    SDF output    : clean 2D coordinates generated automatically
  Summary printed to stdout

Dependencies
------------
  rdkit          pip install rdkit   (molecule I/O and SMARTS matching)

Written by Claude Sonnet 4.6, 2026-02-27
"""

import argparse
import math
import os
import statistics
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

try:
    from rich_argparse import RawDescriptionRichHelpFormatter as _HelpFmt
    _RICH = True
except ImportError:
    _HelpFmt = argparse.RawDescriptionHelpFormatter
    _RICH = False

try:
    import argcomplete
    from argcomplete.completers import FilesCompleter
    _ARGCOMPLETE = True
except ImportError:
    _ARGCOMPLETE = False

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors, Descriptors, Crippen, QED, AllChem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _RDKIT = True
except ImportError:
    _RDKIT = False

try:
    from tqdm import tqdm
    _TQDM = True
except ImportError:
    _TQDM = False

# Default filter-file locations (same directory as this script)
_SCRIPT_DIR = Path(__file__).resolve().parent
_DEFAULT_PAINS   = _SCRIPT_DIR / 'PAINS.txt'
_DEFAULT_REOS    = _SCRIPT_DIR / 'REOS.txt'
_DEFAULT_CUSTOM  = _SCRIPT_DIR / 'custom_filters.txt'


# ─── Logging helpers ──────────────────────────────────────────────────────────

class _Tee:
    """Write to both a file and the original stdout simultaneously."""
    def __init__(self, fh, original):
        self._fh = fh
        self._orig = original
    def write(self, data):
        self._orig.write(data)
        self._fh.write(data)
    def flush(self):
        self._orig.flush()
        self._fh.flush()
    def isatty(self):
        return self._orig.isatty()


def _ok(msg: str):   print(f"  ✓ {msg}")
def _info(msg: str): print(f"    {msg}")
def _warn(msg: str): print(f"  ⚠ {msg}")

def _fatal(msg: str):
    print(f"  ✗ {msg}", file=sys.stderr)
    sys.exit(1)


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s}s"
    h, m = divmod(m, 60)
    return f"{h}h {m}m {s}s"


# ─── Molecule I/O ─────────────────────────────────────────────────────────────

_HEADER_KEYWORDS = {'smiles', 'smi', 'smile', 'mol', 'structure',
                    'canonical_smiles', 'cansmi', 'isomeric_smiles'}


def _split_smiles_line(line: str) -> List[str]:
    """Split a SMILES-file line into stripped fields.

    Splits on tab when present (multi-column TSV); otherwise on commas while
    treating '|...|' CXSMILES enhanced-stereo blocks as a single token so
    embedded commas (e.g. '|&1:3,7,r|') survive into the SMILES field.
    """
    if '\t' in line:
        return [p.strip() for p in line.split('\t')]
    if ',' not in line:
        return [line.strip()]
    parts, buf, in_cx = [], [], False
    for ch in line:
        if ch == '|':
            in_cx = not in_cx
            buf.append(ch)
        elif ch == ',' and not in_cx:
            parts.append(''.join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    parts.append(''.join(buf).strip())
    return parts


def _detect_smiles_col(path: Path) -> int:
    """Scan the first data rows and return the 0-based index of the SMILES column.

    For each row that contains at least one valid SMILES, tallies which columns
    parse successfully.  The column with the most hits wins.
    """
    counts: Dict[int, int] = {}
    n_data = 0
    with path.open(encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('==>') and line.endswith('<=='):
                continue
            parts = _split_smiles_line(line)
            valid = {i: Chem.MolFromSmiles(v) is not None for i, v in enumerate(parts)}
            if not any(valid.values()):
                continue  # header or unparseable row
            n_data += 1
            for i, ok in valid.items():
                if ok:
                    counts[i] = counts.get(i, 0) + 1
            if n_data >= 10:
                break
    if not counts:
        _warn("Could not auto-detect SMILES column; assuming column 1")
        return 0
    detected = max(counts, key=counts.get)
    _info(f"SMILES column auto-detected: column {detected + 1}")
    return detected

def _parse_col_spec(s: str):
    """Parse a column specifier: positive integer string → 0-based int; else column name str."""
    try:
        n = int(s)
        if n < 1:
            _fatal(f"Column index must be ≥ 1, got {n}")
        return n - 1  # convert 1-based user input to 0-based
    except ValueError:
        return s  # column header name


def _resolve_col(spec, header: List[str], label: str) -> int:
    """Resolve a str column name to its 0-based index in *header*. Fatal if not found."""
    if spec not in header:
        _fatal(f"--{label}: column '{spec}' not found in header: {header}")
    return header.index(spec)


def _iter_smiles(path: Path, smiles_col=None, name_col=None,
                 out_header: Optional[List[str]] = None) -> Iterator[Tuple[Chem.Mol, str]]:
    """Yield (mol, name) from a SMILES file (comma- or tab-separated).

    smiles_col : 0-based int or str column name (default: 0 = first column)
    name_col   : 0-based int, str column name, or None (auto = first non-SMILES column)
    out_header : if a list is provided and a header row is detected, it will be populated
                 with column names in output order [smiles, name, *extras] before the
                 first molecule is yielded.

    Extra columns (neither SMILES nor name) are stored as mol properties and
    written back on output.  If a header row is present, column names are used
    as property keys; otherwise col_N (1-based) is used.
    """
    RDLogger.DisableLog('rdApp.*')
    if smiles_col is None:
        smiles_col = _detect_smiles_col(path)
    need_header = isinstance(smiles_col, str) or isinstance(name_col, str)
    col_names: List[str] = []        # populated when a header row is detected
    header_consumed = False
    smi_idx: int = smiles_col if isinstance(smiles_col, int) else 0  # placeholder
    nam_idx: Optional[int] = name_col if isinstance(name_col, int) else None
    auto_name: bool = (name_col is None)
    out_header_written = False

    with path.open(encoding='utf-8', errors='replace') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            # `head`/`tail`-concatenation markers — silently skip wherever they appear
            if line.startswith('==>') and line.endswith('<=='):
                _info(f"Line {lineno}: skipping marker ({line!r})")
                continue
            parts = _split_smiles_line(line)

            # ── Header handling ───────────────────────────────────────────
            if not header_consumed:
                if need_header:
                    col_names = parts
                    smi_idx = _resolve_col(smiles_col, col_names, 'smiles-col')
                    if name_col is not None:
                        nam_idx = _resolve_col(name_col, col_names, 'name-col')
                    header_consumed = True
                    continue
                else:
                    test = parts[smi_idx] if smi_idx < len(parts) else ''
                    if Chem.MolFromSmiles(test) is None:
                        col_names = parts
                        _info(f"Line {lineno}: skipping header row ({parts[smi_idx]!r})")
                        header_consumed = True
                        continue
                    header_consumed = True

            # Repeated column-header row (e.g. from `head`-concatenated files):
            # silently skip if the SMILES field is a known header keyword.
            if smi_idx < len(parts) and parts[smi_idx].lower() in _HEADER_KEYWORDS:
                _info(f"Line {lineno}: skipping column header ({parts[smi_idx]!r})")
                continue

            # ── Auto name column (determined once from first data row) ─────
            if auto_name and nam_idx is None:
                nam_idx = next((i for i in range(len(parts)) if i != smi_idx), None)
                auto_name = False

            # ── Parse SMILES ───────────────────────────────────────────────
            if smi_idx >= len(parts):
                _warn(f"Line {lineno}: fewer columns than expected — skipped")
                continue
            mol = Chem.MolFromSmiles(parts[smi_idx])
            if mol is None:
                _warn(f"Line {lineno}: could not parse SMILES '{parts[smi_idx]}' — skipped")
                continue

            # ── Name ──────────────────────────────────────────────────────
            name = (parts[nam_idx] if nam_idx is not None and nam_idx < len(parts) else '') \
                   or f"mol_{lineno}"
            mol.SetProp('_Name', name)

            # ── Extra columns ─────────────────────────────────────────────
            extra_names: List[str] = []
            for i, val in enumerate(parts):
                if i == smi_idx or i == nam_idx:
                    continue
                cname = col_names[i] if i < len(col_names) else f"col_{i + 1}"
                mol.SetProp(cname, val)
                extra_names.append(cname)
            if extra_names:
                mol.SetProp('_extra_col_names', '\t'.join(extra_names))

            # ── Populate output header once (requires nam_idx to be resolved) ──
            if not out_header_written and out_header is not None and col_names:
                hdr = [col_names[smi_idx] if smi_idx < len(col_names) else 'SMILES']
                if nam_idx is not None and nam_idx < len(col_names):
                    hdr.append(col_names[nam_idx])
                hdr.extend(extra_names)
                out_header.extend(hdr)
                out_header_written = True

            yield mol, name


def _iter_sdf(path: Path) -> Iterator[Tuple[Chem.Mol, str]]:
    """Yield (mol, name) from an SDF file."""
    RDLogger.DisableLog('rdApp.*')
    sup = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    for idx, mol in enumerate(sup):
        if mol is None:
            _warn(f"SDF entry {idx + 1}: could not parse — skipped")
            continue
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"mol_{idx + 1}"
        yield mol, name


def _iter_molecules(path: Path, smiles_col=None, name_col=None,
                    out_header: Optional[List[str]] = None) -> Iterator[Tuple[Chem.Mol, str]]:
    """Auto-detect format from extension and yield (mol, name) pairs."""
    if path.suffix.lower() == '.sdf':
        yield from _iter_sdf(path)
    else:
        yield from _iter_smiles(path, smiles_col=smiles_col, name_col=name_col,
                                out_header=out_header)


def _write_smiles(mol: Chem.Mol, name: str, fh) -> None:
    smi = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    extra = ''
    if mol.HasProp('_extra_col_names'):
        col_names = mol.GetProp('_extra_col_names').split('\t')
        vals = [mol.GetProp(c) if mol.HasProp(c) else '' for c in col_names]
        extra = '\t' + '\t'.join(vals)
    fh.write(f"{smi}\t{name}{extra}\n")


def _make_writer(out_path: Path, header: Optional[List[str]] = None):
    """Return (writer_fn, close_fn) for the given output path.

    header : if provided and the output is a SMILES file, written as the first line.
    """
    if out_path.suffix.lower() == '.sdf':
        w = Chem.SDWriter(str(out_path))
        def _write(mol, name):
            mol.SetProp('_Name', name)
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)
            out_mol = Chem.RWMol(mol)
            if out_mol.HasProp('_extra_col_names'):
                out_mol.ClearProp('_extra_col_names')
            w.write(out_mol.GetMol())
        def _close(): w.close()
    else:
        fh = out_path.open('w', encoding='utf-8')
        if header:
            fh.write('\t'.join(header) + '\n')
        def _write(mol, name): _write_smiles(mol, name, fh)
        def _close():          fh.close()
    return _write, _close


# ─── Preprocessing ────────────────────────────────────────────────────────────

def _largest_fragment(mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
    """Return (largest_fragment, was_stripped).

    Splits on disconnected fragments and returns the one with the most heavy
    atoms.  If the molecule is already a single fragment, returns it unchanged.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) == 1:
        return mol, False
    largest = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    # GetMolFrags strips SD properties — copy them from the original mol
    for prop in mol.GetPropNames():
        largest.SetProp(prop, mol.GetProp(prop))
    return largest, True


_uncharger = None  # lazy-initialised so it's only created when --neutralize is used
_tautomer_canon = None  # lazy-initialised so it's only created when --tautomer-canon is used

def _neutralize_mol(mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
    """Neutralize formal charges using RDKit's Uncharger.

    Returns (neutralized_mol, was_changed).  Quaternary nitrogens and other
    centres that cannot be neutralized without removing atoms are left intact.
    """
    global _uncharger
    if _uncharger is None:
        _uncharger = rdMolStandardize.Uncharger()
    uncharged = _uncharger.uncharge(mol)
    changed = Chem.MolToSmiles(uncharged) != Chem.MolToSmiles(mol)
    # Uncharger creates a new mol object — copy SD properties from the original
    if changed:
        for prop in mol.GetPropNames():
            uncharged.SetProp(prop, mol.GetProp(prop))
    return uncharged, changed


def _canonicalize_tautomer(mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
    """Pick a canonical tautomer using RDKit's TautomerEnumerator.

    Returns (canonical_mol, was_changed).  Collapses e.g. 1H- vs 2H-tetrazole,
    keto/enol, amide/imidic-acid forms so that tautomer pairs share a single
    canonical SMILES for downstream deduplication.
    """
    global _tautomer_canon
    if _tautomer_canon is None:
        _tautomer_canon = rdMolStandardize.TautomerEnumerator()
    # TautomerEnumerator mis-handles explicit-H input (e.g. emits [H]O=C(O)…
    # for a protonated carboxylic acid), producing mols with invalid valence
    # that crash later in QED/Descriptors.  Work on the implicit-H form, then
    # fall back to the original mol if canonicalization still fails sanitize.
    try:
        src = Chem.RemoveHs(mol)
        canon = _tautomer_canon.Canonicalize(src)
        Chem.SanitizeMol(canon)
        canon_smi = Chem.MolToSmiles(canon)
    except Exception:
        return mol, False
    changed = canon_smi != Chem.MolToSmiles(Chem.RemoveHs(mol))
    if changed:
        for prop in mol.GetPropNames():
            canon.SetProp(prop, mol.GetProp(prop))
    return canon, changed


def _preprocess(molecules: Iterator[Tuple[Chem.Mol, str]],
                do_strip: bool,
                do_neutralize: bool,
                do_tautomer: bool,
                do_unique: bool,
                ) -> Tuple[List[Tuple[Chem.Mol, str]], int, int, int, int]:
    """Apply salt stripping, neutralization, tautomer canonicalization, and/or
    deduplication (in that order).

    Returns (processed_list, n_stripped, n_neutralized, n_tautomer, n_duplicates).
    n_stripped    = molecules where at least one fragment was removed
    n_neutralized = molecules where at least one formal charge was neutralized
    n_tautomer    = molecules whose tautomer was changed to the canonical form
    n_duplicates  = extra entries removed during deduplication
    """
    n_stripped = 0
    n_neutralized = 0
    n_tautomer = 0
    n_duplicates = 0

    # seen maps canonical SMILES → (mol, name) for the representative entry
    seen: Dict[str, Tuple[Chem.Mol, str]] = {}
    result: List[Tuple[Chem.Mol, str]] = []

    if _TQDM:
        molecules = tqdm(molecules, desc="  Reading & preprocessing",
                         unit=" mol", file=sys.stderr,
                         disable=not sys.stderr.isatty(), leave=False)

    for mol, name in molecules:
        # 1. Salt stripping
        if do_strip:
            mol, stripped = _largest_fragment(mol)
            if stripped:
                n_stripped += 1

        # 2. Neutralization
        if do_neutralize:
            mol, changed = _neutralize_mol(mol)
            if changed:
                n_neutralized += 1

        # 3. Tautomer canonicalization
        if do_tautomer:
            mol, changed = _canonicalize_tautomer(mol)
            if changed:
                n_tautomer += 1

        if do_unique:
            canon = Chem.MolToSmiles(mol, isomericSmiles=True)
            if canon in seen:
                # Keep the entry with the lexicographically smallest identifier
                _, existing_name = seen[canon]
                if name < existing_name:
                    seen[canon] = (mol, name)
                n_duplicates += 1
            else:
                seen[canon] = (mol, name)
        else:
            result.append((mol, name))

    if do_unique:
        result = list(seen.values())

    return result, n_stripped, n_neutralized, n_tautomer, n_duplicates


# ─── Range parsing ────────────────────────────────────────────────────────────

def _parse_range(s: str, name: str) -> Tuple[Optional[float], Optional[float]]:
    """Parse a range string into (lo, hi), either of which may be None (unbounded).

    Accepted formats:
      '200:500'  →  200 ≤ x ≤ 500
      '200:'     →  x ≥ 200
      ':500'     →  x ≤ 500
      '300'      →  x == 300  (exact value, lo == hi)
    """
    if ':' in s:
        lo_s, hi_s = s.split(':', 1)
        try:
            lo = float(lo_s) if lo_s.strip() else None
            hi = float(hi_s) if hi_s.strip() else None
        except ValueError:
            _fatal(f"--{name}: could not parse range '{s}' — expected format MIN:MAX")
    else:
        try:
            lo = hi = float(s)
        except ValueError:
            _fatal(f"--{name}: could not parse value '{s}'")
    return lo, hi


# ─── Filter loaders ───────────────────────────────────────────────────────────

def _load_pains(path: Path) -> List[Tuple[str, Chem.Mol]]:
    """Load PAINS patterns from a tab-separated file (name\\tSMARTS).

    Returns a list of (name, query_mol) tuples.  Patterns that cannot be
    compiled are warned and skipped.
    """
    patterns = []
    with path.open(encoding='utf-8', errors='replace') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                _warn(f"PAINS line {lineno}: expected name\\tSMARTS — skipped")
                continue
            name, smarts = parts[0].strip(), parts[1].strip()
            q = Chem.MolFromSmarts(smarts)
            if q is None:
                _warn(f"PAINS line {lineno}: could not compile SMARTS '{smarts}' — skipped")
                continue
            patterns.append((name, q))
    return patterns


def _load_reos(path: Path) -> List[Tuple[str, int, str, Chem.Mol]]:
    """Load REOS rules from a tab-separated file (SMARTS\\tmax_count\\tdescription).

    max_count = 0  → reject if any match found
    max_count = N  → reject if more than N matches found

    Returns a list of (smarts_str, max_count, description, query_mol).
    """
    rules = []
    with path.open(encoding='utf-8', errors='replace') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                _warn(f"REOS line {lineno}: expected SMARTS\\tmax_count\\tdescription — skipped")
                continue
            smarts = parts[0].strip()
            desc   = parts[2].strip().strip('"')
            try:
                max_count = int(parts[1].strip())
            except ValueError:
                _warn(f"REOS line {lineno}: non-integer max_count — skipped")
                continue
            q = Chem.MolFromSmarts(smarts)
            if q is None:
                _warn(f"REOS line {lineno}: could not compile SMARTS '{smarts}' — skipped")
                continue
            rules.append((smarts, max_count, desc, q))
    return rules


# ─── Filters ──────────────────────────────────────────────────────────────────
# Each filter is a callable:  filter(mol, args) -> Optional[str]
#   returns None  if the molecule PASSES
#   returns a str (reason) if the molecule FAILS
#
# Add new filters here as named functions, then register them in
# _build_filter_pipeline() below.
# ──────────────────────────────────────────────────────────────────────────────

def _make_pains_filter(patterns: List[Tuple[str, Chem.Mol]]):
    """Return a PAINS filter function closed over the compiled patterns."""
    def _filt(mol, _args):
        for name, q in patterns:
            if mol.HasSubstructMatch(q):
                return f"PAINS: {name}"
        return None
    return _filt


def _make_reos_filter(rules: List[Tuple[str, int, str, Chem.Mol]]):
    """Return a REOS filter function closed over the compiled rules."""
    def _filt(mol, _args):
        for _smarts, max_count, desc, q in rules:
            n = len(mol.GetSubstructMatches(q))
            if n > max_count:
                return f"REOS: {desc}"
        return None
    return _filt


def _make_custom_filter(rules: List[Tuple[str, int, str, Chem.Mol]]):
    """Return a custom filter function; reuses the REOS rule structure."""
    def _filt(mol, _args):
        for _smarts, max_count, desc, q in rules:
            n = len(mol.GetSubstructMatches(q))
            if n > max_count:
                return f"Custom: {desc}"
        return None
    return _filt


def _make_mw_filter(lo: Optional[float], hi: Optional[float]):
    """Return a molecular-weight filter (average MW via RDKit Descriptors)."""
    def _filt(mol, _args):
        mw = Descriptors.MolWt(mol)
        if lo is not None and mw < lo:
            return f"MW {mw:.1f} < {lo}"
        if hi is not None and mw > hi:
            return f"MW {mw:.1f} > {hi}"
        return None
    return _filt


def _make_logp_filter(lo: Optional[float], hi: Optional[float]):
    """Return a Wildman-Crippen logP filter."""
    def _filt(mol, _args):
        lp = Crippen.MolLogP(mol)
        if lo is not None and lp < lo:
            return f"LogP {lp:.2f} < {lo}"
        if hi is not None and lp > hi:
            return f"LogP {lp:.2f} > {hi}"
        return None
    return _filt


def _make_lipinski_filter(strict: bool = False):
    """Return a Lipinski Rule-of-Five filter.

    Rejects molecules that violate more than one of:
      MW  ≤ 500
      HBD ≤ 5   (H-bond donors,    NH + OH)
      HBA ≤ 10  (H-bond acceptors, N + O)
      logP ≤ 5

    strict=False (default): one violation is permitted.
    strict=True: all four rules must pass.
    """
    limit = 0 if strict else 1
    def _filt(mol, _args):
        violations = []
        mw  = Descriptors.MolWt(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        lp  = Crippen.MolLogP(mol)
        if mw  > 500: violations.append(f"MW {mw:.1f}>500")
        if hbd > 5:   violations.append(f"HBD {hbd}>5")
        if hba > 10:  violations.append(f"HBA {hba}>10")
        if lp  > 5:   violations.append(f"logP {lp:.2f}>5")
        if len(violations) > limit:
            return "Lipinski: " + ", ".join(violations)
        return None
    return _filt


def _make_ro3_filter():
    """Return an Astex Rule-of-Three filter for fragment screening.

    Rejects molecules that violate any of:
      MW       ≤ 300
      logP     ≤ 3
      HBD      ≤ 3  (H-bond donors)
      HBA      ≤ 3  (H-bond acceptors)
      RotBonds ≤ 3  (rotatable bonds)
    """
    def _filt(mol, _args):
        violations = []
        mw   = Descriptors.MolWt(mol)
        lp   = Crippen.MolLogP(mol)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if mw   > 300: violations.append(f"MW {mw:.1f}>300")
        if lp   > 3:   violations.append(f"logP {lp:.2f}>3")
        if hbd  > 3:   violations.append(f"HBD {hbd}>3")
        if hba  > 3:   violations.append(f"HBA {hba}>3")
        if rotb > 3:   violations.append(f"RotBonds {rotb}>3")
        if violations:
            return "Ro3: " + ", ".join(violations)
        return None
    return _filt


def _make_qed_filter(lo: Optional[float], hi: Optional[float]):
    """Return a QED (Quantitative Estimate of Drug-likeness) filter.

    QED ranges from 0 (least drug-like) to 1 (most drug-like).
    Typical drug-like threshold: QED ≥ 0.5.
    """
    def _filt(mol, _args):
        score = QED.qed(mol)
        if lo is not None and score < lo:
            return f"QED {score:.3f} < {lo}"
        if hi is not None and score > hi:
            return f"QED {score:.3f} > {hi}"
        return None
    return _filt


def _make_chiral_filter(lo: Optional[float], hi: Optional[float]):
    """Return a chiral-centre count filter (specified + unspecified stereocenters).

    Uses FindMolChiralCenters which handles its own stereo perception — robust
    after the ToBinary roundtrip used to ship mols to worker processes
    (CalcNumAtomStereoCenters raises post-roundtrip because the CIP-perception
    flag is not pickled).
    """
    def _filt(mol, _args):
        n = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        if lo is not None and n < lo:
            return f"Chiral {n} < {int(lo)}"
        if hi is not None and n > hi:
            return f"Chiral {n} > {int(hi)}"
        return None
    return _filt


def _make_tpsa_filter(lo: Optional[float], hi: Optional[float]):
    """Return a topological polar surface area filter (Å²)."""
    def _filt(mol, _args):
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        if lo is not None and tpsa < lo:
            return f"TPSA {tpsa:.1f} < {lo}"
        if hi is not None and tpsa > hi:
            return f"TPSA {tpsa:.1f} > {hi}"
        return None
    return _filt


def _make_hba_filter(lo: Optional[float], hi: Optional[float]):
    """Return an H-bond acceptor count filter."""
    def _filt(mol, _args):
        n = rdMolDescriptors.CalcNumHBA(mol)
        if lo is not None and n < lo:
            return f"HBA {n} < {int(lo)}"
        if hi is not None and n > hi:
            return f"HBA {n} > {int(hi)}"
        return None
    return _filt


def _make_hbd_filter(lo: Optional[float], hi: Optional[float]):
    """Return an H-bond donor count filter."""
    def _filt(mol, _args):
        n = rdMolDescriptors.CalcNumHBD(mol)
        if lo is not None and n < lo:
            return f"HBD {n} < {int(lo)}"
        if hi is not None and n > hi:
            return f"HBD {n} > {int(hi)}"
        return None
    return _filt


def _make_rb_filter(lo: Optional[float], hi: Optional[float]):
    """Return a rotatable bond count filter."""
    def _filt(mol, _args):
        n = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if lo is not None and n < lo:
            return f"RotBonds {n} < {int(lo)}"
        if hi is not None and n > hi:
            return f"RotBonds {n} > {int(hi)}"
        return None
    return _filt


def _make_ha_filter(lo: Optional[float], hi: Optional[float]):
    """Return a heavy-atom count filter."""
    def _filt(mol, _args):
        n = mol.GetNumHeavyAtoms()
        if lo is not None and n < lo:
            return f"HA {n} < {int(lo)}"
        if hi is not None and n > hi:
            return f"HA {n} > {int(hi)}"
        return None
    return _filt


def _mol_props(mol: Chem.Mol) -> dict:
    """Return a dict of calculated properties for a molecule.

    'Chiral' counts both specified and unspecified stereocentres via
    FindMolChiralCenters — CalcNumAtomStereoCenters raises after the ToBinary
    roundtrip used for worker transfer (CIP-perception state is not pickled).
    """
    return {
        'MW':     Descriptors.MolWt(mol),
        'LogP':   Crippen.MolLogP(mol),
        'HBA':    rdMolDescriptors.CalcNumHBA(mol),
        'HBD':    rdMolDescriptors.CalcNumHBD(mol),
        'RB':     rdMolDescriptors.CalcNumRotatableBonds(mol),
        'TPSA':   rdMolDescriptors.CalcTPSA(mol),
        'QED':    QED.qed(mol),
        'Chiral': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        'HA':     float(mol.GetNumHeavyAtoms()),
    }


# ─── Outlier detection ────────────────────────────────────────────────────────

_OUTLIER_PROP_KEY = {
    'mw':     'MW',
    'logp':   'LogP',
    'hba':    'HBA',
    'hbd':    'HBD',
    'rb':     'RB',
    'tpsa':   'TPSA',
    'qed':    'QED',
    'chiral': 'Chiral',
    'ha':     'HA',
}
_VALID_OUTLIER_PROPS = tuple(_OUTLIER_PROP_KEY.keys())


def _compute_outlier_bounds(molecules: List[Tuple[Chem.Mol, str]],
                             props: List[str],
                             k: float) -> Dict[str, Tuple[float, float]]:
    """Compute Tukey IQR fences for each property across the molecule set.

    Returns {prop: (lower_fence, upper_fence)} where
    lower_fence = Q1 - k*IQR,  upper_fence = Q3 + k*IQR.
    """
    prop_vals: Dict[str, List[float]] = {p: [] for p in props}
    for mol, _ in molecules:
        p_dict = _mol_props(mol)
        for p in props:
            prop_vals[p].append(p_dict[p])

    bounds: Dict[str, Tuple[float, float]] = {}
    for p, vals in prop_vals.items():
        if len(vals) < 4:
            _warn(f"Outlier: too few molecules to compute IQR for {p} — skipped")
            continue
        qs = statistics.quantiles(vals, n=4)   # [Q1, Q2, Q3]
        q1, q3 = qs[0], qs[2]
        iqr = q3 - q1
        bounds[p] = (q1 - k * iqr, q3 + k * iqr)
    return bounds


def _make_outlier_filter(bounds: Dict[str, Tuple[float, float]]):
    """Return an IQR outlier filter closed over pre-computed per-property fences."""
    def _filt(mol, _args):
        props = _mol_props(mol)
        for p, (lo, hi) in bounds.items():
            val = props[p]
            if val < lo:
                return f"Outlier {p} {val:.2f} < fence {lo:.2f}"
            if val > hi:
                return f"Outlier {p} {val:.2f} > fence {hi:.2f}"
        return None
    return _filt


# Shared pipeline inherited by worker processes via fork (Linux default).
# Set in main() before the Pool is created.
_shared_pipeline: list = []


_PICKLE_PROPS = (Chem.PropertyPickleOptions.MolProps |
                 Chem.PropertyPickleOptions.PrivateProps)


def _pool_worker(batch_binary: list) -> list:
    """Multiprocessing worker. Inherits _shared_pipeline via fork — no pickling.

    batch_binary : list of (mol_bytes, name)
    Returns      : list of (mol_bytes_or_None, name, reason_or_None, props)

    Mol properties (extra input columns etc.) are preserved across the
    ToBinary / Chem.Mol round-trip via PropertyPickleOptions.
    """
    results = []
    for mol_b, name in batch_binary:
        mol = Chem.Mol(mol_b)
        reason = None
        for _label, filt_fn in _shared_pipeline:
            reason = filt_fn(mol, None)
            if reason:
                break
        props = _mol_props(mol)
        results.append((mol.ToBinary(_PICKLE_PROPS) if reason is None else None, name, reason, props))
    return results


def _build_filter_pipeline(pains_patterns=None, reos_rules=None,
                           custom_rules=None, mw_range=None, logp_range=None,
                           lipinski=False, lipinski_strict=False, ro3=False,
                           qed_range=None, hba_range=None, hbd_range=None,
                           rb_range=None, tpsa_range=None,
                           chiral_range=None, ha_range=None) -> list:
    """Return an ordered list of (label, filter_fn) tuples to apply."""
    pipeline = []
    if pains_patterns is not None:
        pipeline.append(('PAINS',   _make_pains_filter(pains_patterns)))
    if reos_rules is not None:
        pipeline.append(('REOS',    _make_reos_filter(reos_rules)))
    if custom_rules is not None:
        pipeline.append(('Custom',  _make_custom_filter(custom_rules)))
    if mw_range is not None:
        pipeline.append(('MW',      _make_mw_filter(*mw_range)))
    if logp_range is not None:
        pipeline.append(('LogP',    _make_logp_filter(*logp_range)))
    if lipinski:
        pipeline.append(('Lipinski', _make_lipinski_filter(strict=lipinski_strict)))
    if ro3:
        pipeline.append(('Ro3',      _make_ro3_filter()))
    if qed_range is not None:
        pipeline.append(('QED',      _make_qed_filter(*qed_range)))
    if hba_range is not None:
        pipeline.append(('HBA',      _make_hba_filter(*hba_range)))
    if hbd_range is not None:
        pipeline.append(('HBD',      _make_hbd_filter(*hbd_range)))
    if rb_range is not None:
        pipeline.append(('RotBonds', _make_rb_filter(*rb_range)))
    if tpsa_range is not None:
        pipeline.append(('TPSA',     _make_tpsa_filter(*tpsa_range)))
    if chiral_range is not None:
        pipeline.append(('Chiral',   _make_chiral_filter(*chiral_range)))
    if ha_range is not None:
        pipeline.append(('HA',       _make_ha_filter(*ha_range)))
    return pipeline


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    if not _RDKIT:
        _fatal("RDKit is required:  conda install -c conda-forge rdkit")

    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_HelpFmt,
    )

    io = p.add_argument_group('Input / Output')
    io.add_argument('input', metavar='FILE', nargs='?', default=None,
                    help='Input SMILES (.smi/.csv/.tsv) or SDF file')
    io.add_argument('-i', dest='input_flag', metavar='FILE', default=None,
                    help='Input file (alternative to positional argument)')
    io.add_argument('-o', '--output', metavar='FILE',
                    help='Output file (default: <input>_filtered.<ext>)')
    io.add_argument('-j', '--jobs', metavar='N', type=int,
                    default=os.cpu_count(),
                    help=f'Worker processes for filtering and property calculation '
                         f'(default: {os.cpu_count()} — all available CPUs). '
                         f'Use --jobs 1 to disable parallelism.')
    io.add_argument('--smiles-col', metavar='COL', default=None,
                    help='Column containing SMILES: 1-based integer index or column '
                         'header name (default: auto-detect)')
    io.add_argument('--name-col', metavar='COL', default=None,
                    help='Column containing molecule name/ID: 1-based integer index '
                         'or column header name (default: auto = first non-SMILES column)')

    pre = p.add_argument_group('Preprocessing')
    pre.add_argument('--strip', dest='strip', action='store_true', default=True,
                     help='Strip salts and small fragments; keep the largest '
                          'fragment by heavy-atom count (default: on)')
    pre.add_argument('--no-strip', dest='strip', action='store_false',
                     help='Disable salt stripping')
    pre.add_argument('--neutralize', dest='neutralize', action='store_true', default=True,
                     help='Neutralize formal charges where possible '
                          '(e.g. carboxylate → carboxylic acid, ammonium → amine). '
                          'Quaternary N and other centres that cannot be neutralized '
                          'without atom removal are left intact (default: on)')
    pre.add_argument('--no-neutralize', dest='neutralize', action='store_false',
                     help='Disable neutralization')
    pre.add_argument('--tautomer-canon', dest='tautomer', action='store_true', default=True,
                     help='Pick a canonical tautomer (e.g. collapse 1H-/2H-tetrazole, '
                          'keto/enol, amide/imidic-acid) so tautomer pairs share a '
                          'single canonical SMILES for dedup (default: on; slower per '
                          'molecule than the other preprocessing steps)')
    pre.add_argument('--no-tautomer-canon', dest='tautomer', action='store_false',
                     help='Disable tautomer canonicalization')
    pre.add_argument('--unique', dest='unique', action='store_true', default=True,
                     help='Deduplicate on canonical SMILES; for duplicates '
                          'keep the entry with the lexicographically smallest '
                          'identifier (default: on)')
    pre.add_argument('--no-unique', dest='unique', action='store_false',
                     help='Disable deduplication')
    filt = p.add_argument_group('Filters')
    filt.add_argument('--pains', action='store_true',
                      help=f'Reject PAINS (Pan-Assay INterference compoundS); '
                           f'uses {_DEFAULT_PAINS.name} from the script directory')
    filt.add_argument('--pains-file', metavar='FILE', default=None,
                      help='Custom PAINS file (implies --pains; '
                           'format: name<TAB>SMARTS)')
    filt.add_argument('--reos', action='store_true',
                      help=f'Reject REOS (Rapid Elimination Of Swill) unwanted '
                           f'functional groups; uses {_DEFAULT_REOS.name}')
    filt.add_argument('--reos-file', metavar='FILE', default=None,
                      help='Custom REOS file (implies --reos; '
                           'format: SMARTS<TAB>max_count<TAB>description)')
    filt.add_argument('--custom', action='store_true',
                      help=f'Apply custom SMARTS filters from '
                           f'{_DEFAULT_CUSTOM.name} in the script directory '
                           f'(same format as REOS: SMARTS<TAB>max_count<TAB>description)')
    filt.add_argument('--custom-file', metavar='FILE', default=None,
                      help='Custom filter file (implies --custom; '
                           'format: SMARTS<TAB>max_count<TAB>description)')

    prop = p.add_argument_group('Property ranges')
    prop.add_argument('--mw', metavar='RANGE', default=None,
                      help='Molecular weight range (average MW).  '
                           'Format: MIN:MAX, MIN:, :MAX, or exact value.  '
                           'E.g. --mw 150:500  --mw :600  --mw 250:')
    prop.add_argument('--logp', metavar='RANGE', default=None,
                      help='Wildman-Crippen logP range.  '
                           'Format: MIN:MAX, MIN:, :MAX, or exact value.  '
                           'E.g. --logp -2:5  --logp :4.5')
    prop.add_argument('--ro5', action='store_true',
                      help='Reject molecules that violate more than one '
                           'Lipinski Rule of Five (MW≤500, HBD≤5, HBA≤10, '
                           'logP≤5).  One violation is permitted by default.')
    prop.add_argument('--ro5-strict', action='store_true',
                      help='Require all four Lipinski rules to pass '
                           '(zero violations allowed; implies --ro5)')
    prop.add_argument('--qed', metavar='RANGE', default=None,
                      help='Quantitative Estimate of Drug-likeness (0–1, '
                           'higher = more drug-like).  '
                           'Format: MIN:MAX, MIN:, :MAX, or exact value.  '
                           'E.g. --qed 0.5:  --qed 0.4:0.9')
    prop.add_argument('--chiral', metavar='RANGE', default=None,
                      help='Chiral centre count range (specified + unspecified).  '
                           'E.g. --chiral :3  --chiral 0:2')
    prop.add_argument('--tpsa', metavar='RANGE', default=None,
                      help='Topological polar surface area range (Å²).  '
                           'E.g. --tpsa :140  --tpsa 40:130')
    prop.add_argument('--hba', metavar='RANGE', default=None,
                      help='H-bond acceptor count range.  '
                           'E.g. --hba :10  --hba 1:8')
    prop.add_argument('--hbd', metavar='RANGE', default=None,
                      help='H-bond donor count range.  '
                           'E.g. --hbd :5  --hbd 1:3')
    prop.add_argument('--rb', metavar='RANGE', default=None,
                      help='Rotatable bond count range.  '
                           'E.g. --rb :10  --rb 2:8')
    prop.add_argument('--ro3', action='store_true',
                      help='Astex Rule of Three for fragment screening: '
                           'MW≤300, logP≤3, HBD≤3, HBA≤3, RotBonds≤3 '
                           '(all rules must pass)')
    prop.add_argument('--ha', metavar='RANGE', default=None,
                      help='Heavy-atom count range.  '
                           'E.g. --ha 10:35  --ha :40  --ha 5:')

    out = p.add_argument_group('Outlier detection')
    out.add_argument('--outlier', action='store_true',
                     help='Reject statistical outliers using Tukey IQR fences '
                          '(Q1 - k*IQR, Q3 + k*IQR).  Bounds are computed '
                          'from the preprocessed molecule set, so the fence '
                          'values reflect the actual input distribution.')
    out.add_argument('--outlier-iqr', metavar='K', type=float, default=1.5,
                     help='IQR multiplier k (default: 1.5 = standard box-plot '
                          'whiskers; use 3.0 to remove only extreme outliers)')
    out.add_argument('--outlier-props', metavar='PROPS', default='mw,logp,rb',
                     help='Comma-separated properties to check.  '
                          f'Valid: {", ".join(_VALID_OUTLIER_PROPS)}  '
                          '(default: mw,logp,rb)')

    if _ARGCOMPLETE:
        argcomplete.autocomplete(p)

    args = p.parse_args()

    smiles_col = _parse_col_spec(args.smiles_col) if args.smiles_col is not None else None
    name_col   = _parse_col_spec(args.name_col) if args.name_col is not None else None

    t0 = time.time()

    input_file = args.input or args.input_flag
    if not input_file:
        p.error("input file is required (positional FILE or -i FILE)")
    in_path = Path(input_file).resolve()
    if not in_path.exists():
        _fatal(f"Input file not found: {in_path}")

    if args.output:
        out_path = Path(args.output).resolve()
    else:
        out_path = in_path.with_name(in_path.stem + '_filtered' + in_path.suffix)

    log_path = out_path.with_suffix('.log')
    _log_fh  = log_path.open('w', encoding='utf-8')
    sys.stdout = _Tee(_log_fh, sys.__stdout__)

    # ── Config ────────────────────────────────────────────────────────────────
    bar = '─' * 60
    print(bar)
    print("  ligfilter.py — Ligand structural filter")
    print(bar)
    _info(f"Input:         {in_path.name}")
    _info(f"Output:        {out_path.name}")
    _info(f"Strip salts:   {'yes' if args.strip else 'no (disabled)'}")
    _info(f"Neutralize:    {'yes' if args.neutralize else 'no (disabled)'}")
    _info(f"Tautomer canon:{'yes' if args.tautomer else 'no (disabled)'}")
    _info(f"Deduplicate:   {'yes' if args.unique else 'no (disabled)'}")
    _info(f"Output format: {out_path.suffix.lower()[1:].upper()}")
    _info(f"Workers:       {args.jobs}")
    _info(f"SMILES col:    {args.smiles_col if args.smiles_col is not None else 'auto'}")
    _info(f"Name col:      {args.name_col if args.name_col is not None else 'auto'}")
    _info(f"MW range:      {args.mw if args.mw else 'any'}")
    _info(f"LogP range:    {args.logp if args.logp else 'any'}")
    _info(f"QED range:     {args.qed if args.qed else 'any'}")
    _info(f"HBA range:     {args.hba if args.hba else 'any'}")
    _info(f"HBD range:     {args.hbd if args.hbd else 'any'}")
    _info(f"RotBonds range:{args.rb   if args.rb   else 'any'}")
    _info(f"TPSA range:    {args.tpsa   if args.tpsa   else 'any'}")
    _info(f"Chiral range:  {args.chiral if args.chiral else 'any'}")
    _info(f"HA range:      {args.ha if args.ha else 'any'}")
    lip_mode = ('strict' if args.ro5_strict else 'on (1 violation allowed)') if (args.ro5 or args.ro5_strict) else 'off'
    _info(f"Lipinski Ro5:  {lip_mode}")
    _info(f"Rule of Three: {'on' if args.ro3 else 'off'}")
    if args.outlier:
        _info(f"Outlier filter: IQR k={args.outlier_iqr}, props={args.outlier_props}")
    print(bar)

    # ── Load filter data ──────────────────────────────────────────────────────
    pains_patterns = None
    reos_rules     = None

    if args.pains or args.pains_file:
        pains_path = Path(args.pains_file).resolve() if args.pains_file else _DEFAULT_PAINS
        if not pains_path.exists():
            _fatal(f"PAINS file not found: {pains_path}")
        pains_patterns = _load_pains(pains_path)
        _ok(f"PAINS:  {len(pains_patterns)} patterns loaded  ({pains_path.name})")

    if args.reos or args.reos_file:
        reos_path = Path(args.reos_file).resolve() if args.reos_file else _DEFAULT_REOS
        if not reos_path.exists():
            _fatal(f"REOS file not found: {reos_path}")
        reos_rules = _load_reos(reos_path)
        _ok(f"REOS:   {len(reos_rules)} rules loaded  ({reos_path.name})")

    custom_rules = None
    if args.custom or args.custom_file:
        custom_path = Path(args.custom_file).resolve() if args.custom_file else _DEFAULT_CUSTOM
        if not custom_path.exists():
            _fatal(f"Custom filter file not found: {custom_path}")
        custom_rules = _load_reos(custom_path)   # same format as REOS
        _ok(f"Custom: {len(custom_rules)} rules loaded  ({custom_path.name})")

    mw_range   = _parse_range(args.mw,   'mw')   if args.mw   else None
    logp_range = _parse_range(args.logp, 'logp') if args.logp else None
    qed_range  = _parse_range(args.qed,  'qed')  if args.qed  else None
    hba_range  = _parse_range(args.hba,  'hba')  if args.hba  else None
    hbd_range  = _parse_range(args.hbd,  'hbd')  if args.hbd  else None
    rb_range   = _parse_range(args.rb,   'rb')   if args.rb   else None
    tpsa_range   = _parse_range(args.tpsa,   'tpsa')   if args.tpsa   else None
    chiral_range = _parse_range(args.chiral, 'chiral') if args.chiral else None
    ha_range     = _parse_range(args.ha,     'ha')     if args.ha     else None

    pipeline = _build_filter_pipeline(pains_patterns=pains_patterns,
                                      reos_rules=reos_rules,
                                      custom_rules=custom_rules,
                                      mw_range=mw_range,
                                      logp_range=logp_range,
                                      lipinski=args.ro5 or args.ro5_strict,
                                      lipinski_strict=args.ro5_strict,
                                      ro3=args.ro3,
                                      qed_range=qed_range,
                                      hba_range=hba_range,
                                      hbd_range=hbd_range,
                                      rb_range=rb_range,
                                      tpsa_range=tpsa_range,
                                      chiral_range=chiral_range,
                                      ha_range=ha_range)

    # ── Read & preprocess ─────────────────────────────────────────────────────
    out_header: List[str] = []
    raw_stream = _iter_molecules(in_path, smiles_col=smiles_col, name_col=name_col,
                                 out_header=out_header)

    molecules, n_stripped, n_neutralized, n_tautomer, n_duplicates = _preprocess(
        raw_stream, do_strip=args.strip, do_neutralize=args.neutralize,
        do_tautomer=args.tautomer, do_unique=args.unique)
    n_read = len(molecules) + n_duplicates

    # ── Outlier filter (needs population statistics from preprocessed set) ────
    if args.outlier:
        outlier_props = [p.strip() for p in args.outlier_props.split(',')]
        invalid = [p for p in outlier_props if p not in _VALID_OUTLIER_PROPS]
        if invalid:
            _fatal(f"--outlier-props: unknown: {invalid}  valid: {list(_VALID_OUTLIER_PROPS)}")
        outlier_bounds = _compute_outlier_bounds(
            molecules, [_OUTLIER_PROP_KEY[p] for p in outlier_props], args.outlier_iqr)
        _ok(f"Outlier fences (IQR k={args.outlier_iqr}, n={len(molecules)}):")
        for p, (lo, hi) in outlier_bounds.items():
            _info(f"  {p:<8} [{lo:.3g},  {hi:.3g}]")
        pipeline.append(('Outlier', _make_outlier_filter(outlier_bounds)))

    if not pipeline:
        _warn("No preprocessing or filters enabled — all molecules will pass.")
    else:
        _info(f"Active filters ({len(pipeline)}): "
              + ', '.join(label for label, _ in pipeline))
    print(bar)

    # ── Filter ────────────────────────────────────────────────────────────────
    n_pass = n_fail = 0
    rejection_counts: dict = {}
    prop_lists:      Dict[str, list] = {}
    rej_prop_lists:  Dict[str, list] = {}

    # Share pipeline with workers via fork-inherited global (no pickling needed)
    global _shared_pipeline
    _shared_pipeline = pipeline

    # Serialise mols as bytes for inter-process transfer.
    # PropertyPickleOptions preserves all mol properties (including extra input
    # columns stored as string props) across the ToBinary / Chem.Mol round-trip.
    n_jobs = max(1, min(args.jobs, len(molecules)))
    # Aim for ~20 batches per worker so the progress bar updates frequently
    # without ballooning IPC overhead. Cap batch size to keep memory bounded.
    target_batches = max(n_jobs * 20, 50)
    batch_size = max(1, min(2000, math.ceil(len(molecules) / target_batches))) if molecules else 1
    batches_binary = [
        [(mol.ToBinary(_PICKLE_PROPS), name) for mol, name in molecules[i:i + batch_size]]
        for i in range(0, len(molecules), batch_size)
    ]

    write_mol, close_out = _make_writer(out_path, header=out_header or None)

    import multiprocessing
    _fork_available = multiprocessing.get_start_method(allow_none=True) in (None, 'fork')
    if n_jobs > 1 and not _fork_available:
        _warn(f"Multiprocessing requires 'fork' start method (not available on this platform). "
              f"Falling back to --jobs 1.")
        n_jobs = 1

    pool = Pool(processes=n_jobs) if n_jobs > 1 else None
    if pool is not None:
        result_iter = pool.imap_unordered(_pool_worker, batches_binary)
    else:
        result_iter = (_pool_worker(bb) for bb in batches_binary)

    pbar = None
    if _TQDM:
        pbar = tqdm(total=len(molecules), desc="  Filtering",
                    unit=" mol", file=sys.stderr,
                    disable=not sys.stderr.isatty(), leave=False)

    try:
        for batch_results in result_iter:
            for mol_b, name, reason, props in batch_results:
                if reason is None:
                    write_mol(Chem.Mol(mol_b), name)
                    n_pass += 1
                    for prop, val in props.items():
                        prop_lists.setdefault(prop, []).append(val)
                else:
                    rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
                    n_fail += 1
                    for prop, val in props.items():
                        rej_prop_lists.setdefault(prop, []).append(val)
            if pbar is not None:
                pbar.update(len(batch_results))
    finally:
        if pbar is not None:
            pbar.close()
        if pool is not None:
            pool.close()
            pool.join()

    close_out()

    # ── Summary ───────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    print(bar)
    print("  [SUMMARY]")
    _info(f"Molecules read:    {n_read}")
    _info(f"Salts stripped:    {n_stripped}")
    _info(f"Neutralized:       {n_neutralized}")
    _info(f"Tautomer changed:  {n_tautomer}")
    _info(f"Duplicates removed:{n_duplicates}")
    _info(f"Passed filters:    {n_pass}")
    _info(f"Rejected:          {n_fail}")
    if rejection_counts:
        _info("Rejection reasons:")
        for reason, count in sorted(rejection_counts.items(),
                                    key=lambda x: -x[1]):
            _info(f"  {count:>6}  {reason}")
    fmt = {
        'MW':     '{:8.1f}',
        'LogP':   '{:8.2f}',
        'HBA':    '{:8.1f}',
        'HBD':    '{:8.1f}',
        'RB':     '{:8.1f}',
        'TPSA':   '{:8.1f}',
        'QED':    '{:8.3f}',
        'Chiral': '{:8.1f}',
        'HA':     '{:8.0f}',
    }

    def _print_prop_table(label: str, pl: Dict[str, list]):
        print(bar)
        print(f"  [PROPERTY STATISTICS — {label}]")
        _info(f"  {'Property':<10} {'Mean':>8} {'Median':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
        _info(f"  {'-'*10} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
        for prop, vals in pl.items():
            mean   = statistics.mean(vals)
            median = statistics.median(vals)
            std    = statistics.stdev(vals) if len(vals) > 1 else 0.0
            f      = fmt[prop]
            _info(f"  {prop:<10} " +
                  f.format(mean)       + ' ' +
                  f.format(median)     + ' ' +
                  f.format(std)        + ' ' +
                  f.format(min(vals))  + ' ' +
                  f.format(max(vals)))

    if prop_lists:
        _print_prop_table('passing molecules', prop_lists)
    if rej_prop_lists:
        _print_prop_table('rejected molecules', rej_prop_lists)
    _info(f"Time:              {format_time(elapsed)}")
    _info(f"Log:               {log_path.name}")
    print(bar)

    sys.stdout = sys.__stdout__
    _log_fh.close()


if __name__ == '__main__':
    main()
