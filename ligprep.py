#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
Ligand Preparation

Uses CDPKit CONFORGE for 3D conformer generation — a knowledge-based torsion
sampler derived from CSD/PDB experimental structures that outperforms RDKit
ETKDGv3 on crystal-structure reproduction benchmarks and is ~50-100× faster.
CONFORGE handles H addition, MMFF94 minimization, energy window filtering, and
RMSD deduplication internally; no post-hoc chair filter or cregen needed.

Workflow
--------
  Input SMILES / SDF
    → salt strip (largest fragment, RDKit)
    → cxcalc tautomers -D true -n true -H <pH>
        combined tautomer × protomer distribution at target pH
        low-population states filtered by --min-population
    → stereo enumeration (RDKit EnumerateStereoisomers, unspecified centers only)
    → CONFORGE 3D conformer generation (CDPKit)
        internal: MMFF94 minimization + energy window filter + RMSD deduplication
    → dominant/states: 1 best conformer per state
    → ensemble: CONFORGE multi-conformer output (energy-sorted, RMSD-deduplicated)
    → output SDF

Output SDF properties
---------------------
  _Name                molecule name (all modes)
  tautomer_id          tautomer/protomer state index (1-based)
  tautomer_population  fractional population at target pH (0.0–1.0)
  stereo_id            stereo variant index (1-based)
  state_population     tautomer_population × stereo_weight (all modes)
  state_penalty        –log10(state_population); add to docking score (all modes)
  mmff_energy          MMFF94s energy kcal/mol (all modes)
  conf_id              conformer index within tautomer/stereo state (ensemble only)

Dependencies
------------
  Required : rdkit, cdpkit, cxcalc (ChemAxon, on PATH)
  Optional : tqdm, rich-argparse, argcomplete

SMILES input format
-------------------
  Column 1: SMILES string
  Column 2 (optional): molecule name / identifier
  Delimiter: whitespace, comma, or tab

References
----------
  CONFORGE: Seidel et al., J. Chem. Inf. Model. 2023, 63, 5549-5570.
            DOI: 10.1021/acs.jcim.3c00563

Written by Claude, 2026-04-17
"""

from __future__ import annotations

import argparse
import io
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

try:
    from rich_argparse import RawDescriptionRichHelpFormatter as _HelpFmt
except ImportError:
    _HelpFmt = argparse.RawDescriptionHelpFormatter

try:
    from tqdm import tqdm as _tqdm
    _TQDM = True
except ImportError:
    _TQDM = False
    def _tqdm(it, **kw):
        return it

try:
    import argcomplete
    from argcomplete.completers import FilesCompleter
    _ARGCOMPLETE = True
except ImportError:
    _ARGCOMPLETE = False

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from rdkit.Chem.EnumerateStereoisomers import (
        EnumerateStereoisomers,
        StereoEnumerationOptions,
    )
    RDLogger.DisableLog("rdApp.*")
    _RDKIT = True
except ImportError:
    _RDKIT = False

_CXCALC = shutil.which("cxcalc") is not None

try:
    import CDPL.Chem as _CDPLChem
    import CDPL.ConfGen as _CDPLConfGen
    _CDPKIT = True
except ImportError:
    _CDPKIT = False


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

class _Tee:
    """Write to multiple streams simultaneously (stdout + log file)."""
    def __init__(self, *streams):
        self._streams = streams
    def write(self, data: str) -> None:
        for s in self._streams:
            s.write(data)
    def flush(self) -> None:
        for s in self._streams:
            s.flush()


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m{s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m{s:02d}s"


def _strip_salts(mol: "Chem.Mol") -> "Chem.Mol":
    """Return the largest fragment (by heavy atom count)."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return frags[0] if len(frags) == 1 else max(frags, key=lambda m: m.GetNumHeavyAtoms())


# ─────────────────────────────────────────────────────────────────────────────
# Energy evaluation (MMFF94s, no re-minimization)
# ─────────────────────────────────────────────────────────────────────────────

def _eval_mmff_energies(mol_h: "Chem.Mol", cids: list) -> list:
    """Evaluate MMFF94s single-point energy for each conformer (no minimization)."""
    props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
    if props is None:
        return [0.0] * len(cids)
    energies = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=cid)
        energies.append(ff.CalcEnergy() if ff is not None else 0.0)
    return energies


# ─────────────────────────────────────────────────────────────────────────────
# CONFORGE 3D conformer generation
# ─────────────────────────────────────────────────────────────────────────────

def _embed_conforge(
    mol: "Chem.Mol",
    n_confs_max: int = 50,
    ewin: float = 15.0,
    rmsd_thr: float = 0.5,
    name: str = "",
    quiet: bool = False,
    timeout_ms: int = 300000,
    max_sampled: int = 2000,
) -> "tuple[Chem.Mol, list, list, str | None]":
    """Generate 3D conformers with CDPKit CONFORGE.

    CONFORGE handles H addition, MMFF94 minimization, energy window filtering,
    and RMSD deduplication internally.  The RDKit ↔ CDPKit bridge goes via
    temp SDF files (safest cross-toolkit interchange).

    Returns (mol_h, cids, energies, None) on success or (None, None, None, reason) on failure.
    """
    label = name or (Chem.MolToSmiles(mol) if mol else "?")

    if not _CDPKIT:
        reason = "CDPKit not available"
        if not quiet:
            print(f"[ERROR] {reason} — cannot generate conformers for {label}")
        return None, None, None, reason

    try:
        with tempfile.TemporaryDirectory(prefix="lgp3_confgen_") as tmpdir:
            # RDKit mol → SDF → CDPKit
            in_sdf = Path(tmpdir) / "input.sdf"
            w = Chem.SDWriter(str(in_sdf))
            w.write(mol); w.flush(); w.close()

            cdpl_mol = _CDPLChem.BasicMolecule()
            reader = _CDPLChem.MoleculeReader(str(in_sdf))
            if not reader.read(cdpl_mol):
                raise RuntimeError("CDPKit could not read input SDF")
            reader.close()

            # CONFORGE
            gen = _CDPLConfGen.ConformerGenerator()
            gen.settings.energyWindow              = ewin
            gen.settings.minRMSD                   = rmsd_thr
            gen.settings.maxNumOutputConformers     = n_confs_max
            gen.settings.maxNumSampledConformers    = max_sampled
            gen.settings.timeout                   = timeout_ms

            _CDPLConfGen.prepareForConformerGeneration(cdpl_mol)
            status = gen.generate(cdpl_mol)

            _ok = {_CDPLConfGen.ReturnCode.SUCCESS,
                   _CDPLConfGen.ReturnCode.TOO_MUCH_SYMMETRY}
            n_out = gen.getNumConformers()
            if status not in _ok or n_out == 0:
                raise RuntimeError(f"CONFORGE status={status}, n_conformers={n_out}")
            gen.setConformers(cdpl_mol)

            # CDPKit → SDF → RDKit (one SDF record per conformer)
            out_sdf = Path(tmpdir) / "output.sdf"
            _CDPLChem.setMDLDimensionality(cdpl_mol, 3)
            cw = _CDPLChem.MolecularGraphWriter(str(out_sdf))
            if not cw.write(cdpl_mol):
                raise RuntimeError("CDPKit SDF write failed")
            cw.close()

            suppl = Chem.SDMolSupplier(str(out_sdf), removeHs=False, sanitize=True)
            rdmols = [m for m in suppl if m is not None]
            if not rdmols:
                raise RuntimeError("RDKit could not read CDPKit SDF output")

            # Merge into single RDKit mol with multiple conformers
            mol_h = Chem.RWMol(rdmols[0])
            if len(rdmols) > 1:
                from rdkit.Chem import Conformer as _Conf
                from rdkit.Geometry import rdGeometry as _G
                n_atoms = mol_h.GetNumAtoms()
                for rm in rdmols[1:]:
                    if rm.GetNumAtoms() != n_atoms:
                        continue
                    src = rm.GetConformer()
                    c = _Conf(n_atoms)
                    for i in range(n_atoms):
                        p = src.GetAtomPosition(i)
                        c.SetAtomPosition(i, _G.Point3D(p.x, p.y, p.z))
                    mol_h.AddConformer(c, assignId=True)
            mol_h = mol_h.GetMol()

    except Exception as exc:
        reason = str(exc)
        if not quiet:
            print(f"[WARNING] CONFORGE failed for {label} ({reason})")
        return None, None, None, reason

    cids     = [c.GetId() for c in mol_h.GetConformers()]
    energies = _eval_mmff_energies(mol_h, cids)
    return mol_h, cids, energies, None


def _embed_etkdg(
    mol: "Chem.Mol",
    n_confs_max: int = 50,
    name: str = "",
    quiet: bool = False,
) -> "tuple[Chem.Mol, list, list, str | None]":
    """ETKDGv3 + MMFF94s fallback conformer generation (no CDPKit required)."""
    label = name or (Chem.MolToSmiles(mol) if mol else "?")
    try:
        mol_h = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.numThreads    = 1
        params.maxIterations = 50    # 50 attempts × ~14ms/attempt → ≤0.7s worst case
        cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs_max, params=params))
        if not cids:
            return None, None, None, "ETKDGv3: no conformers embedded"
        AllChem.MMFFOptimizeMoleculeConfs(mol_h, mmffVariant="MMFF94s", numThreads=1,
                                          maxIters=500)
        energies = _eval_mmff_energies(mol_h, cids)
        return mol_h, cids, energies, None
    except Exception as exc:
        if not quiet:
            print(f"[WARNING] ETKDGv3 failed for {label}: {exc}")
        return None, None, None, f"ETKDGv3: {exc}"


# ─────────────────────────────────────────────────────────────────────────────
# Conformer selection helpers
# ─────────────────────────────────────────────────────────────────────────────

def _best_conf(mol_h: "Chem.Mol", cids: list, energies: list) -> "tuple[Chem.Mol, float]":
    """Keep only the lowest-energy conformer."""
    best_idx = energies.index(min(energies))
    best_cid = cids[best_idx]
    for cid in cids:
        if cid != best_cid:
            mol_h.RemoveConformer(cid)
    return mol_h, energies[best_idx]


def _mol_to_sdf(mol: "Chem.Mol", conf_id: int = -1) -> str:
    """Serialize a single RDKit mol (one conformer) to SDF text."""
    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    w.write(mol, confId=conf_id)
    w.flush()
    w.close()
    return buf.getvalue()


# ─────────────────────────────────────────────────────────────────────────────
# ProcessPoolExecutor worker
# ─────────────────────────────────────────────────────────────────────────────

def _mol_job(job: dict) -> "tuple[list[str], list[dict]]":
    """Worker: stereo enum + CONFORGE for all tautomer states of one molecule.

    Returns (list_of_sdf_texts, failures) where each failure is
    {"id": state_label, "stage": "conforge", "reason": str}.
    Must be a module-level function so multiprocessing can pickle it.
    """
    mol_id        = job["mol_id"]
    no_stereo     = job["no_stereo"]
    max_stereo    = job["max_stereo"]
    nmax          = job["nmax"]
    ewin          = job["ewin"]
    rmsd_thr      = job["rmsd_thr"]
    mode          = job["mode"]
    max_confs_out = job["max_confs_out"]
    timeout_ms    = job["timeout_ms"]
    max_sampled   = job["max_sampled"]
    sdf_texts = []
    failures  = []
    try:
        for t_idx, mol_bytes, t_pop in job["states"]:
            t_mol   = Chem.Mol(mol_bytes)
            isomers = [t_mol] if no_stereo else _enumerate_stereo(t_mol, max_stereo)
            stereo_weight = 1.0 / len(isomers)

            for s_idx, iso in enumerate(isomers, start=1):
                state_pop     = t_pop * stereo_weight
                state_penalty = -math.log10(state_pop) if state_pop > 0 else 999.0
                state_label   = f"{mol_id} T{t_idx}S{s_idx}"

                mol_h, cids, energies, reason = _embed_conforge(
                    iso, n_confs_max=nmax, ewin=ewin, rmsd_thr=rmsd_thr,
                    name=state_label, quiet=True, timeout_ms=timeout_ms,
                    max_sampled=max_sampled,
                )
                embedder = "conforge"
                if mol_h is None:
                    mol_h, cids, energies, reason2 = _embed_etkdg(
                        iso, n_confs_max=1, name=state_label, quiet=True,
                    )
                    embedder = "etkdg"
                    if mol_h is None:
                        failures.append({"id": state_label, "stage": "conforge",
                                         "reason": f"{reason} | etkdg: {reason2}"})
                        continue

                if mode in ("dominant", "states"):
                    mol_out, best_energy = _best_conf(mol_h, cids, energies)
                    out = Chem.RWMol(mol_out)
                    out.SetProp("_Name", mol_id)
                    out.SetIntProp("tautomer_id", t_idx)
                    out.SetDoubleProp("tautomer_population", round(t_pop, 4))
                    out.SetIntProp("stereo_id", s_idx)
                    out.SetDoubleProp("state_population", round(state_pop, 6))
                    out.SetDoubleProp("state_penalty", round(state_penalty, 3))
                    out.SetDoubleProp("mmff_energy", round(best_energy, 4))
                    out.SetProp("embedder", embedder)
                    sdf_texts.append(_mol_to_sdf(out.GetMol()))
                else:  # ensemble
                    _order   = sorted(range(len(cids)), key=lambda i: energies[i])
                    cids     = [cids[i]     for i in _order]
                    energies = [energies[i] for i in _order]
                    if max_confs_out is not None:
                        cids     = cids[:max_confs_out]
                        energies = energies[:max_confs_out]
                    for c_idx, (conf_id, energy) in enumerate(zip(cids, energies), start=1):
                        out = Chem.RWMol(mol_h)
                        out.SetProp("_Name", mol_id)
                        out.SetIntProp("tautomer_id", t_idx)
                        out.SetDoubleProp("tautomer_population", round(t_pop, 4))
                        out.SetIntProp("stereo_id", s_idx)
                        out.SetDoubleProp("state_population", round(state_pop, 6))
                        out.SetDoubleProp("state_penalty", round(state_penalty, 3))
                        out.SetDoubleProp("mmff_energy", round(energy, 4))
                        out.SetIntProp("conf_id", c_idx)
                        out.SetProp("embedder", embedder)
                        sdf_texts.append(_mol_to_sdf(out.GetMol(), conf_id=conf_id))

    except Exception as e:
        failures.append({"id": mol_id, "stage": "conforge", "reason": str(e)})

    return sdf_texts, failures


# ─────────────────────────────────────────────────────────────────────────────
# cxcalc tautomer enumeration
# ─────────────────────────────────────────────────────────────────────────────

def _cxcalc_batch(
    mols_with_ids: list, pH: float, max_tautomers: int, min_population: float,
) -> dict:
    """Run cxcalc on a list of (mol, mol_id) pairs. Re-indexes internally."""
    idx_to_id = {str(i): mol_id for i, (_, mol_id) in enumerate(mols_with_ids)}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False, prefix="lgp3_") as fh:
        tmp_path = fh.name
        for i, (mol, _) in enumerate(mols_with_ids):
            fh.write(f"{Chem.MolToSmiles(mol)}\t{i}\n")
    try:
        result = subprocess.run(
            ["cxcalc", "tautomers", "-D", "true", "-n", "true",
             "-H", str(pH), "-m", str(max_tautomers), "-f", "sdf", tmp_path],
            capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(f"cxcalc failed (rc={result.returncode}):\n{result.stderr.strip()}")
        sdf_text = result.stdout
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    return _parse_tautomer_sdf(sdf_text, idx_to_id, min_population)


def _run_cxcalc_chunk(
    chunk: list, pH: float, max_tautomers: int, min_population: float,
    retry_batch: int = 50,
) -> dict:
    """Run cxcalc on a chunk with fallback retry for empty results.

    If one molecule crashes cxcalc mid-batch, all subsequent molecules in the
    batch produce no output.  Detected by checking for empty groups after the
    first pass; empties are retried in smaller sub-batches, then individually.
    """
    groups = _cxcalc_batch(chunk, pH, max_tautomers, min_population)

    empties = [(mol, mol_id) for mol, mol_id in chunk if not groups[mol_id]]
    if not empties:
        return groups

    n_recovered = 0
    # Second pass: small sub-batches
    for i in range(0, len(empties), retry_batch):
        sub = empties[i:i + retry_batch]
        sub_groups = _cxcalc_batch(sub, pH, max_tautomers, min_population)
        for _, mol_id in sub:
            if sub_groups[mol_id]:
                groups[mol_id] = sub_groups[mol_id]
                n_recovered += 1

    # Third pass: any still-empty molecules run individually
    still_empty = [(mol, mol_id) for mol, mol_id in empties if not groups[mol_id]]
    for mol, mol_id in still_empty:
        single = _cxcalc_batch([(mol, mol_id)], pH, max_tautomers, min_population)
        if single[mol_id]:
            groups[mol_id] = single[mol_id]
            n_recovered += 1

    if n_recovered:
        print(f"[INFO]  cxcalc fallback recovered {n_recovered}/{len(empties)} empty result(s)")
    return groups


def _run_cxcalc_tautomers(
    mols_with_ids: list, pH: float = 7.4,
    max_tautomers: int = 200, min_population: float = 0.01,
    chunk_size: int = 5000,
) -> dict:
    if not _CXCALC:
        raise RuntimeError("cxcalc not found on PATH")

    chunks = [mols_with_ids[i:i + chunk_size]
              for i in range(0, len(mols_with_ids), chunk_size)]
    n_chunks = len(chunks)

    groups: dict = {}
    for c_idx, chunk in enumerate(chunks, start=1):
        if n_chunks > 1:
            print(f"[INFO]  cxcalc chunk {c_idx}/{n_chunks} ({len(chunk)} molecules)")
        groups.update(_run_cxcalc_chunk(chunk, pH, max_tautomers, min_population))
    return groups


def _parse_tautomer_sdf(sdf_text: str, idx_to_id: dict, min_population: float) -> dict:
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(sdf_text.encode()), removeHs=False)
    groups   = {mol_id: [] for mol_id in idx_to_id.values()}
    skipped  = 0
    for mol in supplier:
        if mol is None:
            skipped += 1; continue
        idx = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        if idx not in idx_to_id:
            skipped += 1; continue
        mol_id = idx_to_id[idx]
        try:
            population = int(mol.GetPropsAsDict().get("TAUTOMER_DISTRIBUTION", 0)) / 100.0
        except (ValueError, TypeError):
            skipped += 1; continue
        if population < min_population:
            continue
        mol.ClearProp("TAUTOMER_DISTRIBUTION")
        if mol.HasProp("_MOLCOUNT"):
            mol.ClearProp("_MOLCOUNT")
        groups[mol_id].append((mol, population))
    if skipped:
        print(f"[WARNING] cxcalc SDF parser: skipped {skipped} unreadable record(s)")
    for mid in groups:
        groups[mid].sort(key=lambda x: x[1], reverse=True)
    return groups


# ─────────────────────────────────────────────────────────────────────────────
# Stereo enumeration + input parsing
# ─────────────────────────────────────────────────────────────────────────────

def _enumerate_stereo(mol: "Chem.Mol", max_isomers: int = 32) -> list:
    opts    = StereoEnumerationOptions(unique=True, onlyUnassigned=True, maxIsomers=max_isomers)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    return isomers if isomers else [mol]


def _read_input(path: Path) -> "tuple[list, list[dict]]":
    """Returns (mols, failures) where failures are {"id", "stage", "reason"} dicts."""
    mols     = []
    failures = []
    if path.suffix.lower() in (".sdf", ".sd"):
        for i, mol in enumerate(Chem.SDMolSupplier(str(path), removeHs=True)):
            if mol is None:
                failures.append({"id": f"record_{i+1}", "stage": "input", "reason": "unparseable SDF record"})
                continue
            name = (mol.GetPropsAsDict().get("_Name") or f"mol_{i}").strip() or f"mol_{i}"
            mols.append((_strip_salts(mol), name))
    else:
        with open(path) as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.replace(",", "\t").split()
                if len(parts) < 1:
                    continue
                # Auto-detect column order: try col0 as SMILES, then col1
                mol = Chem.MolFromSmiles(parts[0])
                if mol is not None:
                    smi, name = parts[0], (parts[1] if len(parts) > 1 else f"mol_{i}")
                elif len(parts) > 1:
                    mol = Chem.MolFromSmiles(parts[1])
                    if mol is not None:
                        smi, name = parts[1], parts[0]
                    else:
                        if i == 0:
                            continue  # silently skip header
                        name = parts[0]
                        failures.append({"id": name, "stage": "input", "reason": f"invalid SMILES: {parts[1]!r}"})
                        continue
                else:
                    if i == 0:
                        continue
                    failures.append({"id": f"line_{i+1}", "stage": "input", "reason": f"invalid SMILES: {parts[0]!r}"})
                    continue
                mol.SetProp("_Name", name)
                mols.append((_strip_salts(mol), name))
    return mols, failures


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=_HelpFmt,
    )
    input_arg = p.add_argument(
        "-i", "--input", type=Path, required=True, metavar="FILE",
        help="Input SMILES (.smi/.csv/.txt) or SDF (.sdf/.sd) file",
    )
    p.add_argument("-o", "--output", type=Path, default=None, metavar="FILE",
                   help="Output SDF (default: <input_stem>_prep.sdf)")
    p.add_argument("--log", type=Path, default=None, metavar="FILE",
                   help="Log file (default: <output_stem>.log)")
    p.add_argument("--pH", type=float, default=7.4, metavar="FLOAT",
                   help="Target pH for tautomer/protomer distribution (default: 7.4)")
    p.add_argument("--min-population", type=float, default=0.10, metavar="FLOAT",
                   help="Minimum fractional population threshold (default: 0.10 = 10%%)")
    p.add_argument("--max-tautomers", type=int, default=200, metavar="INT",
                   help="Max tautomers per molecule passed to cxcalc (default: 200)")
    p.add_argument("--cxcalc-chunk", type=int, default=5000, metavar="INT",
                   help="Molecules per cxcalc batch (default: 5000)")
    p.add_argument("--mode", choices=["dominant", "states", "ensemble"],
                   default="dominant",
                   help="dominant: 1 best conformer; states: 1 conf per state; "
                        "ensemble: multiple conformers per state (default: dominant)")
    p.add_argument("--num-confs", type=int, default=None, metavar="INT",
                   help="Max conformers CONFORGE may output per isomer "
                        "(default: 50 for dominant/states, 100 for ensemble)")
    p.add_argument("--conforge-ewin", type=float, default=5.0, metavar="FLOAT",
                   dest="conforge_ewin",
                   help="CONFORGE energy window kcal/mol for dominant/states "
                        "(default: 5.0; ensemble uses --ewin)")
    p.add_argument("--rmsd-threshold", type=float, default=0.5, metavar="FLOAT",
                   help="RMSD threshold Å for CONFORGE deduplication (default: 0.5)")
    p.add_argument("--ewin", type=float, default=5.0, metavar="FLOAT",
                   help="Energy window kcal/mol for CONFORGE in ensemble mode "
                        "(default: 5.0)")
    p.add_argument("--max-confs-out", type=int, default=50, metavar="INT",
                   help="Max conformers to write per tautomer/stereo state in ensemble mode "
                        "(default: 50); lowest-energy conformers kept first")
    p.add_argument("--no-stereo", action="store_true",
                   help="Skip stereo enumeration (keep input stereo as-is)")
    p.add_argument("--max-stereo", type=int, default=32, metavar="INT",
                   help="Max stereoisomers per molecule for unspecified centers (default: 32)")
    p.add_argument("--conforge-timeout", type=int, default=60, metavar="SEC",
                   dest="conforge_timeout",
                   help="CONFORGE per-molecule timeout in seconds (default: 60); "
                        "on timeout ETKDGv3 fallback is tried automatically")
    p.add_argument("--max-sampled", type=int, default=500, metavar="INT",
                   dest="max_sampled",
                   help="CONFORGE max internal conformers sampled per molecule (default: 500; "
                        "increase for better global minimum coverage of very flexible molecules)")
    p.add_argument("-j", "--jobs", type=int,
                   default=max(1, (os.cpu_count() or 1) - 1),
                   metavar="N",
                   help="Parallel worker processes for CONFORGE (default: cpu_count-1)")

    if _ARGCOMPLETE:
        input_arg.completer = FilesCompleter(allowednames=(".smi", ".csv", ".txt", ".sdf", ".sd"))
    return p


def main() -> None:
    p = build_parser()
    if _ARGCOMPLETE:
        argcomplete.autocomplete(p)
    args = p.parse_args()
    t0   = time.time()

    if not _RDKIT:
        print("[ERROR] RDKit not available", file=sys.stderr); sys.exit(1)
    if not _CXCALC:
        print("[ERROR] cxcalc not found on PATH", file=sys.stderr); sys.exit(1)
    if not _CDPKIT:
        print("[WARNING] CDPKit not available — CONFORGE disabled, using ETKDGv3 for all molecules")

    if args.output is None:
        args.output = args.input.with_name(args.input.stem + "_prep.sdf")

    if args.log is None:
        args.log = args.output.with_suffix(".log")
    _log_fh = open(args.log, "w")
    sys.stdout = _Tee(sys.__stdout__, _log_fh)

    _conforge_ewin = args.ewin if args.mode == "ensemble" else args.conforge_ewin
    _conforge_nmax = args.num_confs if args.num_confs is not None else (
        100 if args.mode == "ensemble" else 50
    )

    print(f"[CONFIG] Input:            {args.input}")
    print(f"[CONFIG] Output:           {args.output}")
    print(f"[CONFIG] pH:               {args.pH}")
    print(f"[CONFIG] Max tautomers:    {args.max_tautomers}")
    print(f"[CONFIG] Mode:             {args.mode}")
    print(f"[CONFIG] Stereo enum:      {'off' if args.no_stereo else f'on (unspecified centers, max {args.max_stereo})'}")
    print(f"[CONFIG] CONFORGE max out: {_conforge_nmax} conformer(s)")
    print(f"[CONFIG] CONFORGE ewin:    {_conforge_ewin} kcal/mol")
    print(f"[CONFIG] RMSD threshold:   {args.rmsd_threshold} Å")
    if args.mode == "ensemble":
        print(f"[CONFIG] Max confs out:    {args.max_confs_out}")
    if args.mode != "dominant":
        print(f"[CONFIG] Min population:   {args.min_population:.0%}")
    print(f"[CONFIG] CONFORGE timeout: {args.conforge_timeout}s")
    print(f"[CONFIG] Max sampled:      {args.max_sampled}")
    print(f"[CONFIG] Jobs:             {args.jobs}")
    print(f"[CONFIG] cxcalc:           {shutil.which('cxcalc')}")

    all_failures: list[dict] = []

    print(f"\n[STEP 1/3] Reading input: {args.input}")
    raw_mols, input_failures = _read_input(args.input)
    all_failures.extend(input_failures)
    if not raw_mols:
        print("[ERROR] No valid molecules in input", file=sys.stderr); sys.exit(1)
    print(f"[INFO]  {len(raw_mols)} molecule(s) read"
          + (f" ({len(input_failures)} failed to parse)" if input_failures else ""))

    seen: dict[str, int] = {}
    mols_with_ids = []
    for mol, name in raw_mols:
        if name in seen:
            seen[name] += 1; uid = f"{name}_{seen[name]}"
        else:
            seen[name] = 0; uid = name
        mols_with_ids.append((mol, uid))

    print(f"\n[STEP 2/3] cxcalc tautomers -D true -n true -H {args.pH}")
    _pop_filter    = 0.0 if args.mode == "dominant" else args.min_population
    tautomer_groups = _run_cxcalc_tautomers(
        mols_with_ids, pH=args.pH,
        max_tautomers=args.max_tautomers, min_population=_pop_filter,
        chunk_size=args.cxcalc_chunk,
    )
    if args.mode == "dominant":
        tautomer_groups = {mid: v[:1] for mid, v in tautomer_groups.items()}
        print(f"[INFO]  dominant state selected for "
              f"{sum(len(v) for v in tautomer_groups.values())} molecule(s)")
    else:
        n_states = sum(len(v) for v in tautomer_groups.values())
        n_empty  = sum(1 for v in tautomer_groups.values() if not v)
        print(f"[INFO]  {n_states} tautomer/protomer state(s) above "
              f"{args.min_population:.0%} threshold")
        if n_empty:
            print(f"[WARNING] {n_empty} molecule(s) produced no states above threshold")

    print(f"\n[STEP 3/3] Stereo enumeration + CONFORGE 3D ({args.jobs} worker(s))")

    jobs = []
    for mol, mol_id in mols_with_ids:
        states = tautomer_groups.get(mol_id, [])
        if not states:
            all_failures.append({"id": mol_id, "stage": "tautomers", "reason": "no states above population threshold"})
            continue
        jobs.append({
            "mol_id":        mol_id,
            "states":        [(t_idx, t_mol.ToBinary(), t_pop)
                              for t_idx, (t_mol, t_pop) in enumerate(states, start=1)],
            "no_stereo":     args.no_stereo,
            "max_stereo":    args.max_stereo,
            "nmax":          _conforge_nmax,
            "ewin":          _conforge_ewin,
            "rmsd_thr":      args.rmsd_threshold,
            "mode":          args.mode,
            "max_confs_out": args.max_confs_out,
            "timeout_ms":    args.conforge_timeout * 1000,
            "max_sampled":   args.max_sampled,
        })

    n_written = 0

    with open(str(args.output), "w") as out_fh:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(_mol_job, job): job for job in jobs}
            for fut in _tqdm(as_completed(futures), total=len(jobs),
                             desc="Molecules", disable=not _TQDM):
                sdf_texts, job_failures = fut.result()
                all_failures.extend(job_failures)
                for text in sdf_texts:
                    out_fh.write(text)
                    n_written += 1

    if all_failures:
        failed_path = args.output.with_name(args.output.stem + "_failed.tsv")
        with open(failed_path, "w") as f:
            f.write("id\tstage\treason\n")
            for rec in all_failures:
                f.write(f"{rec['id']}\t{rec['stage']}\t{rec['reason']}\n")
        print(f"[WARNING] {len(all_failures)} failure(s) — see {failed_path}")

    _written_label = "Conformers written" if args.mode == "ensemble" else "States written"
    elapsed = time.time() - t0
    print(f"\n{'─' * 60}")
    print(f"[SUMMARY] Input molecules:   {len(raw_mols)}")
    print(f"[SUMMARY] {_written_label}:  {n_written}")
    print(f"[SUMMARY] Failed:            {len(all_failures)}")
    print(f"[SUMMARY] Output:            {args.output}")
    print(f"[SUMMARY] Time:              {format_time(elapsed)}")
    print(f"[SUMMARY] Log:               {args.log}")
    sys.stdout = sys.__stdout__
    _log_fh.close()


if __name__ == "__main__":
    main()
