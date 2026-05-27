# ligprepper

Web application for filtering and preparing ligand libraries for *in silico* applications.

`ligprepper` is a FastAPI front-end that wraps two command-line tools shipped in
this repository:

| Script | Role |
|---|---|
| [`ligfilter.py`](ligfilter_README.md) | Multi-stage ligand filter — PAINS / REOS / property ranges / Lipinski / outlier removal |
| [`ligprep.py`](ligprep_README.md) | Open-source LigPrep/Epik equivalent — tautomer + protomer enumeration (cxcalc), 3D conformers (CDPKit CONFORGE), optional xTB refinement |

The UI exposes both scripts as **two independent tabs**. Each tab uploads an
input file, posts a form of CLI options, runs the corresponding script as a
subprocess, and streams stdout to the browser via Server-Sent Events. The
produced output file is exposed for download.

---

## Architecture

```
┌──────────────────────────────────┐
│  templates/index.html            │  two-tab UI, SSE log viewer
└──────────────┬───────────────────┘
               │  POST /run/ligfilter | /run/ligprep
               ▼
┌──────────────────────────────────┐
│  ligprepper_webapp.py            │  FastAPI app
│   ├─ /run/{ligfilter,ligprep}    │  spawn subprocess, return job_id
│   ├─ /jobs/{id}/status           │  poll status
│   ├─ /jobs/{id}/stream           │  SSE stdout stream
│   ├─ /jobs/{id}/download         │  output file
│   └─ /jobs/{id}/log              │  combined log
└──────────────┬───────────────────┘
               │  subprocess
               ▼
        ligfilter.py / ligprep.py
```

Each job gets its own directory under `jobs/<uuid>/` containing the uploaded
input, the script's output, and a combined stdout/stderr log.

---

## Installation

`ligprepper` runs in a dedicated conda environment so that RDKit, CDPKit, and
xTB stay isolated from the system Python.

```bash
conda create -n ligprepper python=3.11
conda activate ligprepper
conda install -c conda-forge rdkit xtb
pip install -r requirements.txt        # CDPKit + FastAPI + helpers
```

`cxcalc` (ChemAxon JChem Suite) must be installed separately and available on
`$PATH` — it is required by `ligprep.py` for tautomer/protomer enumeration and
is not redistributable via conda or pip.

---

## Running

```bash
./run.sh
```

The script activates the `ligprepper` conda env (via `conda run`, so it works
whether or not the env is active in your shell) and launches uvicorn on
`0.0.0.0:5009`.

Override defaults via environment variables:

| Variable | Default | Description |
|---|---|---|
| `LIGPREPPER_ENV` | `ligprepper` | Conda env to activate |
| `LIGPREPPER_HOST` | `0.0.0.0` | uvicorn bind address |
| `LIGPREPPER_PORT` | `5009` | uvicorn port |
| `LIGPREPPER_SCRIPT_PYTHON` | `sys.executable` | Python interpreter used for `ligfilter.py` / `ligprep.py` subprocesses — override if the scripts need a different env than the webapp |

Open <http://localhost:5009> in a browser.

---

## Repository layout

```
ligprepper/
├── ligprepper_webapp.py        FastAPI app
├── run.sh                      conda-aware launcher
├── requirements.txt            pip deps + conda recipe in comments
├── ligfilter.py                CLI filter script
├── ligfilter_README.md         docs for ligfilter.py
├── ligprep.py                  CLI ligand-prep script
├── ligprep_README.md           docs for ligprep.py
├── PAINS.txt                   PAINS SMARTS patterns (auto-discovered)
├── REOS.txt                    REOS SMARTS patterns (auto-discovered)
├── custom_filters.txt          user-defined filters (auto-discovered)
├── templates/
│   └── index.html              two-tab UI
├── static/                     static assets
└── jobs/                       per-job working directories (gitignored)
```

The three filter files (`PAINS.txt`, `REOS.txt`, `custom_filters.txt`) are
auto-discovered by `ligfilter.py` when present in the script directory.

---

## Scripts standalone

`ligfilter.py` and `ligprep.py` are self-contained command-line tools and can
be used without the webapp. See their respective READMEs for full CLI
documentation:

- [`ligfilter_README.md`](ligfilter_README.md)
- [`ligprep_README.md`](ligprep_README.md)
