#!/usr/bin/env python3
"""
ligprepper web app — FastAPI front-end for ligprep.py and ligfilter.py.

Two independent tabs (ligfilter, ligprep). Each tab uploads an input file,
posts a form of CLI options, and the backend runs the corresponding Python
script as a subprocess. Stdout/stderr is tee'd to a per-job log file and
streamed to the browser via Server-Sent Events. The produced output file is
exposed for download.
"""

from __future__ import annotations

import asyncio
import logging
import os
import re
import shutil
import subprocess
import sys
import uuid
from contextlib import asynccontextmanager
from datetime import datetime
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates


# ─── Paths & config ───────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent
JOBS_DIR = ROOT / "jobs"
JOBS_DIR.mkdir(exist_ok=True)

LIGFILTER = ROOT / "ligfilter.py"
LIGPREP = ROOT / "ligprep.py"

# Python interpreter used to run the scripts. By default use the same one that
# runs the webapp (we assume the ligprepper conda env has rdkit/cdpkit/cxcalc
# all available). Override with LIGPREPPER_SCRIPT_PYTHON if a different env
# must run the scripts.
SCRIPT_PYTHON = os.environ.get("LIGPREPPER_SCRIPT_PYTHON", sys.executable)

LOG_FILE = ROOT / "ligprepper_webapp.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(LOG_FILE, encoding="utf-8"),
    ],
)
logger = logging.getLogger("ligprepper")


# ─── Job tracking ─────────────────────────────────────────────────────────────

class Job:
    __slots__ = ("id", "dir", "script", "input_path", "output_path",
                 "log_path", "proc", "returncode", "started", "finished")

    def __init__(self, script: str):
        self.id = uuid.uuid4().hex[:12]
        self.script = script
        self.dir = JOBS_DIR / self.id
        self.dir.mkdir(parents=True, exist_ok=True)
        self.input_path: Optional[Path] = None
        self.output_path: Optional[Path] = None
        self.log_path: Path = self.dir / "run.log"
        self.proc: Optional[asyncio.subprocess.Process] = None
        self.returncode: Optional[int] = None
        self.started: datetime = datetime.utcnow()
        self.finished: Optional[datetime] = None


JOBS: dict[str, Job] = {}


# ─── Helpers ──────────────────────────────────────────────────────────────────

_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


def _safe_basename(name: str, default: str = "input") -> str:
    base = Path(name).name
    base = _SAFE_NAME.sub("_", base).strip("._-")
    return base or default


def _truthy(v: Optional[str]) -> bool:
    return v is not None and str(v).lower() in ("1", "true", "on", "yes")


def _validate_range(s: str) -> str:
    """Allow MIN:MAX / MIN: / :MAX / exact value. Strict regex to avoid shell tricks."""
    if not re.fullmatch(r"-?\d+(?:\.\d+)?(?::-?\d+(?:\.\d+)?)?|:-?\d+(?:\.\d+)?|-?\d+(?:\.\d+)?:", s):
        raise HTTPException(status_code=400, detail=f"invalid range value: {s!r}")
    return s


async def _run_subprocess(job: Job, cmd: list[str]) -> None:
    """Spawn the script and stream stdout+stderr into job.log_path."""
    logger.info("job %s start: %s", job.id, " ".join(cmd))
    log_fh = open(job.log_path, "wb")
    try:
        job.proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            cwd=str(job.dir),
        )
        job.returncode = await job.proc.wait()
    finally:
        log_fh.close()
        job.finished = datetime.utcnow()
        logger.info("job %s done rc=%s", job.id, job.returncode)


# ─── FastAPI app ──────────────────────────────────────────────────────────────

@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("ligprepper webapp starting — script python: %s", SCRIPT_PYTHON)
    if not LIGFILTER.exists() or not LIGPREP.exists():
        logger.warning("ligprep.py or ligfilter.py missing in %s", ROOT)
    yield
    logger.info("ligprepper webapp shutting down")


app = FastAPI(title="ligprepper", lifespan=lifespan)

(ROOT / "static").mkdir(exist_ok=True)
app.mount("/static", StaticFiles(directory=str(ROOT / "static")), name="static")
templates = Jinja2Templates(directory=str(ROOT / "templates"))


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(request, "index.html")


@app.get("/health")
async def health():
    return {
        "status": "ok",
        "script_python": SCRIPT_PYTHON,
        "ligfilter": LIGFILTER.exists(),
        "ligprep": LIGPREP.exists(),
        "jobs": len(JOBS),
    }


# ─── ligfilter endpoint ───────────────────────────────────────────────────────

@app.post("/run/ligfilter")
async def run_ligfilter(
    input_file: UploadFile = File(...),
    # output extension override (smi / sdf / "" = auto)
    output_format: Optional[str] = Form(None),
    # preprocessing
    no_strip: Optional[str] = Form(None),
    no_neutralize: Optional[str] = Form(None),
    no_tautomer_canon: Optional[str] = Form(None),
    no_unique: Optional[str] = Form(None),
    # structural alerts
    pains: Optional[str] = Form(None),
    reos: Optional[str] = Form(None),
    custom: Optional[str] = Form(None),
    # drug-likeness
    ro5: Optional[str] = Form(None),
    ro5_strict: Optional[str] = Form(None),
    ro3: Optional[str] = Form(None),
    # property ranges
    mw: Optional[str] = Form(None),
    logp: Optional[str] = Form(None),
    hba: Optional[str] = Form(None),
    hbd: Optional[str] = Form(None),
    rb: Optional[str] = Form(None),
    tpsa: Optional[str] = Form(None),
    qed: Optional[str] = Form(None),
    chiral: Optional[str] = Form(None),
    ha: Optional[str] = Form(None),
    # outlier
    outlier: Optional[str] = Form(None),
    outlier_iqr: Optional[str] = Form(None),
    outlier_props: Optional[str] = Form(None),
    jobs: Optional[str] = Form(None),
):
    job = Job(script="ligfilter")
    src_name = _safe_basename(input_file.filename or "input.smi")
    in_path = job.dir / src_name
    with open(in_path, "wb") as fh:
        shutil.copyfileobj(input_file.file, fh)
    job.input_path = in_path

    # Determine output extension
    if output_format in ("smi", "sdf"):
        out_ext = "." + output_format
    else:
        out_ext = in_path.suffix.lower() or ".smi"
    out_path = job.dir / (in_path.stem + "_filtered" + out_ext)
    job.output_path = out_path

    cmd: list[str] = [SCRIPT_PYTHON, str(LIGFILTER),
                      "-i", str(in_path), "-o", str(out_path)]

    if _truthy(no_strip):          cmd.append("--no-strip")
    if _truthy(no_neutralize):     cmd.append("--no-neutralize")
    if _truthy(no_tautomer_canon): cmd.append("--no-tautomer-canon")
    if _truthy(no_unique):         cmd.append("--no-unique")

    if _truthy(pains):  cmd.append("--pains")
    if _truthy(reos):   cmd.append("--reos")
    if _truthy(custom): cmd.append("--custom")

    if _truthy(ro5_strict):
        cmd.append("--ro5-strict")
    elif _truthy(ro5):
        cmd.append("--ro5")
    if _truthy(ro3):
        cmd.append("--ro3")

    range_args = {
        "--mw": mw, "--logp": logp, "--hba": hba, "--hbd": hbd, "--rb": rb,
        "--tpsa": tpsa, "--qed": qed, "--chiral": chiral, "--ha": ha,
    }
    for flag, val in range_args.items():
        if val and val.strip():
            cmd += [flag, _validate_range(val.strip())]

    if _truthy(outlier):
        cmd.append("--outlier")
        if outlier_iqr and outlier_iqr.strip():
            float(outlier_iqr)  # validates; raises ValueError on bad input
            cmd += ["--outlier-iqr", outlier_iqr.strip()]
        if outlier_props and outlier_props.strip():
            props = outlier_props.strip()
            if not re.fullmatch(r"[a-z,]+", props):
                raise HTTPException(400, detail="outlier-props must be comma-separated lowercase names")
            cmd += ["--outlier-props", props]

    if jobs and jobs.strip():
        int(jobs)
        cmd += ["--jobs", jobs.strip()]

    JOBS[job.id] = job
    asyncio.create_task(_run_subprocess(job, cmd))
    return {"job_id": job.id, "output_name": out_path.name}


# ─── ligprep endpoint ─────────────────────────────────────────────────────────

@app.post("/run/ligprep")
async def run_ligprep(
    input_file: UploadFile = File(...),
    pH: Optional[str] = Form(None),
    mode: Optional[str] = Form(None),
    min_population: Optional[str] = Form(None),
    max_tautomers: Optional[str] = Form(None),
    num_confs: Optional[str] = Form(None),
    conforge_ewin: Optional[str] = Form(None),
    rmsd_threshold: Optional[str] = Form(None),
    ewin: Optional[str] = Form(None),
    max_confs_out: Optional[str] = Form(None),
    no_stereo: Optional[str] = Form(None),
    max_stereo: Optional[str] = Form(None),
    jobs: Optional[str] = Form(None),
):
    job = Job(script="ligprep")
    src_name = _safe_basename(input_file.filename or "input.smi")
    in_path = job.dir / src_name
    with open(in_path, "wb") as fh:
        shutil.copyfileobj(input_file.file, fh)
    job.input_path = in_path

    out_path = job.dir / (in_path.stem + "_prep.sdf")
    job.output_path = out_path

    cmd: list[str] = [SCRIPT_PYTHON, str(LIGPREP),
                      "-i", str(in_path), "-o", str(out_path)]

    def _float_arg(flag: str, val: Optional[str]):
        if val and val.strip():
            float(val)
            cmd.extend([flag, val.strip()])

    def _int_arg(flag: str, val: Optional[str]):
        if val and val.strip():
            int(val)
            cmd.extend([flag, val.strip()])

    _float_arg("--pH", pH)
    if mode in ("dominant", "states", "ensemble"):
        cmd += ["--mode", mode]
    _float_arg("--min-population", min_population)
    _int_arg("--max-tautomers", max_tautomers)
    _int_arg("--num-confs", num_confs)
    _float_arg("--conforge-ewin", conforge_ewin)
    _float_arg("--rmsd-threshold", rmsd_threshold)
    _float_arg("--ewin", ewin)
    _int_arg("--max-confs-out", max_confs_out)
    if _truthy(no_stereo):
        cmd.append("--no-stereo")
    _int_arg("--max-stereo", max_stereo)
    _int_arg("--jobs", jobs)

    JOBS[job.id] = job
    asyncio.create_task(_run_subprocess(job, cmd))
    return {"job_id": job.id, "output_name": out_path.name}


# ─── Status / log streaming / download ────────────────────────────────────────

@app.get("/jobs/{job_id}/status")
async def job_status(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "unknown job")
    running = job.returncode is None
    return {
        "job_id": job.id,
        "script": job.script,
        "running": running,
        "returncode": job.returncode,
        "output_exists": job.output_path is not None and job.output_path.exists(),
        "output_name": job.output_path.name if job.output_path else None,
    }


@app.get("/jobs/{job_id}/stream")
async def job_stream(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "unknown job")

    async def gen():
        # Wait briefly for the log file to appear.
        for _ in range(50):
            if job.log_path.exists():
                break
            await asyncio.sleep(0.05)
        # Stream the log; tail-follow until process finishes and file is fully drained.
        with open(job.log_path, "rb") as fh:
            while True:
                chunk = fh.read(4096)
                if chunk:
                    text = chunk.decode("utf-8", errors="replace")
                    # SSE: prefix every line with "data: " and end with blank line
                    for line in text.splitlines():
                        yield f"data: {line}\n\n".encode("utf-8")
                else:
                    if job.returncode is not None:
                        # final drain done
                        break
                    await asyncio.sleep(0.25)
        rc = job.returncode if job.returncode is not None else -1
        yield f"event: done\ndata: {rc}\n\n".encode("utf-8")

    return StreamingResponse(gen(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache",
                                      "X-Accel-Buffering": "no"})


@app.get("/jobs/{job_id}/download")
async def job_download(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "unknown job")
    if job.output_path is None or not job.output_path.exists():
        raise HTTPException(404, "output not yet available")
    return FileResponse(str(job.output_path), filename=job.output_path.name,
                        media_type="application/octet-stream")


@app.get("/jobs/{job_id}/log")
async def job_log(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "unknown job")
    if not job.log_path.exists():
        raise HTTPException(404, "log not available")
    if job.output_path is not None:
        download_name = job.output_path.stem + ".log"
    else:
        download_name = f"{job.id}.log"
    return FileResponse(str(job.log_path), filename=download_name,
                        media_type="text/plain")


if __name__ == "__main__":
    import uvicorn
    host = os.environ.get("LIGPREPPER_HOST", "0.0.0.0")
    port = int(os.environ.get("LIGPREPPER_PORT", "5009"))
    uvicorn.run("ligprepper_webapp:app", host=host, port=port, reload=False)
