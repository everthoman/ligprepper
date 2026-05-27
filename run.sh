#!/usr/bin/env bash
# Launch the ligprepper webapp inside the `ligprepper` conda env.
set -euo pipefail

cd "$(dirname "$0")"

# Allow overriding env, host, port from environment
ENV_NAME="${LIGPREPPER_ENV:-ligprepper}"
HOST="${LIGPREPPER_HOST:-0.0.0.0}"
PORT="${LIGPREPPER_PORT:-5009}"

# Use `conda run` so this script works whether or not the user has run
# `conda activate` in the current shell.
CONDA_BASE="$(conda info --base 2>/dev/null || echo "$HOME/Programs/miniconda3")"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

exec python -m uvicorn ligprepper_webapp:app \
    --host "$HOST" --port "$PORT" --log-level info
