#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
python3 -m venv --system-site-packages .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install --no-deps --no-build-isolation -e .
.venv/bin/python -m pytest -q
