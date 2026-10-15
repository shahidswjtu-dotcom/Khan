#!/usr/bin/env bash
set -euo pipefail
cd /workspace/bowen_reproduction
python3 -m venv --system-site-packages .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c 'from bowen_pd.model import create,evaluate; import numpy as np; m=create(); r=evaluate(m,np.zeros_like(m.xyz),np.zeros_like(m.steel_xyz),m.hammer_z_initial,m.p.dt); assert np.max(np.abs(r["fc"]))<1e-7'
