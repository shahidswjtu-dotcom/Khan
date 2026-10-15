"""Extract full-run extrema only from a successfully completed solver run."""
import argparse
import json
import os
from pathlib import Path
import zipfile

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("run", type=Path)
a = p.parse_args()
m = json.loads((a.run / "manifest.json").read_text())
if m.get("exit_code") != 0:
    raise SystemExit("Run has not completed successfully; no full-run summary created")
h = np.genfromtxt(a.run / "history.csv", delimiter=",", skip_header=1, ndmin=2)
s = np.loadtxt(a.run / "final_state.txt")
assert h.shape == (m["steps"], 6)
assert np.isfinite(h).all() and np.isfinite(s).all()
assert int(h[-1, 0]) == m["steps"]
f = int(np.argmax(abs(h[:, 4])))
d = int(np.argmax(abs(h[:, 2])))
report = dict(completed_steps=m["steps"], simulated_time_s=float(h[-1, 1]),
              peak_absolute_vertical_contact_kN=float(abs(h[f, 4])/1000),
              peak_force_time_s=float(h[f, 1]),
              maximum_absolute_original_probe_mm=float(abs(h[d, 2])*1000),
              peak_probe_time_s=float(h[d, 1]),
              final_prescribed_indenter_displacement_mm=float(h[-1, 5]*1000),
              final_maximum_nodal_displacement_mm=float(np.max(np.linalg.norm(s[:, 4:7], axis=1))*1000),
              final_mean_concrete_damage=float(np.mean(s[362:, -1])),
              all_finite=True, bowen_chapter3_model=False,
              corrections=m.get("corrections", []),
              final_maximum_steel_displacement_m=float(np.max(np.linalg.norm(s[:362, 4:7], axis=1))),
              final_maximum_concrete_displacement_m=float(np.max(np.linalg.norm(s[362:, 4:7], axis=1))),
              final_maximum_steel_speed_m_s=float(np.max(np.linalg.norm(s[:362, 7:10], axis=1))),
              final_maximum_concrete_speed_m_s=float(np.max(np.linalg.norm(s[362:, 7:10], axis=1))),
              physical_validity="Requires force, energy, deformation and convergence audit; finite values alone do not establish validity",
              measurement="Original single-node pre-update displacement; post-projection cylindrical reaction")
if report["final_maximum_steel_displacement_m"] > 100:
    report["physical_validity"] = "Rejected: steel displacement exceeds 100 m for a roughly 3.46 m beam"
(a.run / "results.json").write_text(json.dumps(report, indent=2)+"\n")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/bowen-fortran-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/bowen-fortran-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
axes[0].plot(h[:, 1], -h[:, 4]/1000)
axes[0].set_ylabel("Signed downward reaction (kN)")
axes[1].plot(h[:, 1], -h[:, 2]*1000, label="Original probe")
axes[1].plot(h[:, 1], -h[:, 5]*1000, label="Prescribed next-step indenter")
axes[1].set_ylabel("Downward displacement (mm)")
axes[1].set_xlabel("Time (s)")
axes[1].legend()
fig.suptitle("Supplied Fortran model — full loading history")
fig.tight_layout()
fig.savefig(a.run / "results.png", dpi=180)
with zipfile.ZipFile(a.run / "complete_results.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for name in ["manifest.json", "history.csv", "final_state.txt", "results.json", "results.png", "runtime.log", "solver.instrumented.f90"]:
        z.write(a.run/name, name)
    for name in ["diagnostics.csv", "first_yield.txt"]:
        if (a.run/name).exists():
            z.write(a.run/name, name)
print(json.dumps(report, indent=2))
