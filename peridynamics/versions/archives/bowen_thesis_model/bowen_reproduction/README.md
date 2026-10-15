# Bowen RC impact audit and reconstruction

This is a research audit and a restartable diagnostic engine, **not a verified reproduction of Bowen's Chapter 3 results**. See `SOURCE_TRACEABILITY.md` and `REPRODUCTION_REPORT.md`. No material constants were adjusted to fit the target peaks. Original notebooks are preserved under `sources/`; Stage 19 baseline and hashes are under `baselines/stage19/`.

## Run

Python 3.12, NumPy/SciPy, Numba, matplotlib, pytest. The working environment is `.venv`:

```bash
cd /workspace/bowen_reproduction
python3 -m venv --system-site-packages .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pytest -q
.venv/bin/python run.py --out outputs/new_legacy550 --steps 550 --mode legacy
.venv/bin/python run.py --out outputs/new_family --steps 5000 --mode family_projection
# Continue the exact same code/configuration, using an absolute final step target:
.venv/bin/python run.py --out outputs/new_family_continue --steps 6000 --mode family_projection --resume outputs/new_family/final.npz
MPLCONFIGDIR=/tmp/bowen-mpl XDG_CACHE_HOME=/tmp/bowen-cache .venv/bin/python compare.py
.venv/bin/python diagnose.py
```

Run directories should be new; do not overwrite a run under review. `steps` is the absolute final step, not an additional count. Default checkpoints every250 steps are atomic, no pickle, and contain full u/v, hammer position/speed, concrete/interface damage and prior stretch, cached accelerations via force arrays, step, impulse, work and energy histories, source/model hashes. Restart rejects changes to code, timestep, mode or topology. Forces must be restored, not recalculated using committed stretch (which would erase the rate). Checkpoints are tied to the engine revision used for their run. The first four investigations used archived `baselines/engine_v1`; the new energy-enabled engine deliberately rejects these old checkpoints. Their original code is retained for continuation.

## Modes and units

- `legacy`: Stage19 active kernel and Verlet integrator; two lateral interface samples, 10mm bar depth, provisional 2e8N/m per-particle penalty contact. Reproduces original sampled CSV within rounding.
- `projection`: changes bar centre depth to20mm from Fig3-1, replaces penalty by Eqs2-45–47 projection contact, uses thesis explicit velocity/displacement order. Retains old interface topology for isolating changes.
- `family_projection`: additionally reconstructs Console2's 7³ lattice trimmed to nonzero horizon offsets and uses trilinear displacement/transpose force transfer. Full 3D family removes rank-one coupling defect. **It does not recover author interface calibration**; virtual dx³ volume, steel lumped volume, bulk-concrete-derived interface coefficients remain explicit assumptions. No automatic normalization by sample count.
- `AppendixCPublishedKernel` supports elementary `literal` and `irreversible_corrected` comparisons outside the dynamics engine. Its output has source density interpretation, not N. It is never divided by physical nodal mass as N.

Coordinates/displacements m; velocities m/s; physical forces N; masses kg; density kg/m³; micro-modulus N/m⁶; volume m³; energy J. Reaction comes from hammer force or impulse/dt, not a concrete stress proxy. CSV reactions retain sign and every timestep; no positive clipping or sparse peak sampling.

## Diagnostics and limitations

Named probes: reference bottom centre-layer(.300,.035,.0025)m; neutral axis(.300,.035,.050)m; top centre-layer(.300,.035,.0975)m; bottom displacement relative to interpolated support layers; legacy upper-half cross-section average. The original transducer definition remains unavailable; no extrapolation to z=0 or claim of sensor equivalence.

Physical momentum accounts for support impulses and the axial hammer guide's horizontal reaction. Work budget uses trapezoidal internal-force work, actual constraint/projection kinetic-energy changes and penalty potential. New energy diagnostics integrate the actual bilinear envelope: concrete/interface potential, linear rod potential, rupture potential deletion and rate-induced potential exchange are separate. These are **algorithmic accounting**, not a validated partially damaged unloading model or a proof that the discrete mesh dissipates prescribed Gf.

A50mm local displacement or dx/2 penalty penetration triggers a saved diagnostic stop. This is not automatically an instability or physical failure criterion. Original outputs include a fully detached fast particle; diagnose before extending the run or relaxing a threshold. No gravity is included, consistent with initialized impact speed; supports are bilateral constrained particle patches with horizontal anchors, not recovered bearing contacts. Missing thesis surface correction, concrete short-range self-contact and rod yield/rupture preclude predictive claims. Mesh refinement and complete response remain outstanding.

The cloud task is already isolated: use this workspace and existing checkouts; do not create a Git worktree unless the user explicitly requests it. `Khan` remains empty and unchanged; this project lives outside its checkout.
