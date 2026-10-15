# Fresh Bowen Chapter 3 reconstruction

This implementation starts independently from Bowen's documented equations. It does not import Stage 13–19 engines or use their numerical outputs as a regression target. Uploaded originals and previous directories remain untouched.

`pdbeam` contains the mesh, unique bond families, physical-force concrete operator, dynamic thresholds, surface correction, trilinear steel–concrete coupling, corotational elastoplastic rods, spherical hammer projection/recoil, current-neighbour self-contact, supports, explicit time integration, diagnostics and source/configuration-checked checkpoint/restart.

**The interface and several geometry/algorithm choices are reconstructed hypotheses. This is not yet established as Bowen's original executable model.** The [equation and assumption audit](AUDIT.md) specifies what comes from the thesis and what remains unresolved. No parameter is fitted to the target peaks.

Use the installed environment:

```bash
cd /workspace/bowen_from_scratch
PYTHONPATH=. /workspace/bowen_reproduction/.venv/bin/python -m pytest -q
/workspace/bowen_reproduction/.venv/bin/python build.py
/workspace/bowen_reproduction/.venv/bin/python run.py --steps 5000 --out outputs/hooked_dt3e7 --accept-bulk-interface-hypothesis
/workspace/bowen_reproduction/.venv/bin/python run.py --steps 4000 --dt 1.5e-7 --out outputs/hooked_dt1p5e7 --accept-bulk-interface-hypothesis
MPLCONFIGDIR=/tmp/bowen-mpl /workspace/bowen_reproduction/.venv/bin/python analyze.py
```

Fresh installation: create a virtual environment, install `requirements.txt`, and substitute that Python executable. Only runtime packages are reused from the previous workspace; no previous simulation engine is imported.

For continuation, choose a new output directory so existing histories are retained:

```bash
/workspace/bowen_reproduction/.venv/bin/python run.py --steps 10000 --restart outputs/hooked_dt3e7/checkpoint.npz --out outputs/continuation_3ms --accept-bulk-interface-hypothesis
```

`--steps` is the total target step number, not the number of additional steps. Restart requires identical package source hashes, material/geometry values and declared reconstruction choices. Checkpoint histories, hammer state, plastic strain and cumulative diagnostics are restored. A segment's summary peaks cover that segment only; combine histories before reporting full-run maxima.

The initial sphere bottom lies on the geometric beam top, while concrete particles are half a cell below it. The resulting free-flight delay must be distinguished from time since first contact. The bottom-midspan probe is explicitly labelled and is not asserted to equal Bowen's measurement. The 1.5 ms run does not establish the full-event deflection maximum.

`audit_projection_energy.py` measures contact/constraint changes in stored potential separately from the base solver's kinetic-energy accounting, without changing mechanics. A checkpoint-started audit covers that continuation segment only. Its augmented residual retains the prior residual offset; compare changes against the residual at the starting checkpoint. See generated reports before inferring conservation from any single residual column.

The active reconstruction includes the end hooks shown in Bowen's Fig. 3-2, with an explicitly assumed 60 mm rise. Earlier straight-bar results and their exact source are kept under `outputs/impact_*` and `versions/fresh1`, respectively. Atomic-checkpoint correction without topology change is pinned under `versions/fresh1_atomic`. These are versions of the fresh implementation, with no Stage 19 dependency. Do not mix their checkpoints with the current hooked configuration.

Read `RESULTS.md` for completed comparisons and the current reproduction status. The latest 22-test suite passes; Bowen's benchmark has not been reproduced. The compact delivery includes the irreversible final checkpoint and both completed 1.5 ms histories. Generate the reference mesh with `build.py` before rerunning plots. Outputs must use a new directory to preserve complete histories.
