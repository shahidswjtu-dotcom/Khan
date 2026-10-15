# CODEX HANDOFF — Tu Bowen RC Beam Impact Peridynamics Reproduction

**Prepared:** 2026-10-09  
**Project:** Reproduce the **normal-concrete reinforced-concrete beam impact simulation** reported in Tu Bowen's 2025 master's thesis, *基于近场动力学方法的混凝土梁构件冲击损伤研究* (Research on Impact Damage of Concrete Beams Based on Peridynamics).  
**Audience:** Codex or a research software engineer implementing, auditing, testing, and documenting a scientific Python/Numba simulation.  
**Critical rule:** Do **not** claim that Bowen's output has been reproduced until an independently executed, complete simulation with source-supported parameters agrees in force, deflection and crack morphology. Do not tune arbitrary factors to force agreement.

---

## 0. Paste this instruction into Codex

> You are the computational mechanics engineer continuing this Peridynamics project. Read `CODEX_HANDOFF_BOWEN_PD_REPRODUCTION.md` in full, inspect `Bowens Thesis.pdf` (especially Chapters 2–3 and Appendices A–C), inspect `fortran_sources/Console2.f90`, and inspect `Bowen_RC_Chapter3_Stage19_Extended_Impact.ipynb` plus Stages 13–18. Your task is to develop a verified, modular, restartable Python/Numba implementation reproducing Bowen's Chapter 3 normal-concrete RC beam drop-hammer impact, while keeping the **literal printed Appendix C behavior**, **our reconstruction assumptions**, and **physically consistent corrections** separate. Start with a source-to-code matrix and explicit parameter provenance. Do not invent missing Fortran initialization, bond/interface quadrature, reinforcement layout or contact algorithm. If missing, list the exact unresolved definitions and build a transparent configurable alternative. Prioritize unit tests for force units, conservation, irreversible damage, steel nodal mass, contact, energy, and probe selection; then implement checkpoint/resume, execute full-impact runs and compare with 359.8 kN, 22.70 mm, and cracking. Maintain reproducible scripts, tests, run manifests, and figures; describe every material numerical deviation from Bowen.

---

## 1. Scientific objective and definition of success

**Primary objective:** Follow Tu Bowen's computational methodology for a 3D bond-based peridynamic (BBPD) reinforced-concrete beam without changing the benchmark geometry/material inputs for convenience. Implement 1D reinforcement elements and eight-node PDROD steel–concrete interpolation, dynamic concrete damage, support conditions, and rigid hammer impact. Obtain time histories and spatial damage that can be compared with Chapter 3.

**Published benchmark targets** (from prior thesis study; reconfirm figure and table labels directly in PDF):

| Measure | Bowen PD | Experimental comparator |
|---|---:|---:|
| Peak impact contact force | 359.8 kN | 301.0 kN |
| Maximum midspan deflection | 22.70 mm | 23.38 mm |
| Failure morphology | Local crushing/punching and inclined shear cracking | Experimental crack pattern |

Matching two peaks alone is insufficient. Minimum evidence: (1) source traceability, (2) mesh/boundary/contact checks, (3) peak AND force–time/deflection–time response, (4) spatial/time-resolved cracking, (5) energy and momentum audits, (6) validation against experimental reference independently of model fitting. Document the approximately 19.5% impact-force error in Bowen's own model rather than presenting it as perfect agreement.

**Fidelity tiers:**
- **Tier A — printed source equivalence:** match exactly the operations visible in Appendix C / supplied Fortran for a controlled single interaction.
- **Tier B — thesis-method reconstruction:** reproduce documented equations and Chapter 3 setup, documenting assumed details.
- **Tier C — predictive validation:** full dynamics, convergence and agreement with experimental measurements.
These tiers must never be conflated.

## 2. Source inventory (start here)

Core local files available in the provided handoff bundle or user's uploaded workspace:

| File | Description/role |
|---|---|
| `Bowens Thesis.pdf` | Primary scientific specification. Read Ch. 2 theory; Ch. 3 RC impact benchmark; Appendices A, B, C; figures/tables. Chinese PDF, 109 pages. |
| `fortran_sources/Console2.f90` | Supplied 3D Fortran PD + steel interpolation/PDROD implementation (~1,870 original lines). **Different configuration** from Bowen's 600 mm Chapter 3 beam; useful for algorithmic conventions but NOT proven to be the exact benchmark source. |
| `Bowen_RC_Chapter3_Stage19_Extended_Impact.ipynb` | Latest integrated Python/Numba baseline. **128 accumulated cells**; many historical experimental functions and different pathways coexist. The active Stage 19 run uses `stage18_run(mode='conservative')`. |
| `Bowen_RC_Chapter3_Stage19_Executed.ipynb` | Executed variant with previous outputs. |
| `Bowen_RC_Chapter3_Stage18_Conservative_PDROD.ipynb` | Force-balanced interface implementation and checks. |
| `Bowen_AppendixC_Literal_Python_Stage13.ipynb` | Single-interface printed Appendix C literal translation, plus separately labeled irreversible version. |
| `Bowen_RC_Chapter3_Stage14_AppendixC_Beam_Integration.ipynb` | Interface integrated with beam framework; tests. |
| `Bowen_RC_Chapter3_Stage15_Dimensional_Stability_Audit.ipynb` | Diagnostic of interface stiffness and time-step instability. |
| `Bowen_RC_Chapter3_Stage16_Fortran_Force_Density_Audit.ipynb` | Documents important source force-density versus nodal-mass normalization issue. |
| `Bowen_RC_Chapter3_Stage17_Density_Normalized_Impact.ipynb` | Source-inspired acceleration normalization; stable short run, but physical force balance concern. |
| `bowen_stage19_outputs/extended_impact_history.csv` | Saved 550-step (every 5 steps plus first contact) numerical history. |
| `bowen_stage19_outputs/extended_impact_summary.json` | Summary of latest diagnostic run. |
| `bowen_stage19_outputs/extended_impact_diagnostics.png` | Force, midspan and damage diagnostic plots. |
| `bowen_stage15_outputs/appendix_c_dimensional_stability_audit.csv` | Dimensional/local time-step screening data. |

**Other background:** `Review paper about peridynamics.pdf` (Yu et al., 2025, *International Journal of Impact Engineering*, 206, 105457), and independent `Plain Concrete Slab.zip` / `RC Beam Model.zip`. These may provide context or benchmarks but do not override Bowen's thesis.

**Do not assume code path identity from a filename:** The latest notebook includes historical Stage 1–19 definitions. Confirm which function is *actually called* by the active Stage 19 cells before modifying anything.

## 3. Benchmark geometry, physics and numerical settings

Values used in the current Python `Parameters` dataclass (not every material/contact value is fully source-verified):

| Parameter | Current value | Provenance/qualification |
|---|---:|---|
| Beam length, width, depth | 0.600 × 0.070 × 0.100 m | Bowen Ch. 3 benchmark |
| Support span | 0.500 m | Bowen Ch. 3 |
| Concrete particle spacing | 0.005 m | Bowen Ch. 3 |
| Concrete grid dimensions | 120 × 14 × 20 | Computed, 33,600 particles |
| Horizon | 0.015 m = 3Δx | Source benchmark |
| Time step | 3.0×10^-7 s | Source benchmark |
| Concrete E | 25 GPa | Source benchmark |
| Concrete ν | 0.25 | BBPD 3D constraint |
| Concrete ρ | 2,400 kg/m³ | Source benchmark / current model |
| Concrete tensile strength | 4.6 MPa | Source benchmark |
| Concrete compressive strength | 60 MPa | Source benchmark |
| Concrete fracture energy | 105 N/m = 105 J/m² | Source benchmark |
| Steel E | 200 GPa | Current benchmark |
| Steel density | 7,750 kg/m³ | Current notebook; verify against thesis |
| Steel fy | 482 MPa | Benchmark |
| Longitudinal bars | 2 × Ø6 mm | Benchmark; **exact transverse placement/ends uncertain** |
| Hammer mass | 38.4 kg | Benchmark |
| Initial hammer speed | 10.8 m/s | Benchmark |
| Hemispherical hammer radius | 0.060 m | Benchmark |
| Initial impactor gap | **Provisional 1 mm**, geometry-based first contact about step 341 | Not proven to be the original initialization |
| Contact penalty stiffness/short-range law | **Provisional** | Not recovered exactly from source |
| Steel-to-concrete virtual bond offset/location | `COUPLING_OFFSET=dx`, two lateral samples per steel node | **Provisional**, not certified Bowen |
| Interface coefficients `bcc,s01,sc1,vol` | Current `DEMO_INTERFACE` derived from concrete values | **Demonstration only — NOT the original PDROD calibration** |

### Particle and interaction topology

Current concrete grid: **33,600 material points**, nominal **1,755,720 unique pair bonds** (Stage 3 network), neighbor counts ~28 corner to ~122 interior. Steel representation in Stage 4: **226 steel nodes / 224 rod elements / 452 coupling samples**. Confirm exact dimensions and quadrature points in the active cells; do not assume these match Bowen's original reinforcement discretization.

### Known static formulas from thesis (verify equation numbering)

` s0 = ft/E `; ` sc = 5 Gf/(6 E δ s0) `; distance-dependent micromodulus ` c(ξ,δ)=72 E/(πδ^4) (1-(ξ/δ)^2)^2 `; dynamic threshold shifts for ordinary concrete ` s0d=(DIF_t/DIF_E) s0 ` and ` scd=(DIF_t/DIF_E) sc `. Stage 7 also tested `c_d = DIF_E c` separately. The exact strain-rate definitions and Fortran branch logic must be checked against Appendix A. Avoid repeated compounding of thresholds (`s0 *= dif` across time steps); recompute from invariant base values.

A fracture-energy weighted integral consistency test passed for a specific idealized kernel, but **discrete fracture energy over a real finite particle family is not validated**.

### Appendix C steel–concrete coupling workflow

For each interaction sample: interpolate virtual concrete displacement via 8 trilinear shape functions; compute initial/current/previous bond length, bond stretch and rate; evaluate `bc=bcc*(1-(idist/delta)^2)^2`, piecewise bond law with `s01/sc1`; accumulate source-style `pforce33` on steel and weighted opposing contributions on 8 concrete particles; update `failfail` and local damage.

**Critical audit finding:** the *printed* Appendix C force calculation updates `failfail` after evaluating force and does **not visibly check prior `failfail`** when computing force. The literal mode can transmit force after previous failure/reloading. The separately introduced `irreversible_corrected` mode prevents this. The published excerpt alone is insufficient to determine the entire original program's history behavior. Keep both modes explicit and test them.

**Critical normalization finding:** `Console2.f90` updates accelerations from **force density / material density** (`pforce/dens`, `pforce/dengs`). An earlier Python port instead interpreted a force-density contribution as N and divided by steel nodal mass; this generated artificial steel acceleration. The alternate physically conservative Stage 18 implementation converts interface density to nodal force using an **assumed steel control volume**, which restores action–reaction but is *not source-verified*. Distinguish these as two modeling choices, not interchangeable identities.

## 4. Full development record: Stages 1–19

| Stage | Work completed | Main evidence / limitation |
|---|---|---|
| 1 | Reference geometry/particle count, materials, preliminary bar/support/hammer visualization | 33,600 grid confirmed; exact reinforcement coordinates provisional |
| 2 | Bond elastic/softening, rate functions, irreversible damage tests | Unit tests passed; no assembled solver |
| 3 | 3D neighbor network, internal equal-and-opposite pair forces, partial-volume/surface correction tests | ~1.756 M unique bonds; momentum tests pass; correction equivalence unverified |
| 4 | 1D steel rods, 8-node interpolation, prototype interface | 226 steel nodes, 452 samples; **linear penalty interface initially** |
| 5 | Rigid dynamic hammer, provisional contact, supports, velocity-Verlet | Short smoke test; source contact not exact |
| 6 | Verification/source-fidelity audit | Missing constants, rate-law and energy limitations catalogued |
| 7 | Appendix-C-shaped nonlinear PDROD law, rate-dependent micromodulus function | Standalone tests; nonlinear not yet injected into actual solver at that stage |
| 8 | Integrated nonlinear PDROD + concrete rate law into dynamics | Few-step smoke test only |
| 9 | Performance profiling / first-contact diagnostic | Original vectorized full-mesh 355-step attempt timed out |
| 10 | Numba accelerated concrete and PDROD loops | First contact step ~341; short force 22.56 kN by step 365; divergence upon extension |
| 11 | Dynamic instability audit | At step 371, concrete/steel motion explodes; interface damage initiates at ~370 |
| 12 | Concrete-only and weak elastic interface isolation | Both reach step 384; nonlinear coupling implicated |
| 13 | Standalone **literal Appendix C** and corrected irreversible comparison | Printed failure-history discrepancy documented |
| 14 | Literal Appendix C adapter integrated into model; separate options | Force balance and unit tests; provisional material coefficients remain |
| 15 | Dimensional and local stability audit | Provisional `bcc` ~1.132e19 and local estimated dt limit ~2.21 ns vs 300 ns, explaining potential blow-up; NOT a critique of Bowen's actual inputs |
| 16 | Source-level `Console2.f90` force-density vs nodal-mass audit | One-interaction acceleration scaling ratio ~1.46e5 under incorrect interpretation |
| 17 | Source-style density-normalized interface path | Stable 420 steps; ~74.37 kN early recorded max; physically unbalanced effective-volume interface |
| 18 | Force- and moment-conserving interface conversion with assumed steel control volume | Stable 425 steps; force ~74.37 kN; steel motion tiny; SOURCE EQUIVALENCE NOT established |
| **19** | **Extended conservative-interface impact run** | **550 steps, no safety stop, ~208.839 kN sampled max, ~6,252 concrete bonds broken; full response not reached** |

All reported checks are *past local diagnostics*, not an independently certified numerical solution. Re-execute tests in a clean environment before depending on them.

## 5. Current executable state: Stage 19

**Main file:** `Bowen_RC_Chapter3_Stage19_Extended_Impact.ipynb`.  
**Core configuration cell 124 (zero-based indexing):**

```python
RUN_STAGE19 = True
STAGE19_STEPS = 550
STAGE19_SAVE_EVERY = 5
```

**Active call in cell 125:**

```python
stage19_data, stage19_summary = stage18_run(
    steps=STAGE19_STEPS,
    dt=p.dt,
    save_every=STAGE19_SAVE_EVERY,
    mode='conservative',
    limit_mm=50,
    verbose=True,
)
```

The latest notebook contains **128 cells**. Notebook code cell IDs as noted here refer to **zero-based order in the JSON `cells` array**, which Codex can inspect programmatically.

Other important code locations in this notebook:
- Cell **3**: dataclass `Parameters` and benchmark values.
- Cell **42**: `concrete_grid_eight` and provisional interface offset/samples.
- Cell **58**: `apply_supports` and original `run_impact` implementation.
- Cell **74**: `AppendixCInterfaceParameters` / `DEMO_INTERFACE` (**not source verified**).
- Cell **75**: Python Appendix C force translation.
- Cell **94**: `_fast_concrete_kernel`.
- Cell **96**: `_fast_pdrod_kernel` (earlier convention).
- Cell **97**: optimized `fast_system_forces` / `run_fast_impact`.
- Cells **112–113**: source-style force-density normalization (`stage17_*`).
- Cell **118**: `_stage18_pdrod_conservative` (active interface kernel).
- Cell **120**: `stage18_run` (active integrator).
- Cells **124–126**: Stage 19 settings, execution, saved history/plots.

**Warning:** the Stage 19 interface force is multiplied by both `vol` and `steel_cv[sid]` inside cell 118; those are *not proven* to be Bowen's actual Appendix C force-quadrature factors. Merely being stable or conservative does not establish equivalence.

### Stage 19 numerical summary (source file `bowen_stage19_outputs/extended_impact_summary.json`)

```json
{
  "mode": "conservative",
  "steps_completed": 550,
  "first_contact": 341,
  "safety_stop": null,
  "max_recorded_force_kN": 208.83902793749814,
  "max_concrete_mm": 1.584030766835752,
  "max_steel_mm": 0.002218410738845336,
  "elapsed_seconds": 19.05,
  "interpretation": "diagnostic only; interface normalization and material constants not Bowen-verified"
}
```

Previous Stage 19 report also gave **max measured midspan magnitude ~0.0734 mm** and **6,252 failed concrete bonds** at step 550. Confirm these by reading history rather than assuming the monitoring probe matches Bowen's test sensor. Time simulated is 550×3e-7 s = **0.165 ms**. This is a *truncated transient*, not necessarily a final impact peak. First contact occurs about 0.1023 ms (depending on particle discretization and provisional gap). Saved `impact_kN` uses positive clipping of hammer reaction and sampling every fifth step; finer output may change a detected local max.

**Caution about saving/restart:** `stage18_run` initializes all states to zero and resets constitutive history every invocation. It writes CSV but **does not currently save the complete integration state** at configurable checkpoints. A longer run currently starts over, even if `STAGE19_STEPS` is increased. Crucial stored states to checkpoint include u/v of concrete and steel, hammer position and velocity, previous accelerations or their deterministic reconstruction, broken concrete bond flags, previous concrete bond stretch, interface broken flags, previous interface stretch, support/load state, current step/time, and config/source hashes. Verify any internal rate-law or damage state beyond those named.

## 6. High-priority outstanding correctness issues

### P0: Establish genuine parameter provenance, without inventing values

1. From thesis Ch. 2, Ch. 3 and Appendices A–C, compile every equation and parameter with original Chinese term, units, location/page, and Python equivalent.
2. Trace Fortran **declarations/assignments** of `bcc`, `s01`, `sc1`, `vol`, `pforce33`, `dens`, `dengs`, `failfail`, steel nodal/integration volume, rod connectivity and force sums. **Check initialization and use**, not merely variable names.
3. Check whether the supplied `Console2.f90` is actually from Bowen or a different project. It is previously observed to have a different geometry and prescribed moving cylindrical tools. Do **not** attribute its specific constants to the thesis benchmark without evidence.
4. Separate strict **literal Appendix C mode** from **irreversible corrected** and **momentum-conserving reconstructed** modes, with independent tests and clearly named run manifests.

### P0: Units and explicit-integration stability

- Define one consistent physical/numerical convention for every array: N, N/m³, kg, kg/m³, m, m³; use typed metadata or assertions.
- Prove conversion of steel–concrete force density into nodal force, if such a conversion is intended; identify precise integration/control volumes.
- Validate momentum conservation with **physical nodal forces**, not just equal-and-opposite force-density entries.
- Check mass distribution of bar elements and concrete discretization.
- Compute local stiffness eigenvalue or stable time step (PD bonds, steel rods, interface, contact); verify numerical stability independently of displacements being below an arbitrary cutoff.
- Check compression-side piecewise bond and failure history: softening only in tension? Examine Appendix C exactly. No spurious force after rupture in corrected mode.
- Unit-check `dt` in stretch-rate formula and rate-DIF branches; avoid repeated multiplication of critical stretches.

### P0: Midspan deflection probe is probably not equivalent to Bowen's

Current active `stage18_run` uses:

```python
mid=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1.e-12)&(xyz[:,2]>=p.H/2)
midspan_mm=-1e3*np.mean(uc[mid,2])
```

This averages **all particles in a central cross-sectional half-depth**, not the bottom midspan point or measurement at the location of the displacement transducer. Define explicit probes: bottom-center concrete displacement, top-center, neutral-axis, midspan cross-section mean, and support-relative deflection; use trilinear interpolation at exact sensor coordinates if necessary. Match Figure 3-4's sensor definition only after checking the thesis/experimental source. Include probe metadata in every CSV.

### P0: Contact and boundary conditions

- Confirm falling hammer hemisphere versus full-sphere effective representation and whether the hammer is rigid with full mass and **free deceleration** after contact.
- Recover exact first-contact position and surface force algorithm, contact stiffness, contact search threshold, correction of penetrating particles and reaction force extraction.
- Current initial gap/penalty law is **provisional**. Do not tune contact stiffness to target peak without independent calibration.
- Confirm simple-support positions, actual bearing/contact width, support displacement/rotation constraints, treatment of upward lift-off, horizontal anchor, and whether support particles are artificially clamped.

### P1: Reinforcement behavior and bond geometry

- Confirm steel endpoint positions, cover distance, number/spacing of rod nodes, 6 mm bars, endpoints and ties.
- Check 8-node stencil selection for each virtual concrete point; sample quadrature and interpolation partition-of-unity / moment checks.
- Resolve why Stage 19 steel displacement (~0.0022 mm) is much smaller than localized concrete displacement (~1.584 mm); could be realistic early transient, weak coupling, or missing force scaling.
- Compare steel and concrete acceleration, local relative slip, interface energy, concrete-to-steel force and support reactions in both elastic and damaged phases.
- Confirm rebar plasticity/yield/hardening and any rate response match thesis.

### P1: Validation, energy and numerical convergence

- Add energies: hammer KE, concrete/steel KE, recoverable bond energy, rod axial energy, interface stored energy, cumulative irreversible fracture energy, contact potential/dissipation, gravity/external work and support work. Check residual carefully; rate-dependent and path-dependent constitutive laws require explicit accounting.
- Log contact force at **every** integration step near first contact, peak and separation; record reaction definition/sign convention.
- Distinguish peak **particle displacement** from **support-relative midspan deflection**.
- Plot spatiotemporal damage: per-particle fraction of failed bonds; shear wedge/punching crack morphology; location of damage before and after peak load.
- Repeat with Δt/2 and grid refinement; when Δx changes, distinguish fixed-horizon from horizon/spacing-ratio refinement. The thesis's 4,5,6 mm comparisons vary horizon together with spacing.
- Include independent benchmark for elastic 3D BBPD, 1D rod, force quadrature, rate law, impactor conservation, and multi-body coupling before production runs.

### P1: Performance and software reproducibility

- Refactor a 128-cell cumulative research notebook into a documented Python package with separate modules. Do not alter the baseline without regression tests.
- Use Numba precompiled loops for ~1.76 million unique concrete bonds.
- Implement deterministic checkpoint/restart and verify **restart equivalence** versus a continuous run.
- Save partial results at least every fixed number of steps; use atomic checkpoint writes.
- Record Git commit/source hash, package versions, user inputs, random seeds, simulation wall time and hardware details.
- Make run and plotting scripts operate without manual editing of many notebook cells.

## 7. Execution plan for Codex (do in this order)

### Task 1 — Inventory and baseline preservation

1. Confirm expected files exist and inspect notebook cells via `nbformat`/JSON.
2. Run static notebook code checks and a minimal Stage 18 kernel sanity check (avoid repeatedly executing full 128-cell notebook unless explicitly intended).
3. Preserve Stage 19 source/history/summary under `baselines/stage19/` with checksums.
4. Write `SOURCE_TRACEABILITY.md` with exact PDF page and Fortran line ranges for each algorithm.

**Deliverable:** manifest + reproducible baseline and source-to-code comparison table.

### Task 2 — Minimal modular engine and tests

Recommended structure:

```
bowen_pd/
  README.md
  pyproject.toml
  config/reference_bowen.yaml
  config/assumptions.yaml
  src/bowen_pd/
    geometry.py
    neighborlist.py
    concrete.py
    rod.py
    pdrod_literal.py
    pdrod_conservative.py
    impact_contact.py
    integrator.py
    diagnostics.py
    io.py
  tests/
    test_geometry.py
    test_bond_law.py
    test_rate_dif.py
    test_pdrod_interpolation.py
    test_force_units_and_balance.py
    test_rod_mass_and_forces.py
    test_contact.py
    test_restart.py
  scripts/
    run_reference.py
    compare_bowen.py
    export_figures.py
  baselines/stage19/
  outputs/  # gitignored
```

**Deliverable:** unit-tested reusable code with 'literal', 'source_density', 'conservative' experimental modes controlled from config and per-mode provenance.

### Task 3 — Critical physics corrections and tests

1. Resolve or flag each missing `bcc/s01/sc1/vol` variable and nodal mass normalization; never silently supply unrelated bulk-concrete values as steel–concrete calibration.
2. Prove interface force/moment conservation in physical units.
3. Check time-step limits and correct contact/support loading.
4. Correct and document sensor probes.
5. Verify source-history vs irreversible-history options. Do not default to an unmarked correction.

**Deliverable:** test reports, dimensional map, scientific gap list; no full validation claim yet.

### Task 4 — Checkpointed full transient

1. Reproduce first contact (~341 steps in current provisional geometry) and Stage 19 550-step baseline under identical legacy settings as a regression check.
2. Advance beyond 550 to cover rising contact force, first local maximum, unloading and major deflection, with periodic checkpoint saves and restarts.
3. Monitor at every step for NaN, energy residual, contact penetration, velocity/acceleration spikes, unsupported negative energy, abrupt damage reversal, and physically impossible displacements.
4. Save all time histories and spatial snapshots in HDF5/NPZ; one immutable JSON manifest per run.

**Deliverable:** completed reference transient **or a reproducible failure report pinpointing the exact root cause**.

### Task 5 — Quantitative comparisons and independent validation

Produce figures/table with: peak contact force and timing, midspan deflection peak and timing, force–time and deflection–time curve overlays with digitized experimental/PD data where legally obtained, crack/damage maps at matching instants, energy residual, and mesh/time-step sensitivity. Compute clear percent errors with sign; distinguish predicted model vs experiment.

**Deliverable:** `REPRODUCTION_REPORT.md` + publication-quality source-traceable figures and machine-readable CSV.

## 8. Acceptance checklist (must be demonstrated, not asserted)

- [ ] Exact beam geometry, mesh count and bar/support/hammer locations verified against thesis.
- [ ] Reconstructed parameter provenance includes all interface values and loading constants or explicitly labels those irrecoverable.
- [ ] Code equivalence to Appendix C on deterministic elementary force cases; original broken-bond flag behavior separately tested.
- [ ] Concrete pair-force momentum conservation; steel–concrete physical force and moment conservation.
- [ ] Rigid translation and constant-strain patch tests; correct dynamic DIF and fracture-energy tests.
- [ ] Hammer free flight, impactor/contact Newton's-third-law and energy tests.
- [ ] Explicit stability demonstrated under mesh-specific stiffness estimates and dt refinement.
- [ ] Checkpoint/restart gives the same trajectory within floating-point tolerance.
- [ ] Full impact simulation completes with physical reaction, stable energy diagnostics and meaningful cracks.
- [ ] Displacement definition matches Bowen measurement point, not the current half-depth average.
- [ ] Results compared to Bowen's 359.8 kN and 22.70 mm **and** experiments (301.0 kN, 23.38 mm).
- [ ] Parameter tuning choices are predeclared and separated from verification/validation data.
- [ ] All deviations from Bowen's published method are documented and visible in plots/report.

## 9. Known red flags and things NOT to do

1. **Do not** call Stage 19's 208.839 kN sample a fully converged impact peak; it is only a 0.165 ms early response.
2. **Do not** compare 0.0734 mm 'midspan' to 22.70 mm without matching probe location/definition.
3. **Do not** claim `Console2.f90` is the author's exact Chapter 3 solver solely because it contains PDROD-like operations; its original geometry/load differ.
4. **Do not** hide missing interface inputs with `DEMO_INTERFACE` while labeling the simulation 'exact'.
5. **Do not** fix instability by arbitrarily reducing stiffness, increasing displacement thresholds, turning off rupture, or modifying mass without physical justification and documentation.
6. **Do not** conflate `pforce` **N/m³** with physical force **N**; always check masses and volumes.
7. **Do not** count array-value equality for density forces as full physical momentum conservation.
8. **Do not** trust passing single-interaction tests as proof of correct full nonlinear dynamics.
9. **Do not** infer fracture energy from the global broken-bond fraction; calculate actual dissipated energy.
10. **Do not** assume surviving 550 steps proves the time step is stable across severe post-peak fracture.

## 10. Desired response from Codex after starting

Codex should produce:
1. A short source inventory listing every file actually read and any missing files.
2. A table mapping thesis equation / Appendix C Fortran statement / current Python function / unit / fidelity status.
3. A ranked list of *proven bugs* vs *hypotheses* vs *missing-source ambiguities*.
4. The first concrete code patch with tests and an executable command to reproduce results.
5. Clear next run instructions and computational budget.

**Final guiding principle:** Accurate, auditable reproduction of documented science is more important than a visually plausible force–time curve or fitting the two reported peak values.
