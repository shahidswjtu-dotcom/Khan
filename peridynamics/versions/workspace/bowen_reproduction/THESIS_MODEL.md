# Integrated Bowen-method reconstruction

The user requested a renewed audit and a model matching Bowen. All27 modules in the preserved `bowen_pd` and new `thesis_pd` engines were parsed; the full existing and new verification suite was executed. Original Stage19 files and the old engine remain unchanged. Results below distinguish implemented equations from unprovided original definitions. This package must not be labelled the author's exact solver.

## Runnable model

```bash
cd /workspace/bowen_reproduction
.venv/bin/python run_thesis.py --out outputs/my_thesis_run --steps 5000
# compatible restart; steps is an absolute final target
.venv/bin/python run_thesis.py --out outputs/my_thesis_continue --steps 10000 --resume outputs/my_thesis_run/final.npz
.venv/bin/python -m pytest -q
```

`thesis_pd` integrates:

- Chapter3 geometry,5mm spacing,33,600 concrete points,0.015m horizon and3e-7s step.
- Eq2-12 weighted3D micro-modulus, Eq2-5 bilinear concrete tension envelope, invariant base thresholds, theory Eq2-13–15 DIF and persistent complete rupture.
- Eq2-39 partial-volume weighting and Eq2-40 directional energy surface correction. Symmetric endpoint angular blending follows the supplied Fortran; unlike the prior run, correction multiplies actual pair forces and their energy.
- Fig3-1 reinforcement centre depth20mm and full3D122-sample horizon neighbourhood per rod node, trilinear displacement and transpose force transfer.
- Eq2-45–47 particle projection and reaction from physical impulse/dt, with free axial hammer deceleration.
- Eq2-48–49 concrete short-range repulsion, with current-coordinate neighbour search including initially distant particles and broken bonds. Repulsion potential appears separately in energy budgets. Search is rebuilt each force evaluation; no stale family or arbitrary contact stiffness is used.
- Per-step signed force, explicit displacement probes, support/guide impulse and kinetic changes, concrete/interface/rod/self-contact potentials, rate-potential exchange and rupture deletion.
- Atomic full-state checkpoints with source/topology hashes and cached forces. Last contact reaction/count/penetration are now also retained, fixing immediate post-load diagnostic observation equivalence.

The current run starts from undeformed geometry. Repulsive force is not retroactively applied to already penetrated earlier endpoints: those endpoints have artificial overlap potentials larger than initial impact energy under the selected law.

## Explicit remaining reconstruction choices

The benchmark values are not fitted. However, unavailable author definitions cannot be replaced by claiming exactness:

1. Interface `bcc,s01,sc1` are still the inherited **bulk-concrete-derived demo inputs**, and Vh=dx³, Vs=lumped steel volume. Eq2-34 supports the volume product; it does not identify all numeric author quadrature inputs.
2. Rate law follows theory Eq2-13–15 rather than the differing AppendixA formula. Signed tensile rate, DIF_E multiplication of stiffness, and persistent rupture remain choices; supplied Fortran visible force law ignores prior flags.
3. Self-contact uses constant static c=12E/(πδ⁴) with cs=1.5c. The thesis does not establish how a reference-distance weighted micro-modulus should be assigned to newly contacting pairs beyond the original horizon. This convention is declared, not silently extrapolated or tuned.
4. Rods remain small-strain elastic axial elements; supplied source has a different yield cap, while the benchmark's complete rod loading/unloading/rupture law is not supplied. Yield activation is to be audited from output, not automatically guessed.
5. Exact hook geometry, lateral cover, bilateral supports versus bearing lift-off, axial hammer guide and contact splitting, initial1mm gap, and experimental transducer remain unverified.

The original equivalent Chapter3 computation is therefore not uniquely identified by the available inputs. Neither a fitted peak nor software test success settles this. A completed transient, converged force/deflection, source-matched sensor and matching shear morphology are required for a reproduction claim.

## Re-audit findings

The full suite passed22 tests (JUnit:outputs/reaudit_tests.xml). Tests include full physical interface force/moment, mass units, literal versus irreversible history, energy derivative consistency, current-neighbour self-contact and exact restart continuation/observation after contact. Stage19's112 recorded samples and immutable checksums still pass.

The actual assembled surface-corrected strain patches have mean energy ratios1.013(x),0.979(y),0.991(z) relative to continuum energy; individual boundary nodes range0.714–1.291. Multiplying a node's calibration ratio by its original energy is exact by construction, but the subsequent symmetric angular bond blending does not force every boundary node to match. This distinction is quantified rather than hidden by a trivial patch assertion.

Reference-direction linear rods produce about7067N under a0.05rad rigid rotation, demonstrating their finite-rotation limitation. A corotational rod would correct that mechanics limitation but would change the supplied linear-FEM convention; it has not been silently substituted in a model intended to reproduce the original author computation.

The matched-time0.24ms timestep comparison changes raw peak146.071→196.597kN (Δt halved), impulse6.369→5.838kg·m/s, bottom-layer displacement0.14134→0.12722mm, broken bonds32497→31166. Energy residual decreases42.302→19.481J. Thus the early response is not timestep-converged. Reaction is computed as impulse/dt, with no smoothing/filter fitted to obtain the target. Detailed values:outputs/thesis_timestep_comparison.json and outputs/assembled_methods_audit.json.

## Detached-flight diagnostic continuation

The inherited all-particle50mm trigger stopped at4111 steps (1.2333ms). The triggering particle, originally(.2725,.0375,.0975)m, had no surviving nonzero-stiffness PD bonds and lay outside the reference steel-coupling region. Its kinetic energy was0.937J; largest still-PD-attached displacement was28.56mm, steel2.93mm, and diagnostics remained finite. This trigger is therefore evidence of unbound fragment travel, not sufficient evidence of bulk instability.

`continue_thesis.py` resumes the unchanged full-state checkpoint and tracks detached flight separately, retaining the50mm trigger for concrete with surviving nonzero-stiffness PD bonds and for steel. It does not increase the threshold, modify mass/stiffness/damage/contact or select a timestep to fit peaks. The driver hash and diagnostic-policy change are explicit in its run manifest; source/model hashes and core force trajectory remain unchanged. The inherited stopped run remains preserved. This revised diagnostic is not proof of physical correctness; energy and timestep checks still apply. Results of the continuation will be recorded separately.

## Final executed interval and source-history alternative

Continuation completed5000 steps (1.5ms) without a further stop. The main observed force peak is146.071kN at0.1956ms; bottom-centre-layer support-relative deflection reaches4.07444mm. These do not reproduce359.8kN or22.70mm. The force and deflection are observations within this interval, not proven completed-event maxima, and the probe remains unmatched to the experimental sensor. Final144,879 concrete bonds and2830 interface samples are flagged broken. Envelope energy residual46.795J and momentum residual1.41e-12kg·m/s coexist: momentum conservation does not establish energy accuracy or source identity.

Final largest particle displacement71.089mm includes detached travel; largest displacement among nodes with surviving nonzero-stiffness PD bonds41.961mm; steel4.363mm. Final linear axial rod strain remains below the given yield strain, so missing plasticity was not activated under this specific linear-strain formulation. Finite-rotation validity remains a separate concern. The inspected reference-coordinate damage map shows central damage extending to the bottom with a broad upper crushed region; this alone does not establish the author's inclined shear crack opening/morphology. Figure:data in outputs/thesis_complete5000.

`printed_pd` and `run_printed.py` are an additional **printed force-history** reconstruction, not a claim of complete original-program identity. Unlike irreversible force gating, the active concrete/interface kernels compute force from current stretch without consulting prior rupture flags, retaining persistent flags only for damage reporting. The corresponding capped envelope potential remains stored beyond the threshold instead of falsely assigning irreversible fracture release to a source that reloads. Two new tests verify prior-flag reloading and persistent damage bookkeeping. The350-step contact smoke run reaches first contact341 and completes; its full impact response is unrun.

```bash
.venv/bin/python run_printed.py --out outputs/my_printed_history --steps 550
```

All other ambiguities remain explicit for this alternative: theory-vs-AppendixA rate formula, DIF stiffness, interface initialization and physical volume convention, supports/rod geometry and sensor. It is not an automatically preferred model or a parameter-fit branch. The full24 tests were covered as22 re-audit tests plus2 new printed-history tests; JUnit reports retain each executed run separately.
