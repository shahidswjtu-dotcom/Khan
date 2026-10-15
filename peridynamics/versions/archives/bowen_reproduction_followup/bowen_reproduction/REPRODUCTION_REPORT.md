# Bowen Chapter 3 reproduction investigation

**Status: not reproduced.** Stage19 regression verified; modular engine and separately labelled scientific corrections implemented. No material parameters tuned. Diagnostic extended runs have not established full impact/deflection/crack agreement. See `SOURCE_TRACEABILITY.md` for source equation and page references.

## Verified software results

Original Stage19 notebook, CSV, summary and plot were copied with SHA-256 checksums. The active path is cell125 → `stage18_run` → cell94 concrete kernel + cell118 conservative PDROD kernel, not the obsolete notebook functions. Seven Stage13–19 notebooks parsed without syntax errors. Baseline regression:112 stored samples, peak208.8390279446 vs208.8390279375kN, maximum absolute force difference3.46e-8kN; bond damage counts match exactly. First contact step341; at550:6252 failed bonds, max particle1.58403mm, max steel0.00221841mm. All-step recording finds244.731811kN at0.1506ms; sparse sampling missed this spike. Neither number establishes the final peak.

The original midspan metric is0.073428mm at550; reference bottom-centre-layer is0.041113mm. These must not be interchanged or compared as a matched experimental sensor.

Modules separate geometry/network, concrete, rod, PDROD, contact, integrator, probes, energy and atomic checkpoint IO. Unit tests verify geometry/masses, force units, rigid translation, concrete momentum, interface force/moment, interpolation and linear-field probes, rupture history, contact potential/impulse, energy-law derivatives, full-rank reconstructed family and restart. Checkpoints restore cached force history to preserve nonzero rate; code/model/timestep mismatches rejected. Versioned dependencies are in requirements.txt. **16 tests passed** in the final run; JUnit results:outputs/tests.xml.

## Runs and limitations

Machine-readable values are in outputs/comparison.csv; plots in outputs/comparison.png and outputs/comparison.pdf. All peaks/deflections below are **observed within truncated diagnostic intervals**, not converged final-response predictions. The bottom-centre-layer probe has not been matched to the transducer. Percent differences in the CSV are arithmetic comparisons of these partial responses, not validation errors of a completed impact.

| Run | Last time(ms) | Observed force max(kN) | Max bottom probe(mm) | Outcome |
|---|---:|---:|---:|---|
| Legacy Stage19 regression | .1650 | 244.732(all steps);208.839(original samples) | .04111 |550 steps completed |
| Legacy extension | .4029 |244.732 including initial550 | .24183 |50mm local-displacement diagnostic stop at1343 |
| Legacy Δt/2 | .4500 |243.203 | .26180 |3000-step diagnostic interval completed |
| Bar depth + projection contact, two samples | .4929 |649.216 | .42242 |50mm diagnostic stop at1643 |
| Full3D-family projection reconstruction | .4938 |720.025 | .40574 |50mm diagnostic stop at1646 |

Reference Table3-1: Bowen359.8kN,22.70mm; experiment301.0kN,23.38mm. Bowen's own relative force error is+19.535%; deflection−2.908%. The delivered calculations do not reproduce these targets or the complete punching/shear pattern.

The original legacy extension resumes its550-step checkpoint (the peak in that extension's summary applies to the continuation only; comparison.py merges the initial history). No live process is assumed to survive cloud publication. All four initial runs use archived baselines/engine_v1; subsequent energy-enabled code has its own hashes and intentionally rejects v1 checkpoints.

## Diagnosis

At the legacy stop, the largest-motion point is initially(.3025,.0325,.0875)m, fully detached (damage1), and has speed≈208.58m/s. Current-coordinate search finds750 pairs closer than the thesis's short-range threshold; smallest pair distance≈0.120mm vs5mm initial spacing. The projection alternative similarly has1227 overlapping pairs and detached speed≈159.79m/s. This proves the missing self-contact is active and material interpenetration is occurring. It does **not** by itself prove that all high fragment speeds are integration instability.

Legacy force peaks atΔt andΔt/2 differ≈0.625%, but damage trajectories, steel motion and local displacement differ substantially. A stable initial elastic step bound is not nonlinear convergence. Conservative undeformed stiffness bounds in outputs/stability.json give a combined sufficient Verlet step≈3.80e-7s for the two-sample model, before rate enhancement and compression geometric stiffness. This is not proof that3e-7s is safe for severe impact or explicit projection.

Legacy momentum residual≈1.5e-12kg·m/s, but small momentum residual does not prove constitutive fidelity. Its work residual at the stop≈1.019J; projection work residual≈46.779J. Projection splitting, constitutive rate changes, rupture and support constraints are separated in the enhanced energy logs. No dissipated fracture energy is inferred from counts. Dynamic instantaneous envelope potentials and rupture deletion are not a recovered author's unloading model.

## Verified corrections and outstanding scientific requirements

1. Physical PDROD force remains N=law·V_h·V_s, with transpose interpolation and kg masses. This corrects the earlier density/mass mistake and preserves action–reaction. Source Eq2-34 supports volume-product structure; the numeric integration volumes and interface calibration are still assumed.
2. Named bottom/top/neutral/support-relative probes and every-step signed reaction replace ambiguous measurement and peak sampling.
3. Fig3-1's80mm depth from top gives20mm from bottom; corrected modes use this without overwriting legacy geometry. Lateral cover and hook lengths remain unknown.
4. Eqs2-45–47 projection/impulse contact replaces penalty in corrected modes. No stiffness parameter fitted. Axial hammer and recoil timing remain reconstruction assumptions.
5. Console2's full3D horizon sample convention is reconstructed separately; the original two lateral samples have rank-one stiffness and no first-order vertical/longitudinal load transfer. Coefficients are still labelled demo; sample-count normalization is not fitted.

A faithful full run remains blocked by unavailable benchmark initialization (`bcc,s01,sc1,vol`, exact volumes/family, supports and loading/probe definitions) and unresolved source inconsistencies in rate law/integration. Missing surface energy correction, concrete self-contact and rod yield/rupture must be implemented and verified before predictive validation. Mesh refinement and complete response remain unrun. Supplying the author's Chapter3 initialization/full solver and the primary Fu experimental measurement definition would resolve key ambiguities; do not replace them by fitting target peaks.

## Environment persistence

Pinned install_script and startup/test/run instructions saved to the cloud configuration draft. No services or secrets required. The running machine and scientific outputs are prepared; saving the draft does not publish a snapshot. Review/save environment settings, then publish to retain this filesystem for later tasks. Fresh-task restoration has not been independently verified. The empty Khan repository remains unchanged.

## Final full-family and energy audit

The `family_projection` run used27,572 virtual interactions (122 per steel node), with the same demo coefficients as Stage19. It reached1646 steps, not the requested5000, and stopped at50.020mm local motion. Maximum steel displacement0.4100mm and69 broken interface samples show stronger vertical load transfer than two samples; this is not calibration evidence. Reaction peak720.025kN exceeds Bowen by100.118% and the experiment by139.211%; bottom-layer maximum0.405735mm is far below the reported maximum, but the run is incomplete and probe matching unresolved. Its momentum residual5.3e-13kg·m/s is small while envelope energy budget residual42.358J is not negligible. The instantaneous envelope contains10.523J concrete,0.0721J interface and0.1601J rod potential; cumulative rupture deletion14.535J and rate-potential exchange78.353J are distinct. Projection changes total kinetic energy by−83.156J. These terms expose the energetically consequential rate/splitting choices rather than hiding them in a fitted parameter.

The enhanced-energy legacy550 rerun preserves the old trajectory exactly:244.731811kN all-step peak,6252 broken bonds. Its current concrete envelope13.826J, rupture deletion4.428J and cumulative rate-potential exchange53.494J demonstrate why damage count cannot be called fracture energy. Envelope energy budget residual−0.555J (initial hammer energy2239.488J). The bilinear envelope still lacks author-defined unloading history; thermodynamic predictive validity remains unproven.

No completed full-impact run is delivered. A reproducible diagnostic stop, overlap evidence and quantitative conservation/energy results are delivered instead. A50mm fragment displacement is not enough to identify a unique instability cause, and this report does not claim one. Continuing an unverified model to an arbitrary later time would not establish scientific reproduction. Required next work is source resolution and verification of surface/self-contact/rod history, followed by matched-time timestep convergence and full response.

The final family damage map, inspected in reference coordinates, shows concentrated upper-midspan crushing and localized damage below it; it does not show the full through-depth inclined shear/punching pattern reported by Bowen. Thus crack morphology also fails the reproduction criterion at the saved diagnostic endpoint.

## Follow-up on exact reproduction (available files only)

The user clarified that no original author solver/input files are available. The goal remains author-result reproduction, but exact equivalence is not uniquely verifiable from the supplied incomplete specification. `EXACT_REPRODUCTION_REQUIREMENTS.md` defines the unavailable details and the distinction between rounded reported values and raw-run/bitwise comparison.

Two missing source methods were implemented independently in `source_methods.py`, with4 additional tests passing (`outputs/source_methods_tests.xml`). Directional surface-energy correction follows Eq2-40 and the supplied Fortran angular blending. Short-range repulsion follows Eq2-48–49 with all-current-neighbour search and an explicitly declared constant micro-modulus convention. They have not been inserted into the prior runs, and **no new improved impact result is claimed**. Quantitative endpoint audit:outputs/missing_methods_audit.json. The unchanged core engine remains compatible with its existing final checkpoints.
