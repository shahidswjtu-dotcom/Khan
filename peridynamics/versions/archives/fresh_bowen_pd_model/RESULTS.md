# Reproduction assessment

**Bowen's original Chapter 3 results have not been reproduced.** Numerical component verification passes, but the available implementation definitions do not identify the original Chapter 3 executable uniquely. Neither tested force-history interpretation matches the published force peak or establishes the published maximum deflection and crack pattern.

## Fixed-window results

Both fresh hooked models ran to 5,000 steps (1.5 ms) at the stated 0.30 µs timestep. Material constants were held at the published values. The interface and undimensioned geometry choices are explicitly recorded hypotheses in each manifest.

| Quantity | Bowen | Experiment cited by Bowen | Irreversible reconstruction | Printed-force-history interpretation |
|---|---:|---:|---:|---:|
| Peak impact force, kN | 359.8 | 301.0 | 220.663 | 252.045 |
| Maximum deflection / observed bottom-midspan displacement, mm | 22.70 | 23.38 | 16.547 | 16.691 |
| Global damage index | 0.1701 | — | 0.30381 | 0.31811 |

The computed displacement columns are **observed window maxima of the bottom-layer midspan mean**, not established complete-event maxima or a matched original transducer. Deflection was still increasing at the window end. The original damage-index reporting time is not fully identified; the computed indices are at 1.5 ms and follow Eq. 2-10's family-volume definition. These qualifications prevent interpreting the table as a successful validation.

The observed force peaks are respectively 38.67% and 29.95% below the thesis value. No material parameter was optimized or rescaled against the targets. The history alternatives were selected to investigate the visible omission of prior failure flags in the printed force calculation, not to fit the peaks.

## Numerical evidence

- All **22 tests pass** in the newly installed Khan environment. Checks cover force units, momentum and torque, rigid motion, interpolation, yield history, force/energy gradients, rate and rupture accounting, contact, full assembled checkpoint continuation and atomic-write failure recovery.
- Initial elastic stability bound: sufficient timestep about 0.868 µs; used 0.300 µs. This bound excludes evolving nonlinear self-contact.
- Irreversible branch, matched 0.6 ms comparison: halving dt changes impulse by about 0.157%, bottom-midspan displacement by about 0.116%, and the isolated raw force peak by **8.485%**. Raw peak convergence remains unresolved. The printed-history branch has not had an equivalent half-timestep study; its cross-variant plot is labelled accordingly.
- Maximum absolute closed-system material energy drift over the sampled histories: approximately 2.03 J (irreversible) and 2.59 J (printed history), against initial hammer energy 2,239.488 J. Maximum momentum residuals are below 9×10⁻¹⁰ kg m/s.
- Contact/constraint kinetic substep changes and projection potential changes are separate diagnostics. The approximately 808/817 J kinetic-substep ledger residuals are not the closed-system energy drift; interpreting them that way double-counts part of the internal contact process.
- All 28 bottom-midspan probe particles remain in supported components at the final state under each branch's actual force-history rule. The irreversible run has detached fragments whose motion is much larger than the probe displacement. In the printed branch, previous failure flags do not prevent restored force transmission; its ever-failed damage map therefore cannot by itself establish permanent crack separation.

## Source-supported corrections

The fresh package includes physical pair-force normalization; consistent lumped rod mass; conservative, torque-balanced trilinear steel–concrete force transfer; three-dimensional coupling samples; the documented 20 mm longitudinal bar height; vertical end hooks visible in Fig. 3-2; finite-mass hammer recoil; source-equation surface correction and self-contact; explicit measurements and diagnostic definitions; and atomic checkpoint/restart.

Hook presence is verified. The assumed 60 mm rise and other undimensioned positions remain hypotheses. A directly inspected damage field has a broad central failed region; agreement with Bowen's punching/shear crack morphology is not established.

## Unresolved original-model definitions

The inspected material does not establish Chapter 3 interface `bcc`, `s01`, `sc1`, control volumes, exact coupling quadrature, support bearing/constraint regions, complete reinforcement bends, rate/history policy, hammer initialization, measurement node or full-event termination. Chapter 2 and Appendix A also disagree in their rate-law exponent. The supplied Fortran uses another beam/loading case and cannot identify these benchmark inputs merely by being supplied alongside the thesis.

The chosen initial sphere position gives about a 0.2415 ms free-flight interval before cell-centre contact. The original time origin and initialization remain unknown. Results cannot be aligned by silently selecting a gap or replacing the displacement statistic with whichever number is closest to 22.70 mm.

A limited public GitHub repository search found no verified author implementation. General search/API requests were blocked by the runtime's domain policy; this is not proof that no public package exists. The large `Peridynamic articles.zip` remains uninspected because the upload download tool limits files to 32 MiB. This limits claims of complete coverage of the supplied archive collection.

Further runs of a declared hypothesis can characterize that hypothesis, but cannot establish the unknown original executable definition. Matching two rounded peaks by selecting unreported inputs would not demonstrate exact model reproduction.

## Files and restart

The current independent source is in `pdbeam/`, with `AUDIT.md` and `SOURCE_PROVENANCE.json`. `outputs/hooked_report/` contains the irreversible comparison; `outputs/printed_report/` contains the printed-history comparison. CSV histories, manifests, native checkpoints and diagnostic snapshots remain in the workspace. The compact delivery package includes the current source, tests, reports, complete comparison histories and the irreversible run's final native checkpoint. Run `build.py` to regenerate the reference mesh before plotting.

Use a new output directory for continuation; `run.py` refuses to overwrite an existing history. Its step target is absolute:

```bash
.venv/bin/python run.py --steps 10000 --restart outputs/hooked_dt3e7/checkpoint.npz --out outputs/new_continuation --accept-bulk-interface-hypothesis
```

The restarted calculation remains an explicitly labelled reconstruction. This command does not promise the original author's peaks or establish a full-event maximum.
