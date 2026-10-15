# Equation-to-code audit

Read the supplied handoff completely before edits. Document instructions are research context; the user's request authorizes the work. Original uploads and extracted source notebooks remain unchanged. Baseline copies and SHA-256 checksums are in `baselines/stage19`.

Page references below give **printed thesis page / 1-based PDF page** (offset 17). Thesis: Tu Bowen, 2025, 基于近场动力学方法的混凝土梁构件冲击损伤研究. The primary experimental reference is Fu et al., IJIE 137 (2020) 103432, not Wang's 2022 numerical study.

## Inventory and reading

- Handoff Markdown: read in full, 30,444 bytes.
- Thesis: Chapter 2 theory and implementation, Chapter 3 baseline and damage chronology, Appendices A–C extracted and studied. Figure 3-1 independently rendered and inspected.
- All seven supplied Stage 13–19 notebooks statically parsed (zero syntax errors); definition/cell inventory in `sources/notebook_inventory.json`. Only the selected active baseline definitions executed. Old experimental cells are not executed on import.
- `Console2.f90`: declarations, material initialization, reinforcement, sample topology, interpolation, force assembly, acceleration and contact reviewed. Source from `5(成功).zip` extracted for comparison. `1.zip` contains a different `Console322.f90` (background, not baseline).
- Wang et al. 2022 supplied PDF: Sections 2.2–2.3 and 3.2 studied, including contact equations and RC beam validation. Confirms projection contact and benchmark inputs; does not supply Bowen's 3D PDROD quadrature.
- Review PDF downloaded; background only, not used to define benchmark constants.
- `Peridynamic articles.zip` exceeds the tool's **32 MiB transfer limit** and could not be inspected. This does not prevent inspection of separately supplied primary thesis/source/Wang paper. The executed Stage 19 notebook mentioned in the handoff is not present in the bundle; its saved CSV/JSON baseline is present.

## Equation / statement mapping

| Source | Physical quantity and units | Baseline code | Audit status / action |
|---|---|---|---|
| Eq 2-1; 2-38, p10/27 and p24/41 | acceleration from force density, N/m³ divided by kg/m³ | Stage19 cell94 `_fast_concrete_kernel`, `model.evaluate` | Pair kernel produces N by integrating over **both** point volumes; divide by nodal kg. Correct convention for Eq 2-30. |
| Eq 2-5–2-8, p11–12/28–29 | bilinear tension law; s0=ft/E; sc=5Gf/(6Eδs0) | cell94 | Persistent full rupture preserved; partial-softening unloading is an instantaneous envelope, not a separately specified damage law. Compression remains linear and unbounded. |
| Eq 2-12, p13/30 | c=72E/(πδ⁴)(1-r²/δ²)², **N/m⁶** | cell94 | Distance kernel agrees. Stage19 additionally multiplies c by DIF_E; not visible in Appendix A/B. Mark as reconstruction choice. |
| Eq 2-13–15, p14–15/31–32 | strain rate s⁻¹, DIF, invariant-base thresholds | cell94 | Uses theory exponent 1/(10+6fc/10MPa). Signed positive rate convention. No compounding. Appendix A uses a different empirical exponent and gives no resetting context. |
| Eq 2-26–27, p19/36 | concrete lumped mass ρdx³; rod half-element masses, kg | cells56/117; `geometry` | Mass total and endpoint half masses tested. No steel-density-to-nodal-mass substitution. Concrete volume is not reduced for steel occupation; retain and disclose. |
| Eq 2-30–31, p19/36 | pair stiffness cViVj/r and rod EA/l, N/m | cells94/44 | Concrete finite geometry; rods use reference-direction small-strain axial force. No rod yield/failure despite supplied fy. |
| Eq 2-33–37, p20–21/37–38 | QᵀTᵀKihTQ, Kih=cViVh/r | cells42/118 | Eight-node extension supported by Appendix C. Partition of unity, linear field and deformed moment conservation verified. **Two lateral samples are not the thesis's full family.** |
| Appendix C p89–91/106–108 | interpolated displacement, s, bc, dforce=bc·law·vol·unit | cells104/118, `interface` | Literal kernel retained separately. dforce is density if vol has m³. `source_density` cancellation in array values is not physical momentum balance at unequal volumes. |
| Fortran lines 93–97, 126–127, 894–922 | vol=dx³; dens=2500,dengs=7850; acc=pforce/density | `Console2.f90` | Different configuration, dx=19 mm, dt=3 μs, 26 mm main bars, damping cn=3e6. Not transferable Chapter3 inputs. |
| Fortran lines 240–278 and 434–455 | 7³ offset lattice trimmed to spherical horizon, trilinear stencil | cell42 | **Recoverable algorithmic convention:** many nonzero 3D offsets, not just ±y. Benchmark truncation/volumes still require declaration. |
| Appendix C force followed by failfail update, p90–91/107–108 | rupture history | cell104 literal vs corrected; cell118 irreversible | Printed force ignores prior flag: reloading transmits force. Test demonstrates it. Baseline already corrects this; cannot label corrected history as literal. |
| Eq 2-39 p25/42 | horizon partial-volume factor, dimensionless | cell31 | Agrees on truncated r≤δ family. |
| Eq 2-40–41 p25–26/42–43 | directional energy surface correction | cell31 disabled | Verified missing method; occupancy ratio is not this correction. Do not substitute occupancy calibration silently. |
| Eq 2-42–43 p27/44 | explicit velocity then displacement | cell120 Verlet | Integration differs from printed method. Projection alternative uses semi-implicit explicit velocity/displacement; rate-law and splitting conventions remain disclosed. |
| Eq 2-45–47 p27–28/44–45 | relocation, corrected v, hammer reaction=-mΔv/dt in N | cell54 penalty | **Verified mismatch corrected separately:** `impact_contact.project_contact`, no tuned stiffness. Axial hammer guide and explicit recoil timing are declared assumptions. |
| Eq 2-48–49 p28/45 | short-range repulsion, csh=1.5c, ds=min(.9r0,1.35dx) | absent | Verified missing concrete self-contact. Critical when particles detach/overlap. Requires careful new-neighbour search, not just original surviving bonds. |
| Fig3-1 and §3.1.1 p30–31/47–48 | 600×70×100 mm; span500; 2Ø6; bar depth80 from top | cell7 | **Verified diagram mismatch:** baseline z=10 mm; reconstruction uses z=20 mm. End hooks are drawn but exact lengths and lateral positions not dimensioned. |
| Fig3-4; Table3-1 p32/49 | 359.8 kN, 22.70 mm vs 301.0 kN,23.38 mm | cell120 central upper-half mean | Probe is not an identified experimental transducer. New probes labelled by reference coordinates; no claim bottom centre layer equals original sensor. Fig3-4 is **force–deflection**, not force–time; x-axis 0–5 mm differs from tabulated max22.70, and caption cites [79] while table cites [78]. Do not fabricate time histories. |
| Fig3-8 p35–36/52–53 | damage snapshots .15,.30,.60,1.20 ms; plateau by4000 steps | Stage19 550 steps=.165ms | Baseline is truncated, far earlier than damage completion. No reproduction claim. |

## Ranked findings

**Proven implementation/method discrepancies**

1. Earlier density contribution divided by steel nodal mass has wrong acceleration units. Stage18 converts via steel control volume and balances physical N; this fixes units but not unknown sample quadrature.
2. Penalty contact is not Eq2-45–47; separately implemented projection has no arbitrary stiffness tuning.
3. Bar centre-depth differs from Fig3-1. Corrected alternative changes 10→20 mm without changing original baseline.
4. Two y-axis bonds provide zero first-order z and x coupling stiffness: local K∝Σk n⊗n has rank one. Fortran constructs a 3D horizon family; recoverable topology convention, exact benchmark family unspecified.
5. Surface correction, concrete self-contact, rod yielding/failure are absent. These are not remedied by adjusting material E or interface stiffness.
6. Legacy half-depth probe and sparse force sampling do not measure all-step reaction and a named local displacement.
7. Original solver cannot restart. New checkpoints include cached forces and committed stretch; re-evaluation on restart would produce erroneous zero rates.

**Missing source definitions**

Benchmark assignments of `bcc,s01,sc1,vol`, steel/virtual control volumes, exact 3D sample placement and weights, bar hooks/end positions/lateral cover, bearing widths and lift-off/anchor constraints, initial gap, hammer guide/rotation, concrete self-contact search, softening unloading, rate convention and exact transducer. Appendix C references `bcc,s01,sc1` but never initializes them. The provided Fortran uses `bc,s0,sc` and different values/geometry; do not invent identity.

**Source inconsistencies / unresolved hypotheses**

Appendix A low-rate exponent =1.016/(1+0.8×60) vs theory1/(10+6×60/10); high-rate theta also differs. Absolute rate is printed in Appendix C but unused in its law; it does not settle Appendix A rate sign. Multiplication `s0=s0*dif` is only a bug if base thresholds are not reset in omitted context. Large detached-node motion may be a compression/contact/model problem rather than integration blow-up; timestep checks and force/energy diagnostics must distinguish them.

## Correction policy

Keep legacy results immutable. All interface demo coefficients remain labelled, no parameter fitting. Conservative force convention: F_s = c·law·V_h·V_s·n (N), F_ca=-N_aF_s, masses in kg. This agrees with the volume-product structure of Eq2-34, but identifying V_h and V_s with dx³ and lumped rod volume remains a reconstruction, not recovered author initialization. `literal` returns source density entries and is only an elementary equivalence test; it is never passed to a physical-mass integrator as N.

## Implemented lattice alternative

`family_projection` implements the supplied Fortran's nonzero integer-offset spherical neighbourhood, giving122 samples per node (27,572 total) on the present reconstructed straight bars. Exact-horizon samples are allowed with the kernel's zero coefficient; the original Python input validator's strict `<delta` excluded these despite the source's `<=delta`. The validator now accepts the closed horizon within roundoff tolerance. Bulk-derived demo coefficients and V_h=dx³ remain unchanged and labelled. Full-rank elementary stiffness is verified; author's Chapter3 calibration remains missing.

## Follow-up: target original-results equivalence

The user requested exact Bowen results and confirmed no additional original source is available. Requirements/provenance are recorded in `EXACT_REPRODUCTION_REQUIREMENTS.md`. Further inspection found that `Console2.f90` concrete force evaluation also omits a visible prior-failure gate (lines654–668, failure bookkeeping672–675). The original Chapter3 behavior cannot be inferred from imposing irreversible damage in a physically corrected port.

`source_methods.py` now implements the independently specified missing energy-based surface correction (Eq2-40 with supplied Fortran directional blending) and current-neighbour short-range repulsion (Eq2-48–49). Four isolated verification tests pass. Surface factors on the5mm beam range0.969–3.177 per direction; disabling correction is therefore a material method deviation. Static constant-c self-contact is one explicitly labelled convention; a weighted c for newly contacting pairs outside the reference horizon is not uniquely specified. These methods are not yet integrated into the existing simulation; previous outputs and checkpoint-compatible core code remain unchanged.

Evaluation of repulsion on the saved legacy/family endpoints gives potential4883/8063J, respectively, under the declared convention, exceeding initial hammer energy2239.488J. This does not measure energy actually dissipated or created by a solver (the force was not applied); it demonstrates that the saved penetrated states cannot be reused as self-contact-consistent initial states. Any self-contact reconstruction must restart from undeformed geometry and undergo timestep/energy verification.
