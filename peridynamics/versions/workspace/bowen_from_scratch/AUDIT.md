# Source equations and reconstruction decisions

This is an independent implementation in `pdbeam`, with no imported notebook kernels or Stage 19 results. The primary sources are Bowen's Chapter 2–3 equations, Figure 3-1 and Appendices A–C. The supplied Fortran is useful implementation evidence but describes a different beam/load, so its constants are not transferred to this benchmark.

## Verified corrections relative to errors identified in the earlier notebook audit

1. **Units:** force-density micromodulus is N/m^6. A unique pair's physical force is c × stretch × Vi × Vj × corrections, in N. Divide by nodal mass, not density, when advancing that force. Steel nodal mass follows half of adjacent rod-element volumes. Interface reactions use the same physical force with opposite signs and shape weights; differently sized steel and concrete control volumes never cancel by accident.
2. **Interface geometry:** eight-node trilinear interpolation reproduces affine displacement fields. Three-dimensional spherical sample offsets replace the rank-deficient two-sample family identified in the old implementation. This establishes interpolation consistency, but does not establish the unknown benchmark interface constants or sampling quadrature.
3. **Bar height:** Figure 3-1 places the centre 80 mm below the top of a 100 mm beam: 20 mm above the bottom.
4. **Displacement:** the reported probe is the mean vertical displacement of bottom-layer particles at midspan, across the width. It is labelled that way; the thesis's original transducer/node definition remains unknown. Maximum particle motion is a separate diagnostic.
5. **Contact:** force is computed from the change in physical nodal momentum required by spherical projection, and the opposite axial impulse is applied to the finite-mass hammer. Horizontal impulse is accounted as the axial guide's external reaction. Contact energy change is reported explicitly; projection is not presumed energy conserving.

## Explicit unresolved choices in the exploratory run

| Item | Current reconstruction | Status |
|---|---|---|
| Rate law | Eqs. 2-13–15, signed bond extension rate; thresholds reset from static values each step | Appendix A uses a different empirical exponent; meaning of rate sign in original Chapter 3 is unresolved |
| Dynamic modulus | DIFE divides threshold enhancement; c remains static | Mirrors visible threshold-only snippets; alternative multiplying c is not assumed |
| Failure | Irreversible concrete and interface flags gate subsequent force | Supplied snippets do not visibly gate by previous flag; author-run equivalence unresolved |
| Steel | Corotational axial rod, elastic-perfectly-plastic return mapping, 482 MPa yield | Objective reconstruction; supplied different-case source uses reference-linear rod force and a cap, not this plastic history |
| Bar layout | Two 550 mm longitudinal legs, x25–575 mm, y20/50 mm, z20 mm; vertical end hooks rising 60 mm | Hook presence is verified from Fig. 3-2; rise, ends and lateral centres are reconstruction assumptions. 60 mm rise gives 20 mm top and bottom bar-centre cover |
| Interface | Concrete bulk weighted micromodulus and static concrete thresholds; lumped rod volume × dx³ sample volume | **Hypothesis**, not a recovered definition of bcc/s01/sc1/vol |
| Interface samples | Integer offsets inside a 15 mm sphere; eight-node interpolation; boundary samples omitted | Pattern from supplied Fortran; exact Chapter 3 quadrature unknown; zero stiffness at outer radius |
| Supports | Bilateral vertical constraints on 10 mm bottom patches at x50/550 mm; minimal in-plane anchors | Width, uplift, bearing contact and anchors unresolved |
| Hammer | Sphere R60 mm, axial guide; bottom initially at geometric beam top, no gravity | Initial contact with cell-centre particles is delayed by the half-cell surface offset; original placement unresolved |
| Surface correction | Three uniaxial energy calibrations using source 0.6 E ε², pair ellipsoidal directional blend | Matches method seen in supplied source, not proven original benchmark implementation; not exact pointwise continuum matching |
| Self-contact | Eq. 2-48–49, current-coordinate neighbour search, constant c=12E/(πδ⁴) | Assignment of weighted c for new neighbours is unresolved |
| Integrator | Velocity then displacement, Eqs. 2-42–43 | Contact prediction then hammer velocity recoil; exact original coupling order not available |

No parameter is selected by minimizing differences from 359.8 kN or 22.70 mm. The CLI requires `--accept-bulk-interface-hypothesis` so an exploratory run cannot silently be mistaken for the fully specified author model.

## Diagnostics and numerical verification

Elementary checks cover rigid-motion invariance, affine interpolation, physical force units, action/reaction and interface torque, rod yielding history, contact impulse, analytical force versus energy gradients, rate-potential exchange and rupture energy deletion. A full assembled contact state is continued both directly and through a checkpoint; states are compared exactly. Checkpoints contain displacement, velocity, previous displacement, concrete/interface histories, rod plastic strain, hammer state, threshold arrays and cumulative diagnostics. Manifest equality includes parameters, choices and all package source hashes.

Stored energy is sampled every 100 steps at the end state, with the force step's thresholds frozen. The ledger includes hammer/concrete/steel kinetic energy, bilinear bond/interface potentials, rod elastic energy, short-range potential, deleted rupture energy, plastic dissipation, rate-induced changes in potential, contact kinetic-energy change and constraint kinetic-energy change. The remaining residual is reported, not forced to zero. Explicit integration, projection changes to potential, constraint treatment and plastic return can contribute numerical residual. A separate internal-path-work residual is also reported. Unsampled stored-energy/residual CSV fields are NaN, not zero.

Momentum residual compares all nodal and hammer momentum against initial momentum plus accumulated support and guide impulses. Damage flags count failed bonds; this is not automatically proof of a correctly reproduced punching/shear crack. Damage images must be assessed separately from numerical peaks.

The 5,000-step run spans 1.5 ms, sufficient to compare a damage interval shown in the thesis but not to establish the full-event maximum deflection. A half-timestep run to 0.6 ms tests early impact discretization sensitivity. Mesh convergence and complete-event stopping remain additional work; exact original results cannot be inferred from a short run or two rounded targets.

## Projection energy follow-up

The base residual subtracts contact/constraint **kinetic** changes. Projection also changes configuration and potential, so this residual alone is not a closed-system physical-energy error. A separate checkpoint-continuation audit computes frozen-parameter potential changes over all affected concrete/interface pairs and the union of short-range neighbours before/after projection. It does not change force, damage or motion; an independent comparison verifies identical mechanical states with and without the audit. The adjusted segment residual still contains the pre-segment offset and must be compared to its start value.

`analyze.py` additionally reports the closed-system material energy drift: kinetic + stored + deleted rupture + plastic dissipation − initial kinetic − rate-potential exchange. Hammer contact is internal to this system; its kinetic substep change is not an external energy supply. This diagnostic includes numerical effects of contact, supports and integration rather than subtracting one piece of contact energy and implying conservation. Neither diagnostic is forcibly normalized.

`analyze_connectivity.py` identifies components attached to vertical support particles through surviving nonzero-stiffness concrete bonds, surviving positive-weight interface links and rods. This distinguishes detached flight from supported material. It is a connectivity diagnostic, not the author's displacement-sensor definition.

The assembled surface correction is checked by applying independent 10⁻⁶ uniaxial strain fields and summing actual corrected pair energies. `outputs/surface_patch.json` reports global ratios and pointwise extrema against the 3D ν=0.25 continuum density. This check evaluates the assembled operator and exposes its approximation rather than asserting exactness from the calibration step.

## Verified topology and checkpoint corrections

The reinforcement figure on printed page 31 (Fig. 3-2, PDF page 48) visibly contains vertical end hooks. The current builder therefore includes both hooks on both bars. Their presence is source-supported; their rise is not dimensioned, so the exploratory configuration explicitly assumes 60 mm rise (20 mm bar-centre cover at top and bottom). Piecewise straight rod segments represent the undimensioned bends. Tests verify connectivity, total rod length and lumped mass. Straight-bar outputs are retained as a diagnosed earlier reconstruction, not as an accepted author-equivalent result.

A checkpoint reader encountered a partially written file during the initial continuation. Current writes use a temporary file in the same directory, flush/fsync, then atomic replace. Tests inject a write failure and verify that the previous checkpoint survives unchanged and temporary files are removed. Original completed checkpoints retain their pinned source versions. The explicit atomic-I/O-only migration tool rejects physical/configuration changes; it is not a way to restart the straight-bar state as a hooked model.

Failed irreversible bonds are skipped before rate/force evaluation to reduce cost. Their force and potential remain zero; failure flags remain false. This leaves failed-bond threshold metadata frozen, which cannot affect a run with fixed irreversible-history policy. A source/configuration mismatch still rejects restart under a different history policy.

## Printed history comparison

`--instantaneous-history` evaluates the scalar force law regardless of earlier failure flags, consistent with the visible Appendix C / supplied Fortran force evaluation. Failure flags still record whether a bond ever exceeded the threshold, as in the printed bookkeeping; under this option they do not imply permanent force removal. The concrete use of this interpretation in the original Chapter 3 executable remains unverified. It changes no published material parameter. The comparison distinguishes reproducing printed force behavior from imposing irreversible fracture. The capped primitive is used as the instantaneous potential and rupture deletion is zero under this interpretation; this is not the same dissipative fracture model as the irreversible branch.

`verify_stability.py` gives a sufficient bound for the initial elastic concrete/interface/rod assembly through Rayleigh and Jensen inequalities. It excludes evolving self-contact and nonlinear geometry. At the current mesh/parameters the bound is about 0.868 µs, while the used timestep is 0.300 µs. This is initial linear stability evidence, not force-peak convergence or experimental validation.

Additional interface discretization choices are explicit: the exploratory interface coefficient includes the same radial partial-volume factor as the concrete operator, but Appendix C does not visibly multiply by this factor. Whether its omitted `vol` initialization already embeds quadrature corrections is unknown. Interface surface correction is not applied. Neither decision is verified as the Chapter 3 initialization. The concrete mesh occupies the full beam volume while rod mass is added through the hybrid mass matrix; concrete mass displaced by reinforcement is not subtracted, consistent with the described embedded PD/FEM assembly, but the original mass initialization is still unavailable.

## Damage comparison metric

The primary report now follows Eq. 2-10: failed reference-family neighbour volumes divided by total reference-family volumes. With uniform dx³ this is the failed-neighbour fraction, including outer-radius zero-micromodulus bonds present in the family. Global damage is the mean over all concrete particles. Earlier force-carrying, partial-volume-weighted images are a separate diagnostic and are not the thesis's damage definition. Connectivity continues to exclude zero-stiffness bonds because they cannot mechanically attach a particle to a support. Report images show damage in reference coordinates; crack openings and deformed geometry require separate interpretation.

Connectivity uses current stretch and the last force step's terminal threshold. Previous failure flags gate links only in the irreversible branch. The instantaneous branch can restore a link despite a recorded prior failure, so removing every ever-failed link would misclassify most of that branch's still-connected material as detached. Damage and mechanical connectivity are distinct recorded quantities.
