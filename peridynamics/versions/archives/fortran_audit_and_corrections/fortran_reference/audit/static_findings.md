# Independent audit of the supplied Fortran model

Audited source: `../Console2.original.f90`, with line numbers referring to that
file, and the existing completed run in `../runs/literal_100000`. The source was
not modified for this audit. No material parameters were fitted. A finite run
and agreement between two early elastic runs do not establish physical validity.

## Main finding: a force-component error explains the steel runaway

Lines **746–747** assign the second steel node's transverse components from its
**axial** component:

```fortran
pforce22(ccnode,2)=pforce22(ccnode,1)+0.0d0
pforce22(ccnode,3)=pforce22(ccnode,1)+0.0d0
```

These are in the steel tensile-yield branch. For an isolated yielded element,
its left node receives `(F,0,0)` and its right node receives `(-F,-F,-F)` rather
than `(-F,0,0)`. Physical action and reaction no longer balance. At a shared node
the transverse value depends on the current accumulated axial force, so the
effect also depends on element traversal order. This is a proved component-index
error, rather than an uncertain interface parameter or mesh sensitivity.

The magnitude of the source's tensile cap (lines 742–745) is

```text
F = 416 MPa × 530.66 mm² = 220754.56 N
V_s = dx × 530.66e-6 = 1.008254e-5 m³
q = F/V_s = 2.1894736842105263e10 N/m³
rho_s = 7850 kg/m³, cn = 3e6 kg/(m³ s)
```

The source's steel velocity update (lines 912–926) applies

`rho_s * dv/dt = -q - cn*v` to each incorrectly forced component. Its terminal
velocity is therefore

`v_inf = -q/cn = -7298.245614035088 m/s`,

and its damping time is only `rho_s/cn = 0.002616667 s`. The computed initial
spurious acceleration is **2.789138e6 m/s²**. Thus hundreds of metres of motion
are an expected result of the erroneous equation when active for approximately
0.1 s. Reducing the timestep alone does not remove this forcing.

The existing full-run state confirms the mechanism directly:

| Observation | Completed-run value |
|---|---:|
| Node 2 final y velocity | −7298.245614034692 m/s |
| Node 2 final z velocity | −7298.245614034692 m/s |
| Relative difference from predicted terminal magnitude | 5.42e−14 |
| Steel transverse components matching terminal magnitude | 720 of 724 |
| Maximum displacement node | Steel node 181 |
| Node 181 displacement vector | (−749.145, −751.627, −751.182) m |
| Maximum steel displacement magnitude | 1300.168099603453 m |
| Maximum concrete displacement magnitude | 1.4943543281281102 m |
| Maximum steel speed | 12640.932209624805 m/s |
| Maximum concrete speed | 2.315952100228636 m/s |
| Final steel kinetic energy using declared nodal mass | 1.5261075146e9 J |
| Final concrete kinetic energy | 0.5575029858 J |

The literal result is dominated by pathological steel forcing. This proves an
explanation for the approximately 1300 m displacement. It does **not** prove
that the entire contact-force peak comes from this error or that the remaining
model is valid after changing two lines. A corrected rerun is needed.

The minimal component correction replaces the right-hand sides with their
matching y and z accumulators. It must be labeled a corrected-source variant;
the original and its invalid completed run remain an immutable reference.

## Other proved mathematical and implementation issues

### 1. Steel–concrete exchange conserves density, not physical force

The interface uses `dforce = bc*s*vol*direction` at lines 779–786, which has units
N/m³. It adds that density to a steel node and subtracts weighted copies of the
same density from eight concrete nodes (lines 792–818). The interpolation weights
sum to one. However, steel and concrete representative volumes differ:

```text
V_c = dx³ = 6.859e-6 m³
V_s = 1.008254e-5 m³
V_s/V_c = 1.469972299168975
```

For one interface interaction, physical nodal forces sum to

`q*V_s - sum(N_a*q*V_c) = q*(V_s - V_c)`,

which is nonzero. The acceleration denominators are the respective densities
(lines 896–914); multiplying accelerations by their physical nodal masses does
not cancel this residual. For the declared source masses, a minimal conservative
transfer multiplies every concrete return density by `V_s/V_c`, or accumulates
physical force first and divides by each endpoint's volume. This is a units and
momentum correction, not a free parameter. A different definition of steel
representative mass would need explicit justification.

### 2. Steel yield criterion and force are inconsistent at finite deformation

The criterion uses current three-dimensional element length at lines 719–721,
but the elastic force uses a fixed x-direction stiffness matrix constructed at
lines 551–581 and applied at 722–740. The yielded force also remains in reference
x, regardless of current element orientation. There is no plastic strain or
return map. Only tensile geometric stretch is capped; negative strain is not.

For a rod undergoing a rigid 90-degree rotation, current length is unchanged
and the geometric stretch is zero, but the fixed stiffness produces a nonzero
axial force of magnitude `E_s*A ≈ 110.433 MN`. For zero axial displacement and
a transverse difference of just over **1.202 mm** between consecutive 19 mm
nodes, the nonlinear length criterion crosses 0.002; force switches from zero
reference-axial elastic force to approximately **220.755 kN** reference-axial
cap. This discontinuity is a proved mismatch between the criterion and force.

A small-rotation reference-linear rod may have been intended, but cannot be
trusted after large rotations or inversion. Replacing it with a current-direction
corotational rod is a declared model correction, not an exact literal copy of
the source. The component typo can be corrected separately to isolate its effect.

### 3. Damage flags do not gate force

Concrete force is calculated from current stretch at lines 654–670; `fail(i,j)`
is set to zero at 673–676 and used only in the damage sum at 679. Likewise,
interface `failfail` at 820–826 does not gate the force at 779–818. Consequently,
a previously failed bond can transmit force again if its stretch returns below
the failure threshold. A displayed persistent crack/damage pattern and the
active force network therefore mean different things.

This is proved source behavior. Whether it is an accidental omission or the
intended instantaneous envelope law is not established by this source alone.
An irreversible variant must be explicitly named and independently compared;
it must not be silently described as the original model.

### 4. Support displacement and velocity constraints disagree

At lines 931–952, the code sets constrained x or z displacement to zero after
the dynamic update but leaves corresponding velocity and acceleration unchanged.
Contact subsequently may overwrite constrained values. A physically fixed
degree of freedom should maintain both its displacement and velocity and record
the corresponding reaction. The current operation can leave fictitious kinetic
state and makes energy/momentum closure undefined without an external constraint
impulse ledger. This is a proved kinematic inconsistency; it has not been shown
to be the cause of the 1300 m motion.

### 5. Contact uses concrete mass for every node

Both cylindrical contact loops cover all nodes (977 and 1019), including steel,
but reactions use `dens*vol` (1007–1009 and 1049–1051). If a steel node contacts
an indenter, the recorded reaction uses approximately 0.0171475 kg instead of
its source mass 0.079147939 kg. It also resets y velocity/displacement despite
computing contact geometry only in x–z. Contact should either deliberately
exclude steel or use each particle's actual mass and preserve free axial
motion for a frictionless cylinder. Whether steel actually contacted an
indenter during the literal run requires event instrumentation; the branch
error is proved but its contribution to the reported peak is not yet measured.

The code divides by `dpend`/`pdpend` without a zero-distance guard (986–988,
1028–1030). No zero-distance failure is observed in the completed finite state.

### 6. Angle evaluation divides by zero for legitimate bonds

At line 646, `atan(abs(dy)/abs(dx))` receives `dx=0, dy!=0` for legitimate
nonvertical 3D bonds. IEEE infinity gives the limiting angle, but trapping
division by zero stops execution. `atan2(abs(dy),abs(dx))` removes this
portability error without changing the intended angle. The existing 200-step
comparison confirms only early elastic execution equivalence, not nonlinear
model validity. This does not explain the demonstrated steel terminal velocity.

## Measurement, loading, and diagnostics

- The source's top reaction sum at lines 689–693 is over an **empty range**:
  extra top nodes are commented out and `totint == tottop`. `loadxx`, `loadyy`,
  `loadzz` therefore stay zero; the nonzero reaction is `forceld` from contact.
- The original probe (699–703) lies at `(0.0095,0,-0.0095)` m, not the bottom
  midspan fibre. It is sampled before the current integration/contact update,
  while contact force is sampled afterward. A force–displacement graph mixes
  different update states by one step.
- `dispdispdisp` is evaluated at 956–958 **before contact projection**, then
  written alongside post-contact displacement, so its output norm can be stale.
  The high-precision final-state norm independently confirms the 1300 m result;
  it is not an artifact of this stale norm.
- Prescribed quintic loading (1063–1071) moves two 10 mm diameter cylinders
  downward 10 mm in 0.3 s. There is no freely moving finite-mass impactor or
  impact-momentum equation. The source is a different beam/loading model from
  Bowen Chapter 3. Fixing code defects cannot by itself change this fact.
- Physical densities are used by the active update. `massvec` (539–549) belongs
  to an unused, commented adaptive-dynamic-relaxation scheme (868–890); changing
  `massvec` has no effect on the active calculation.
- The active update is velocity-first semi-implicit Euler with explicit damping.
  Its damping-only factors are stable (`cn*dt/rho_c=0.0036`, steel≈0.0011465).
  This audit has not established complete stiffness/contact timestep stability.
  The analytical terminal-velocity evidence makes an arbitrary timestep change
  an inadequate explanation for the current failure.
- No energy, physical momentum, reaction/support impulse, plastic work, or
  constitutive-branch counters are maintained by the original program. Finite
  arithmetic is therefore not a physical acceptance criterion.

## Initialization: what is proved, and what remains unverified

The primary displacements, velocities, body-force arrays, force accumulators,
damage flags, family counts, and point labels are initialized explicitly. The
active old displacement/velocity arrays are assigned during every integration
update before contact reads them. `cn` is explicitly set to 3e6. The supplied
source's declared but unused rod matrices and relaxation scalars cannot be
blamed for the active result merely because they are uninitialized.

Interpolation connectivity `node`/`nodenode`, shape-function arrays, and the probe
array have no explicit initialization or validity sentinels. They rely on the
geometry searches finding every match (315–338, 418–442, 699–703). The present
regular geometry appears to cover the intended domain and the finite completed
run provides no positive evidence of a failed lookup. Instrumented assertions
should check eight valid corner indices, bounded/finite weights, partition of
unity, affine reproduction, nonzero reference lengths, and exactly one probe.
Undefined reads are a risk to verify, **not a proved cause** of this run.

## Recommended correction and verification sequence

1. Correct only the two yielded steel component assignments and add branch,
   maximum-motion, momentum/force residual, finite-value, and energy diagnostics.
   Compare to the literal reference before and at first plastic activation.
2. Correct physical interface action/reaction using the declared endpoint
   volumes; verify a one-interaction physical momentum balance and assembled
   zero-external-force balance. Keep this change separately identifiable.
3. Enforce support velocities consistently, use type-specific contact masses,
   and align force/displacement sampling. Verify reaction impulses and work.
4. Decide and document the finite-deformation rod formulation and failure-history
   semantics. These are larger scientific model choices than the index/units
   fixes. Test rigid motions, yield transitions and monotonic/unloading behavior.
5. Run short identical physical windows with timestep refinement, then complete
   the same source loading schedule once diagnostics pass. Do not tune material
   constants to recover an expected force.

This is a defect audit of the supplied reference, not evidence that this source
is Bowen's Chapter 3 input or that a corrected-source variant reproduces that
chapter's experimental/numerical results.
