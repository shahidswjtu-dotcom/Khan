# Trace of the completed supplied-Fortran calculation

This is read-only analysis of `runs/literal_100000`. The original source,
instrumented solver, histories, and nodal dumps have not been altered.
`trace_failure.py` can regenerate all results below. The final dump supplies
16-digit displacement and velocity data; earlier snapshots retain the source's
five-digit printed precision. Vector norms are recomputed from their components,
not taken from the scalar norm column, which was calculated before contact.

The largest movements are caused by a sustained spurious force in the yielded
steel branch. The terminal velocities provide a direct quantitative fingerprint
of that error; this is more specific than describing the outcome simply as an
unstable integration. This analysis does not establish what first triggered a
rod's yield branch because the full run did not record that event every step.

## Observed transition

| Step | Time (s) | Steel maximum displacement (m) | Concrete maximum displacement (m) | Rods currently above `sy` |
|---:|---:|---:|---:|---:|
| 60,000 | 0.180 | 0.0069828 | 0.0092512 | 0 |
| 64,000 | 0.192 | 0.0076368 | 0.0104181 | 0 |
| 65,000 | 0.195 | 1.1521580 | 0.5359694 | 229 |
| 66,000 | 0.198 | 20.5330993 | 1.4531427 | 360 |
| 70,000 | 0.210 | 162.5827314 | 1.4944034 | 360 |
| 100,000 | 0.300 | 1300.1680996 | 1.4943543 | 360 |

All rows except the final row use low-precision snapshots. The last row uses
the high-precision final state. Zero above-threshold rods in a snapshot does
not establish that the branch was never executed on an earlier step.

The vertical contact reaction reaches **1,059,806.44 N at step 64,881**
(0.194643 s), inside the interval where the model abruptly becomes unphysical.
The maximum probe reading, 10.602894 mm at step 99,720, describes one concrete
node; it does not reveal the catastrophic motion elsewhere.

## Exact force-to-velocity fingerprint

In `Console2.original.f90`, lines 746–747 assign the second endpoint's
transverse force components from its **axial** force component:

```fortran
pforce22(ccnode,2)=pforce22(ccnode,1)+0.0d0
pforce22(ccnode,3)=pforce22(ccnode,1)+0.0d0
```

For the all-yielded final state, replaying the rod loop in its original order
gives the following artificial force density:

- The first node of each bar has positive axial force and zero transverse force.
- Each interior node has zero final axial force but negative transverse forces.
  The next element cancels its axial accumulator after the previous element has
  copied the accumulator into the transverse components.
- The last node of each bar has negative force in all three components.

The source uses

\[
F_y=416\times530.66=220754.56\ \mathrm{N},\qquad
V_s=0.019\times530.66\times10^{-6}=1.008254\times10^{-5}\ \mathrm{m^3}.
\]

Therefore the magnitude of the imposed force density is
`21,894,736,842.105263 N/m³`. With damping `cn = 3,000,000`, its terminal
velocity is

\[
|v_\infty|=\frac{F_y/V_s}{cn}=7298.245614035088\ \mathrm{m/s}.
\]

Every component of the final velocity of all **362 steel nodes** agrees with
the exact replay prediction to within **3.96 × 10⁻¹⁰ m/s**. The artificial
net rod force is zero axially and **−79.4716416 MN in each transverse direction**.
The resulting final steel kinetic energy, using the mass actually assigned by
the source, is **1.5261 GJ**. Opposite internal rod forces should sum to zero;
these transverse forces fail that basic requirement.

The largest-displacement steel node is node 181. Its displacement is
`[-749.1451561, -751.6274571, -751.1822602] m`; all three velocity components are
`−7298.2456140 m/s`. Thus the late growth follows persistent forced drift at
the source's damped terminal velocity. Runtime checks do not reject this result
because these values remain finite and the indexing stays within allocated arrays.

## Concrete is also corrupted

The 1,300 m maximum is steel, but the failure is not confined to steel:
**5,028 concrete nodes exceed 100 mm displacement**, and the largest concrete
displacement is **1.4943543 m at node 5171**. Its displacement components are
`[-0.1841301, -1.0339645, -1.0630655] m`. Its final velocity is practically zero
because of damping and its fully failed damage state; a small final velocity
does not restore the physical validity of its displaced position.

The existing dumps prove the sustained-drift mechanism and its timing. A
controlled run with only the two component indices corrected is needed to
separate this error from other model deficiencies and to find the first trigger.
The 1,059.81 kN spike and final probe value cannot currently serve as validated
predictions or as evidence of matching Bowen's Chapter 3 model.

Artifacts: `trace_summary.json`, `trace_snapshots.csv`, and `trace_failure.png`.
