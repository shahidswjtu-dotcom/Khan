# Force units, rod assembly, constraints, and stability audit

This audit applies to the supplied `Console2.original.f90`, with its actual 19 mm spacing, 3 μs step, 26 mm longitudinal rods, 2500/7850 kg/m³ concrete/steel densities, and prescribed two-cylinder loading. It does not substitute Chapter 3 material or mesh values.

## Proven cause of the metre-scale runaway

The yielded rod branch at original lines 746–747 copies the right node's **axial** force-density accumulator into its **y and z** accumulators. In the simplest zero-transverse-deformation case, the compiled original source produces physical nodal forces

- left node: `(220754.56, 0, 0)` N;
- right node: `(-220754.56, -220754.56, -220754.56)` N.

Thus one internal rod has a net artificial transverse force of −220.75456 kN in each transverse direction. Shared-endpoint accumulation also suffers assignment rather than component-preserving accumulation. The verified correction preserves previously accumulated y/z contributions; it does not indiscriminately zero them.

The active source uses damping force density `cn*v`, with `cn=3e6`, and steel nodal volume `vol1=0.019*530.66e-6=1.008254e-5 m³`. A persistently wrong yielded-endpoint transverse force therefore approaches

`v_terminal = -220754.56 / (3e6*1.008254e-5) = -7298.245614035088 m/s`.

The completed run's node 2 has `vy=vz=-7298.245614034692 m/s`, agreeing within `4e-10 m/s`. **360 steel nodes have both transverse velocities at this wrong terminal value.** Node 2's final transverse displacements are approximately −747 m; the largest displacement norm is 1300.1681 m. This matches the deterministic faulty rod force and damping equilibrium. The extreme motion is not an unexplained floating-point failure.

## PDROD physical force and work transfer

Original force-density variables have units N/m³ after neighbor-volume integration. The concrete sample volume is `vol=dx³=6.859e-6 m³`; the declared steel rod volume is `vol1=1.008254e-5 m³`, a ratio of **1.4699722992**.

For one interpolated interface link with steel density contribution `d` and partition-of-unity interpolation weights `N_a`, the original adds `d` to steel and `-N_a*d` to concrete. This is adjoint in unweighted density space. If nodal masses/forces use the declared material volumes, its physical internal resultant is

`d*vol1 - sum(N_a*d)*vol = d*(vol1-vol)`,

which is nonzero. Preserving the original steel-side constitutive force, each concrete-side density must instead be `-N_a*d*(vol1/vol)`. This restores both physical force balance and virtual work. The compiled native transfer test reduces its example resultant from `(3223.54, -6447.08, 9670.62)` N to numerical roundoff.

The normalization statement explicitly adopts the steel rod mass `rho_steel*vol1`, consistent with the source's axial rod force divided by `vol1`. Interpreting all nodes as having the common `vol` would conserve this interface force, but would scale the physical rod stiffness/cap by `vol/vol1` and give a rod mass inconsistent with its stated area. The source mixes these two conventions rather than declaring a consistent choice.

The diameter-derived `A1=pi*(26 mm)²/4=530.929158 mm²` and the volume/cap area `530.66 mm²` differ by about 0.051%. This small original input inconsistency is documented; changing it is unnecessary to explain the runaway.

## Constraint and contact checks

The original four kinematic constraints reset displacement only (lines 933, 938, 945, 950). Their constrained velocities remain nonzero after force integration. This is a velocity-constraint inconsistency; it does **not** imply that an exactly reset support position drifts. A compiled extracted-block check verifies that setting only the corresponding velocity components to zero preserves all other components and original displacement constraints.

Both original contact loops use `dens*vol` for reaction momentum, even though they loop over steel as well as concrete nodes. The verified mass correction uses `dengs*vol1` for a contacted steel node; the concrete reaction is unchanged. For the native example, the original steel reaction is only 21.66513% of the mass-consistent value. This matters if steel reaches the cylindrical contact region, as it can after the runaway.

The x–z cylinder's normal has no y component, yet the original projection sets both `disp(i,2)` and `vel(i,2)` to zero and records no y reaction. For a frictionless cylinder, these tangential components must remain unchanged. The compiled whole-contact-loop check verifies the selected correction preserves y displacement/velocity, retains the identical x–z projected point and velocity, and computes the reaction from the appropriate nodal momentum change. Interpreting the original y reset as an intentional lateral clamp would require documenting and recording that additional constraint reaction; the supplied source does not do so.

## Stability: what the checks establish

The active integration is semi-implicit Euler with explicit damping, not the commented adaptive dynamic-relaxation update. `massvec` is initialized but is not used by the active stepping equations.

For a scalar initial elastic mode with frequency squared `lambda` and damping `gamma=cn/rho`, the source's update is

`v_new=(1-gamma*dt)*v-lambda*dt*u`, `u_new=u+dt*v_new`.

Its stability condition is `lambda*dt² + 2*gamma*dt < 4` (for a nonzero elastic mode). Here `gamma_c*dt=0.0036` and `gamma_s*dt=0.00114650`, so the force-free damping recurrence is stable.

For an isolated uniform axial rod chain with the source's `EA1/dx` stiffness and `rho_s*vol1` masses, the high-frequency limiting eigenvalue is `4*EA1/(rho_s*vol1*dx)=2.9374243608e11 s⁻²`. Its critical step is **3.68887155 μs**, above the source's 3 μs. This is only a rod-chain check; it does **not** establish the critical step of the assembled concrete/interface/contact system. A full tangent/evolving-contact study is still required before declaring the complete simulation converged. Reducing the time step cannot fix the spurious nonzero internal rod resultant.

The reference-linear truss matrix also produces approximately −551.7 kN under an exact 0.1 rad rigid rotation with unchanged rod length. The source mixes a finite-length stretch to choose its branch with a fixed reference-axis small-displacement force matrix. This confirms loss of finite-rotation objectivity; it may exceed the author's intended small-rotation use rather than being an assignment typo. Changing the rod kinematics would be a separately documented model correction.

## Verification artifacts

- `units_stability_check.py` and its JSON: constants, analytical checks, and direct completed-run velocity corroboration.
- `test_native_operators.py`: seven tests compile exact operator blocks extracted from original or explicitly corrected complete Fortran source. They cover yielded rod resultant, accumulator preservation, unchanged preyield/compression, interface force/work transfer, reaction mass, velocity constraints, and whole-cylinder tangential projection.
- `native_operator_verification.json`: numerical outcomes plus original-byte and normalized-source hashes. Generated harnesses are under `native_operator_build/`.

All seven native tests passed. These are implementation checks, not a validation of the supplied benchmark against Bowen's Chapter 3 results. No material parameter was fitted or changed by this audit.
