# Supplied Fortran model

`Console2.original.f90` is a byte-for-byte copy of the supplied Console2 source.
This model is independent of the previous Python reconstruction and Stage 19.
Executing the original Fortran is the reference for literal model equivalence;
it avoids introducing new constitutive or interface assumptions through a port.

The supplied program is **not Bowen's Chapter 3 impact input**. It has 36,036
concrete nodes, 362 steel nodes, two straight 26 mm rods, and a roughly
3.458 × 0.190 × 0.342 m concrete grid domain. The concrete coordinates span
3.439 × 0.190 × 0.323 m. Domain thickness/width beyond the outer node centres
is not explicitly declared, so those lengths describe the grid construction.
It uses two 10 mm diameter cylindrical indenters with a prescribed downward
10 mm displacement over 100,000 steps (0.3 s), rather than a freely falling
38.4 kg hemispherical hammer. These calculations cannot validate Chapter 3's
359.8 kN or 22.70 mm results.

## Run

On this Debian x86_64 workspace with system GCC 14:

```sh
cd /workspace/Khan/fortran_reference
python bootstrap_compiler.py
python run_reference.py --steps 200 --out runs/new_literal_200
```

With a compiler installed elsewhere, pass `--compiler /path/to/gfortran`.
Use a new output directory for each run; existing output is never overwritten.
For the original complete loading schedule use `--steps 100000`. Short runs
stop the time loop early while retaining the original `nt=100000` and the
literal `tt/100000` loading function. Short runs do not accelerate loading.
The complete 100,000-step calculation has not been completed here.

`run_reference.py` adds a high precision CSV and final-state dump, and records
source hashes, compiler and flags. It leaves the original mechanical operations
unchanged in literal mode. It does not add checkpoint/restart to the Fortran.
CSV displacement is the original single-node, pre-update observation at
(0.0095, 0, −0.0095) m; CSV contact force is the original post-projection
reaction. Their sampling order differs and is deliberately preserved.

## Runtime portability correction

The untouched source divides by zero at line 646 for bonds with Δx=0,
Δy≠0 and Δz≠0 when evaluating `atan(abs(Δy)/abs(Δx))`. With IEEE infinities,
the limiting angle is π/2; trapping division by zero stops the program.
This is an angle-evaluation portability issue, not evidence explaining the
Chapter 3 numerical discrepancy.

```sh
python run_reference.py --steps 200 --portable-angles --strict-zero --out runs/new_portable_200
```

This explicitly recorded variant substitutes `atan2(abs(Δy),abs(Δx))`
at both source occurrences. It changes no material input or loading input.
The literal variant traps invalid/overflow and checks bounds but permits IEEE
divide-by-zero. The portable variant can additionally trap divide-by-zero.

Both 200-step runs completed with finite history and state values. Their z
contact reaction and original probe histories agree exactly at written precision;
maximum nodal displacement difference is 6.21e−25 m. The maximum z reaction
was only 2.526 N during this early loading window, so this verifies execution
and the angle correction in the early elastic response, not fracture/plasticity
or the peak load. Full reference results are in `VERIFICATION.json`.

```sh
../.venv/bin/python verify_reference.py runs/literal_200 runs/portable_200 --report comparison.json
```

## Source-specific behavior retained

| Source lines | Literal behavior |
| --- | --- |
| 80–147, 201–204 | dx=19 mm, δ=3.015dx, dt=3 µs, Ec=30.5 GPa, Es=208 GPa, densities 2500/7850, damping coefficient 3e6 |
| 81–83 | Fixed thresholds s0=5.00456e−5, sc=1.25114e−3, steel yield stretch 0.002 |
| 206–231 | Two straight steel paths at y=±50 mm, z=−110 mm; 181 nodes each |
| 255–289, 414–438 | Virtual steel neighbours and eight-node trilinear concrete interpolation |
| 617–680 | Constant concrete micromodulus, partial volume and directional surface correction; failed flags affect damage but do not gate force |
| 709–750 | Reference-direction rod stiffness; tensile cap 416 MPa rather than a plastic return map |
| 744–745 | Yield branch assigns the second rod node's y/z forces from its x-force; apparent component-index error retained literally |
| 763–825 | Interface uses concrete bc, concrete vol and concrete s0/sc directly, without surface/partial-volume correction; equal opposite **force-density** interpolation |
| 686–697 | Top reaction sum is over an empty node range because extra top layers are commented out; these load columns are consequently zero |
| 916–949 | Supports reset selected displacement components, without matching velocity reset |
| 971–1050 | Sequential cylindrical projection, z/x reaction based on concrete density even when iterating steel nodes |
| 1061–1070 | Prescribed quintic displacement applied to the next step; no finite-mass hammer equation |

The interface's density transfer uses equal coefficients for steel and concrete,
although the declared representative steel and concrete volumes differ. A port
that instead conserves physical nodal force using different volumes changes this
source behavior. Likewise, correcting the yield component assignments changes
the literal solver. Such changes must be named corrected variants and compared
against this baseline, rather than presented as exact copies.

Compiler-specific initialization and floating-point behavior can affect results.
The original has no author compiler settings or reference output supplied.
Literal source execution is established here; equality to an author's historical
binary is not established.
