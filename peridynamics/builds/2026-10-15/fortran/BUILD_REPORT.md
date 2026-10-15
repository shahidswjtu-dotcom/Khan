# Fortran build verification — 2026-10-15

All six distinct recovered `.f90` files compiled successfully, unchanged, using
GNU Fortran (Debian 14.2.0-19) 14.2.0. All compiler logs are empty. Seven archived
native operator regression tests were rerun successfully. This is fresh build
and operator verification; no complete beam simulation was performed.

| Variant | SHA-256 prefix | Result |
| --- | --- | --- |
| Original two-dimensional Console322 | `6f8b5b825b08` | Compiled |
| Original three-dimensional Console2 | `99d6a1b962a8` | Compiled |
| Corrected three-dimensional Console2 | `a5308e8fb313` | Compiled |
| Corrected 200-step instrumented source | `73d1c52cd1d8` | Compiled |
| Literal 200-step instrumented source | `d221de76665e` | Compiled |
| Portable-angle 200-step instrumented source | `fbcf6409e665` | Compiled |

`manifest.json` records full hashes and every location sharing those exact bytes.
`compile_results.json` records each command, compiler version, source and
executable hashes, duration and exit code. `variants/*/compile.log` records
compiler output. Files with identical bytes were built once and every origin
remains listed in the manifest.

The build used `-std=legacy -ffree-line-length-none -O2 -g -fcheck=all
-fbacktrace -ffpe-trap=invalid,overflow`. No source portability patch was needed.
The preserved portable-angle version already contains the explicitly recorded
`atan2` substitution; it is a distinct recovered source, not a new change.

The compiler was installed only in `reference/toolchain/` by the recovered
`bootstrap_compiler.py`. All four pinned Debian packages passed their archived
SHA-256 checks. The system installation was not changed. GNU compiler binaries,
Debian packages, Python dependencies and generated executable files are local
build artifacts and need not be committed to the source repository.

To reproduce in this workspace:

```sh
python reference/bootstrap_compiler.py
python compile_all.py --toolchain reference/toolchain
PYTHONPATH=python-deps python reference/audit/test_native_operators.py
```

The last command requires NumPy and pytest. The seven native tests compile small
harnesses from exact original/corrected Fortran operator blocks. They cover
yielded rod force balance, preservation of accumulated transverse force,
unchanged pre-yield/compression behavior, interface force and work balance,
contacted-node reaction mass, constrained support velocity, and frictionless
cylinder tangential motion/reaction. Results and hashes are recorded in
`reference/audit/native_operator_verification.json`; stdout is `native_tests.log`.

These checks establish compilation and the tested operator behavior. They do
not establish the physical validity, convergence or benchmark equivalence of
the full peridynamic models. Historical simulation reports remain historical
evidence and were not represented as new runs.
