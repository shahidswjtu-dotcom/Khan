# Repeating the builds

Run these commands from the repository root. Build products belong in `.build/`;
the preserved version folders are inputs and should remain unchanged.

## Python

Use Python 3.12 for the verified dependency versions. The recovery was checked
with Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Numba 0.68.0,
llvmlite 0.50.0, Matplotlib 3.10.8 and pytest 9.1.1.

```sh
python3.12 -m venv .build/venv
.build/venv/bin/python -m pip install -r peridynamics/versions/workspace/bowen_from_scratch/requirements.txt setuptools wheel
.build/venv/bin/python peridynamics/tools/verify_sources.py
.build/venv/bin/python peridynamics/tools/build_python.py --out .build/python-check
```

The build tool copies both Python projects into the output folder, runs their
existing test suites, builds the independent reconstruction wheel offline using
the installed build tools, builds the reference mesh, checks the saved Stage 19
sample regression, and checks the four main CLI entry points. Source hashes are
verified before and after. It records commands, logs and a JSON summary. Choose
a new output directory for each invocation. `--skip-tests` is available for a
build when the unchanged suites have already been verified separately.

For a new simulation, work from a copy under `.build/python-check/sources/`, read
that project's README and use a new output directory. The independent model's
driver explicitly requires selecting its documented interface hypothesis.
Full impact runs are separate from software builds and can be expensive.

## Fortran

GNU Fortran 14.2.0 compiled all six distinct recovered sources without edits:

```sh
python3 peridynamics/tools/build_fortran.py --compiler gfortran --out .build/fortran-check
```

The tool groups identical files by SHA-256, compiles every distinct source once,
and records all version aliases and the resulting executable hashes. It uses the
historical GNU legacy/free-line flags with bounds checking and floating-point
traps. It does not run the full solvers.

If GNU Fortran is unavailable on a compatible Debian x86_64 host, the recovered
bootstrap installs pinned packages into a local directory, with SHA-256 checks:

```sh
mkdir -p .build
cp -a peridynamics/versions/archives/fortran_audit_and_corrections/fortran_reference .build/fortran-reference
python3 .build/fortran-reference/bootstrap_compiler.py
python3 peridynamics/tools/build_fortran.py --toolchain .build/fortran-reference/toolchain --out .build/fortran-check
.build/venv/bin/python .build/fortran-reference/audit/test_native_operators.py
```

That bootstrap requires the system GCC 14 runtime/linker and `dpkg-deb`; use a
normal GNU Fortran installation on other platforms. Native operator tests need
NumPy and pytest. Their copied test driver uses the local bootstrapped compiler.

## Recorded verification

The fresh reconstruction has 22 existing tests and the earlier reproduction has
24. Seven native Fortran tests were also rerun. All 331 recovered Python files
(144 distinct contents) parse successfully, and every recovered notebook is
valid JSON. Full evidence is under [builds/2026-10-15/](builds/2026-10-15/).

These are software checks. The existing research audits document unresolved
physical assumptions and failure to reproduce the complete experimental
benchmark. The later `pdbeam_v4` package named in the PDF was not recovered.
