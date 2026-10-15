# Build and verification evidence — 2026-10-15

The recovered software builds successfully. No scientific model source was
changed to obtain these results.

| Check | Result |
| --- | --- |
| Recovered file integrity | All 473 source/document/fixture files match their recorded SHA-256 |
| Python source syntax | All 331 Python files parse (144 distinct contents) |
| Independent reconstruction | 22 existing tests pass |
| Earlier reproduction engines | 24 existing tests pass |
| Stage 19 saved-sample regression | All 112 saved samples match |
| Python distribution | Wheel built offline, installed separately, and smoke-tested |
| Mesh and installed-wheel smoke | 33,600 particles, 1,755,720 bonds; one integration step and exact checkpoint roundtrip |
| Fortran builds | All six distinct sources compile unchanged with GNU Fortran 14.2.0 |
| Native Fortran operator checks | Seven existing tests pass |
| Integrated build tools | Repository-relative Python and Fortran build tools pass |

The two Python suites were checked against the original projects and the
recovered repository copies. These are repeated checks of the same 46 Python
tests, not 92 distinct tests. Fortran source aliases with identical hashes were
compiled once per distinct source; all aliases are recorded.

The installable Python wheel is in [python/wheels/](python/wheels/). Python
commands, dependency versions, JUnit reports, source hashes and logs are under
[python/](python/). Fortran commands, compiler identity, source/executable hashes
and native-test results are under [fortran/](fortran/). The `integrated-*` folders
record fresh runs of the portable repository build tools; their Python run used
`--skip-tests` because the unchanged recovered suites had already passed.
`ARTIFACTS.json` records hashes of the collected evidence files.

The evidence includes absolute paths from this attached workspace as execution
provenance. Use [the repository build instructions](../../BUILD.md) to reproduce
the checks in another checkout. Local compiler packages, executable binaries,
large simulation outputs, caches and test checkpoints are not included here.

No full impact simulation or experimental validation was performed during this
recovery. The later `pdbeam_v4` package referenced by the existing PDF remains
unavailable; the preserved earlier code is labelled accordingly.
