# Python build and software verification

Build and software verification of recovered Oct 9 Python projects; no full impact rerun or scientific benchmark validation.

Original source files were not edited. Caches, temporary files, builds and smoke checkpoints stayed under `work/python-build`.

| Check | Result | Duration |
|---|---|---|
| fresh-cli | Passed | 0.925 s |
| fresh-mesh | Passed | 1.922 s |
| fresh-tests | 22 passed, 0 failed, 0 errors | 13.031 s |
| fresh-wheel | Passed | 1.218 s |
| printed-cli | Passed | 1.095 s |
| recovered-fresh-tests | 22 passed, 0 failed, 0 errors | 18.572 s |
| recovered-regression | Passed | 0.165 s |
| recovered-reproduction-tests | 24 passed, 0 failed, 0 errors | 83.829 s |
| reproduction-cli | Passed | 1.196 s |
| reproduction-tests | 24 passed, 0 failed, 0 errors | 84.933 s |
| thesis-cli | Passed | 1.005 s |
| wheel-install | Passed | 0.265 s |
| wheel-smoke | Passed | 11.857 s |

The installed wheel assembled 33,600 concrete particles, 270 steel nodes and 1,755,720 bonds, advanced one step, and passed an exact checkpoint roundtrip. The mesh starts at zero force. The recovered Stage 19 fixture regression passes all 112 retained samples.

## Reproduce in the attached workspace

```bash
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py fresh-tests
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py reproduction-tests
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py recovered-fresh-tests
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py recovered-reproduction-tests
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py regression
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py wheel
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py mesh
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py cli
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py wheel-smoke
/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/collect_summary.py
```

The driver records every executed command, working directory, environment override, exit status and duration in `results/*.json`. JUnit XML and detailed test durations accompany the logs. The build command is `python -m pip wheel --no-deps --no-build-isolation --no-index --wheel-dir <wheel-directory> <source-copy>`.

Python 3.12.14. Runtime dependencies match the pinned requirements; detailed versions and hashes are in `build-summary.json`.

These checks establish executable software and checkpoint consistency. The models remain declared reconstructions; they do not establish Bowen author equivalence, full-event agreement or experimental validation.
