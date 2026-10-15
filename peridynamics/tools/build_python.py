#!/usr/bin/env python3
"""Build and test copies of the recovered Python projects without editing snapshots."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from verify_sources import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(".build/peridynamics"))
    parser.add_argument("--skip-tests", action="store_true",
                        help="Build only; use when the unchanged suites were already checked separately")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error("choose a new empty build directory")
    if out.is_relative_to(repo / "peridynamics/versions"):
        parser.error("build output must be outside the preserved versions")
    out.mkdir(parents=True, exist_ok=True)
    integrity = verify(repo)
    if not integrity["passed"]:
        print(json.dumps(integrity, indent=2))
        return 1
    env = os.environ.copy()
    for key, name in (("MPLCONFIGDIR", "matplotlib"), ("NUMBA_CACHE_DIR", "numba"),
                      ("XDG_CACHE_HOME", "xdg")):
        path = out / "cache" / name
        path.mkdir(parents=True, exist_ok=True)
        env[key] = str(path)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLBACKEND"] = "Agg"
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[key] = "1"
    results = []

    def run(name, command, cwd):
        started = time.monotonic()
        log = out / (name + ".log")
        with log.open("w") as stream:
            process = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                     stderr=subprocess.STDOUT, check=False)
        record = {"name": name, "command": command, "exit_code": process.returncode,
                  "elapsed_seconds": round(time.monotonic() - started, 3),
                  "log": log.name}
        results.append(record)
        print(json.dumps(record), flush=True)
        return process.returncode == 0

    fresh = None
    older = None
    for name in ("bowen_from_scratch", "bowen_reproduction"):
        source = repo / "peridynamics/versions/workspace" / name
        target = out / "sources" / name
        shutil.copytree(source, target)
        if name == "bowen_from_scratch":
            fresh = target
        else:
            older = target
        if not args.skip_tests:
            run(name + "_tests", [sys.executable, "-m", "pytest", "-q", "-p",
                "no:cacheprovider", "--junitxml=" + str(out / (name + "_tests.xml"))], target)

    wheels = out / "wheels"
    wheels.mkdir(exist_ok=True)
    run("wheel", [sys.executable, "-m", "pip", "wheel", "--no-deps",
        "--no-build-isolation", "--no-index", "--wheel-dir", str(wheels), "."], fresh)
    run("mesh", [sys.executable, "build.py"], fresh)
    run("stage19_saved_sample_regression", [sys.executable, "verify_regression.py"], older)
    for project, entrypoints in ((fresh, ("run.py",)),
                                 (older, ("run.py", "run_thesis.py", "run_printed.py"))):
        for entrypoint in entrypoints:
            run(project.name + "_" + entrypoint.replace(".", "_") + "_help",
                [sys.executable, entrypoint, "--help"], project)
    final_integrity = verify(repo)
    summary = {"python": sys.version, "source_integrity_before": integrity,
               "source_integrity_after": final_integrity, "checks": results,
               "test_suites_skipped": args.skip_tests,
               "passed": all(r["exit_code"] == 0 for r in results)
               and final_integrity["passed"],
               "scope": "Software builds and checks; no full impact simulation or experimental validation."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
