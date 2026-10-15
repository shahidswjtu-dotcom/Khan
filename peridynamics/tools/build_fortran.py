#!/usr/bin/env python3
"""Compile every distinct recovered Fortran source, keeping all version aliases."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from verify_sources import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=shutil.which("gfortran"))
    parser.add_argument("--toolchain", type=Path,
                        help="Optional root of the archived Debian GNU 14 compiler extraction")
    parser.add_argument("--out", type=Path, default=Path(".build/fortran"))
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error("choose a new empty build directory")
    if out.is_relative_to(repo / "peridynamics/versions"):
        parser.error("build output must be outside the preserved versions")
    toolchain = args.toolchain.resolve() if args.toolchain else None
    compiler = args.compiler or (str(toolchain / "usr/bin/gfortran-14") if toolchain else None)
    if not compiler:
        parser.error("install gfortran or supply --compiler or --toolchain")
    integrity = verify(repo)
    if not integrity["passed"]:
        print(json.dumps(integrity, indent=2))
        return 1
    manifest = json.loads((repo / "peridynamics/RECOVERY_MANIFEST.json").read_text())
    variants = {}
    for record in manifest["files"]:
        if Path(record["destination"]).suffix.lower() in (".f", ".f90", ".f95", ".for"):
            variants.setdefault(record["sha256"], []).append(record["destination"])
    if not variants:
        parser.error("no recovered Fortran sources found")
    out.mkdir(parents=True, exist_ok=True)
    version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    flags = ["-std=legacy", "-ffree-line-length-none", "-O2", "-g", "-fcheck=all",
             "-fbacktrace", "-ffpe-trap=invalid,overflow"]
    if toolchain:
        # Same Debian x86_64 GCC-14 runtime prerequisites as bootstrap_compiler.py.
        flags += ["-B/usr/lib/gcc/x86_64-linux-gnu/14/",
                  "-L" + str(toolchain / "usr/lib/x86_64-linux-gnu"),
                  "-Wl,-rpath," + str(toolchain / "usr/lib/x86_64-linux-gnu")]
    results = []
    for sha, aliases in sorted(variants.items()):
        folder = out / sha[:12]
        folder.mkdir()
        binary = folder / "solver"
        command = [compiler, *flags, str(repo / aliases[0]), "-o", str(binary)]
        start = time.monotonic()
        with (folder / "compile.log").open("w") as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        record = {"source_sha256": sha, "source_paths": aliases, "command": command,
                  "exit_code": process.returncode,
                  "elapsed_seconds": round(time.monotonic() - start, 3),
                  "log": sha[:12] + "/compile.log"}
        if process.returncode == 0:
            record["executable_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest()
        results.append(record)
        print(json.dumps(record), flush=True)
    summary = {"compiler": version, "variants": results, "full_simulations_run": False,
               "source_integrity_after": verify(repo),
               "passed": all(r["exit_code"] == 0 for r in results)}
    summary["passed"] = summary["passed"] and summary["source_integrity_after"]["passed"]
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
