#!/usr/bin/env python3
"""Check every recovered source byte against its recorded SHA-256."""
import hashlib
import json
from pathlib import Path
import sys


def verify(root):
    root = Path(root).resolve()
    manifest = json.loads((root / "peridynamics/RECOVERY_MANIFEST.json").read_text())
    failures = []
    destinations = set()
    for record in manifest["files"]:
        relative = record["destination"]
        target = (root / relative).resolve()
        if not target.is_relative_to(root) or relative in destinations:
            failures.append({"path": relative, "error": "invalid or duplicate path"})
            continue
        destinations.add(relative)
        if not target.is_file():
            failures.append({"path": relative, "error": "missing file"})
            continue
        data = target.read_bytes()
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            failures.append({"path": relative, "error": "SHA-256 mismatch"})
    return {"files_checked": len(destinations), "failures": failures,
            "passed": not failures}


if __name__ == "__main__":
    result = verify(Path(__file__).resolve().parents[2])
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["passed"] else 1)
