"""Install a pinned GNU Fortran into this workspace without root access.

Debian x86_64 only; requires the system GCC 14 runtime/linker and dpkg-deb.
Other platforms should install gfortran and use run_reference.py --compiler.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request

root = Path(__file__).resolve().parent
target = root / "toolchain"
target.mkdir(exist_ok=True)
for package in json.loads((root / "compiler_packages.json").read_text()):
    path = target / package["url"].rsplit("/", 1)[1]
    if not path.exists():
        urllib.request.urlretrieve(package["url"], path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != package["sha256"]:
        raise RuntimeError("Compiler package checksum mismatch: " + str(path))
    subprocess.run(["dpkg-deb", "-x", str(path), str(target)], check=True)
subprocess.run([str(target / "usr/bin/gfortran-14"), "--version"], check=True)
