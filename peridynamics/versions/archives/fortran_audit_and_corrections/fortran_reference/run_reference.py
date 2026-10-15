"""Build and execute the supplied solver, with explicit, recorded instrumentation.

Short runs retain nt=100000 and its original prescribed displacement schedule.
No mechanics or material parameters are changed in literal mode.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from corrections import apply_corrections, add_diagnostics, AVAILABLE

ROOT = Path(__file__).resolve().parent


def instrument(source, steps, portable_angles=False):
    if not 1 <= steps <= 100000:
        raise ValueError("steps must be between 1 and 100000")
    assert source.count("do tt=1,nt") == 1
    source = source.replace("do tt=1,nt", f"do tt=1,{steps}")
    if portable_angles:
        old = "datan(dabs(coord(cnode,2)-coord(i,2))/dabs(coord(cnode,1)-coord(i,1)))"
        new = "datan2(dabs(coord(cnode,2)-coord(i,2)),dabs(coord(cnode,1)-coord(i,1)))"
        assert source.count(old) == 2
        source = source.replace(old, new)
    # Diagnostic files supplement the original outputs. Measurement is the
    # original pre-update single-node displacement, not a new bottom probe.
    source = source.replace("do tt=1,", 'open(91,file="history.csv",status="replace")\n'
        'write(91,\'(A)\') "step,time_s,probe_z_m,contact_x_N,contact_z_N,prescribed_z_next_m"\n'
        "do tt=1,", 1)
    source = source.replace("   if(tt.eq.nt)then", 'write(91,\'(I0,5(A,ES25.16E3))\') tt,",",tt*dt,",",dispdispld(tt),",",forceld(tt,1),",",forceld(tt,3),",",penuz\n'
        "   if(tt.eq.nt)then", 1)
    source = source.replace("end program main", 'close(91)\n'
        'open(92,file="final_state.txt",status="replace")\n'
        'do i=1,totnode\n'
        'write(92,\'(I0,10(1X,ES25.16E3))\') i,coord(i,:),disp(i,:),vel(i,:),dmg(i,1)\n'
        'enddo\nclose(92)\nend program main')
    return source


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--portable-angles", action="store_true",
                   help="Replace atan(y/x) with atan2(y,x); recorded portability correction")
    p.add_argument("--strict-zero", action="store_true",
                   help="Trap division by zero; original literal code fails on axis-aligned bonds")
    p.add_argument("--corrections", nargs="*", choices=sorted(AVAILABLE), default=[])
    p.add_argument("--diagnostics", action="store_true")
    p.add_argument("--compiler", default=shutil.which("gfortran") or str(ROOT / "toolchain/usr/bin/gfortran-14"))
    a = p.parse_args()
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    original = (ROOT / "Console2.original.f90").read_bytes()
    source = instrument(original.decode(), a.steps, a.portable_angles)
    source = apply_corrections(source, a.corrections)
    if a.diagnostics:
        source = add_diagnostics(source)
    generated = out / "solver.instrumented.f90"
    generated.write_text(source)
    flags = ["-std=legacy", "-ffree-line-length-none", "-O2", "-g", "-fcheck=all",
             "-fbacktrace", "-ffpe-trap=invalid,overflow" + (",zero" if a.strict_zero else "")]
    if str(ROOT / "toolchain") in a.compiler:
        flags += ["-B/usr/lib/gcc/x86_64-linux-gnu/14/",
                  "-L" + str(ROOT / "toolchain/usr/lib/x86_64-linux-gnu"),
                  "-Wl,-rpath," + str(ROOT / "toolchain/usr/lib/x86_64-linux-gnu")]
    command = [a.compiler, *flags, str(generated), "-o", str(out / "solver")]
    manifest = dict(original_sha256=hashlib.sha256(original).hexdigest(),
                    generated_sha256=hashlib.sha256(source.encode()).hexdigest(),
                    steps=a.steps, original_nt=100000, dt_s=3e-6,
                    portable_angles=a.portable_angles, command=command,
                    corrections=a.corrections, diagnostics=a.diagnostics,
                    compiler=subprocess.check_output([a.compiler, "--version"], text=True).splitlines()[0])
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    subprocess.run(command, check=True)
    with (out / "runtime.log").open("w") as log:
        result = subprocess.run([str(out / "solver")], cwd=out, stdout=log, stderr=subprocess.STDOUT)
    manifest["exit_code"] = result.returncode
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if result.returncode:
        raise SystemExit("Reference execution failed; inspect " + str(out / "runtime.log"))
    print("Reference run completed: " + str(out))


if __name__ == "__main__":
    main()
