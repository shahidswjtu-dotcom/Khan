"""Monitor the already-running native solver and summarize on completion."""
import json
import os
from pathlib import Path
import subprocess
import time
import argparse
root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run", nargs="?", type=Path, default=root / "runs/literal_100000")
run = parser.parse_args().run.resolve()
while True:
    try:
        m = json.loads((run / "manifest.json").read_text())
    except (OSError, json.JSONDecodeError):
        time.sleep(10)
        continue
    if "exit_code" in m:
        if m["exit_code"] != 0:
            (run / "progress.json").write_text(json.dumps({"status":"failed", "solver_exit_code":m["exit_code"]},indent=2)+"\n")
            raise SystemExit("Solver failed; inspect runtime.log")
        subprocess.run([str(root.parent / ".venv/bin/python"), str(root / "summarize_run.py"), str(run)], check=True)
        (run / "progress.json").write_text(json.dumps({"status":"completed", "completed_steps":m["steps"],"target_steps":m["steps"],"solver_exit_code":0},indent=2)+"\n")
        break
    history = run / "history.csv"
    last_step = 0
    if history.exists():
        with history.open("rb") as stream:
            stream.seek(max(0, history.stat().st_size - 4096))
            data = stream.read()
            lines = data[:data.rfind(b"\n")+1].splitlines()
            for line in reversed(lines):
                try:
                    last_step = int(line.split(b",")[0])
                    break
                except ValueError:
                    pass
    (run / "progress.json").write_text(json.dumps({"status":"running", "completed_steps_at_least":last_step,"target_steps":m["steps"],"updated_unix_s":time.time()},indent=2)+"\n")
    time.sleep(10)
