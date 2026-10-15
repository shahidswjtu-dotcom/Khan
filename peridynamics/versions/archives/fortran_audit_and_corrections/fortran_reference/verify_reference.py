"""Compare independent literal/portable runs without hiding nonfinite results."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("literal", type=Path)
    parser.add_argument("portable", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    ma = json.loads((args.literal / "manifest.json").read_text())
    mb = json.loads((args.portable / "manifest.json").read_text())
    assert ma["exit_code"] == mb["exit_code"] == 0
    assert ma["original_sha256"] == mb["original_sha256"]
    assert ma["steps"] == mb["steps"]
    assert not ma["portable_angles"] and mb["portable_angles"]
    history = [np.genfromtxt(p / "history.csv", delimiter=",", skip_header=1, ndmin=2)
               for p in (args.literal, args.portable)]
    state = [np.loadtxt(p / "final_state.txt", ndmin=2)
             for p in (args.literal, args.portable)]
    assert history[0].shape == history[1].shape == (ma["steps"], 6)
    assert state[0].shape == state[1].shape == (36398, 11)
    for a, b in (history, state):
        assert np.isfinite(a).all() and np.isfinite(b).all()
        np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-11)
    report = dict(original_sha256=ma["original_sha256"], steps=ma["steps"],
                  duration_s=ma["steps"] * ma["dt_s"],
                  all_finite=True, concrete_nodes=36036, steel_nodes=362,
                  history_max_absolute_difference_by_column=np.max(abs(history[0]-history[1]), axis=0).tolist(),
                  state_max_absolute_difference_by_column=np.max(abs(state[0]-state[1]), axis=0).tolist(),
                  full_loading_completed=ma["steps"] == 100000,
                  bowen_chapter3_equivalence=False)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print("Reference comparison passed")


if __name__ == "__main__":
    main()
