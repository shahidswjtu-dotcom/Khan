"""Read-only postprocessing of the completed supplied-Fortran run.

Snapshot displacements have the original E12.5 precision. Final displacements
and velocities use the instrumented ES25.16E3 dump. Mechanics are not modified.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "literal_100000"
OUT = Path(__file__).resolve().parent
N_STEEL = 362
DT = 3.0e-6
SY = 2.0e-3


def rods(steel_positions: np.ndarray) -> np.ndarray:
    return np.concatenate((np.diff(steel_positions[:181], axis=0),
                           np.diff(steel_positions[181:], axis=0)))


def group_state(a: np.ndarray) -> dict:
    # columns: node ID, X,Y,Z, uX,uY,uZ, vX,vY,vZ, damage.
    u = a[:, 4:7]
    v = a[:, 7:10]
    magnitude = np.linalg.norm(u, axis=1)
    i = int(np.argmax(magnitude))
    return {"nodes": len(a), "maximum_displacement_m": float(magnitude[i]),
            "maximum_node_one_based": int(a[i, 0]),
            "maximum_node_reference_xyz_m": a[i, 1:4].tolist(),
            "maximum_node_displacement_xyz_m": u[i].tolist(),
            "maximum_node_velocity_xyz_m_per_s": v[i].tolist(),
            "minimum_displacement_xyz_m": u.min(axis=0).tolist(),
            "maximum_displacement_xyz_m": u.max(axis=0).tolist(),
            "minimum_velocity_xyz_m_per_s": v.min(axis=0).tolist(),
            "maximum_velocity_xyz_m_per_s": v.max(axis=0).tolist(),
            "nodes_over_100mm": int(np.count_nonzero(magnitude > 0.1)),
            "mean_damage": float(np.mean(a[:, 10]))}


def main() -> None:
    final = np.loadtxt(RUN / "final_state.txt")
    history = np.genfromtxt(RUN / "history.csv", names=True, delimiter=",")
    assert final.shape == (36398, 11)
    assert np.array_equal(final[:, 0], np.arange(1, 36399))
    assert len(history) == 100000
    assert np.all(np.isfinite(final))

    rows = []
    for path in sorted(RUN.glob("coord_disp_pd_*.txt"),
                       key=lambda p: int(p.stem.rsplit("_", 1)[1])):
        step = int(path.stem.rsplit("_", 1)[1])
        a = np.loadtxt(path)
        assert a.shape == (36398, 8)
        u = a[:, 3:6]
        # Recompute after-contact norms from the vector components. The source
        # scalar displacement column was evaluated before contact projection.
        norms = np.linalg.norm(u, axis=1)
        current_rods = rods(final[:N_STEEL, 1:4] + u[:N_STEEL])
        stretch = np.linalg.norm(current_rods, axis=1) / 0.019 - 1.0
        force_window = history[:step]
        peak = int(np.argmax(np.abs(force_window["contact_z_N"])))
        rows.append({"step": step, "time_s": step * DT,
                     "steel_max_displacement_m": float(norms[:N_STEEL].max()),
                     "concrete_max_displacement_m": float(norms[N_STEEL:].max()),
                     "steel_max_abs_x_m": float(np.abs(u[:N_STEEL, 0]).max()),
                     "steel_max_abs_y_m": float(np.abs(u[:N_STEEL, 1]).max()),
                     "steel_max_abs_z_m": float(np.abs(u[:N_STEEL, 2]).max()),
                     "maximum_rod_stretch": float(stretch.max()),
                     "rods_above_sy": int(np.count_nonzero(stretch > SY)),
                     "steel_mean_damage": float(a[:N_STEEL, 7].mean()),
                     "concrete_mean_damage": float(a[N_STEEL:, 7].mean()),
                     "max_abs_contact_z_N_so_far": float(abs(force_window["contact_z_N"][peak])),
                     "contact_peak_step_so_far": int(force_window["step"][peak])})
    with (OUT / "trace_snapshots.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fy = 416.0 * 530.66  # MPa times mm^2 = N, exactly as printed.
    volume = 0.019 * 530.66e-6
    force_density = fy / volume
    damping = 3.0e6
    predicted_terminal = force_density / damping
    # Replay only the all-yielded rod branch in its actual element order.
    # The right endpoint's y/z are assigned from its current axial accumulator.
    # Its axial accumulator is subsequently cancelled by the next element, but
    # those previously copied y/z entries persist.
    defective_rod_density = np.zeros((N_STEEL, 3))
    for bar_start in [0, 181]:
        for left in range(bar_start, bar_start + 180):
            right = left + 1
            defective_rod_density[left, 0] += force_density
            defective_rod_density[right, 0] -= force_density
            defective_rod_density[right, 1] = defective_rod_density[right, 0]
            defective_rod_density[right, 2] = defective_rod_density[right, 0]
    expected_velocity = defective_rod_density / damping
    measured_velocity = final[:N_STEEL, 7:10]
    stretched = np.linalg.norm(rods(final[:N_STEEL, 1:4] + final[:N_STEEL, 4:7]), axis=1) / 0.019 - 1.0
    pk = int(np.argmax(np.abs(history["contact_z_N"])))
    def_pk = int(np.argmax(np.abs(history["probe_z_m"])))
    summary = {"source": str(RUN), "final_dump_precision": "ES25.16E3",
               "snapshot_precision": "E12.5; rod stretches reconstructed from vector components",
               "steel": group_state(final[:N_STEEL]),
               "concrete": group_state(final[N_STEEL:]),
               "sampled_first_rod_above_sy_step": next(r["step"] for r in rows if r["rods_above_sy"]),
               "sampled_last_no_rod_above_sy_step": max(r["step"] for r in rows if not r["rods_above_sy"]),
               "first_actual_yield_step": None,
               "first_actual_yield_caveat": "Not recorded; snapshot threshold crossing brackets the sampled transition, not first execution of that branch.",
               "contact_peak": {"step": int(history["step"][pk]), "time_s": float(history["time_s"][pk]),
                                "force_z_N": float(history["contact_z_N"][pk])},
               "probe_peak": {"step": int(history["step"][def_pk]), "time_s": float(history["time_s"][def_pk]),
                              "displacement_z_m": float(history["probe_z_m"][def_pk])},
               "final_rods_above_sy": int(np.count_nonzero(stretched > SY)),
               "yield_force_N": fy, "steel_nodal_volume_m3": volume,
               "yield_force_density_N_per_m3": force_density,
               "predicted_terminal_velocity_m_per_s": predicted_terminal,
               "max_abs_final_velocity_minus_defective_branch_terminal_m_per_s": float(np.abs(measured_velocity - expected_velocity).max()),
               "defective_rod_net_force_xyz_N": (defective_rod_density.sum(axis=0) * volume).tolist(),
               "final_steel_kinetic_energy_J_using_source_nodal_mass": float(0.5 * 7850 * volume * np.sum(measured_velocity**2)),
               "final_concrete_kinetic_energy_J_using_source_nodal_mass": float(0.5 * 2500 * 0.019**3 * np.sum(final[N_STEEL:, 7:10]**2))}
    (OUT / "trace_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    t = np.array([r["time_s"] for r in rows])
    for name, label in [("steel", "Steel"), ("concrete", "Concrete")]:
        axes[0, 0].semilogy(t, [r[f"{name}_max_displacement_m"] for r in rows], "o-", label=label)
    axes[0, 0].set(ylabel="Maximum displacement magnitude (m)", xlabel="Time (s)")
    axes[0, 0].legend()
    for direction in "xyz":
        axes[0, 1].semilogy(t, [r[f"steel_max_abs_{direction}_m"] for r in rows], "o-", label=direction.upper())
    axes[0, 1].set(ylabel="Maximum absolute steel displacement (m)", xlabel="Time (s)")
    axes[0, 1].legend()
    axes[1, 0].plot(t, [r["rods_above_sy"] for r in rows], "o-")
    axes[1, 0].set(ylabel="Rods currently above sy = 0.002", xlabel="Time (s)", ylim=(-5, 365))
    mask = (history["step"] >= 62000) & (history["step"] <= 67000)
    axes[1, 1].plot(history["time_s"][mask], history["contact_z_N"][mask] / 1000.0)
    axes[1, 1].set(ylabel="Contact vertical force (kN)", xlabel="Time (s)")
    axes[1, 1].axvline(history["time_s"][pk], color="red", linestyle="--", label=f"Peak step {int(history['step'][pk])}")
    axes[1, 1].legend()
    for ax in axes.flat:
        ax.grid(alpha=0.3)
    fig.suptitle("Completed literal Fortran run: sampled transition and artificial rod drift")
    fig.savefig(OUT / "trace_failure.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
