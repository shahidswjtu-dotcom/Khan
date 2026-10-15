import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from numba import njit
import time, csv, json

# Stage 19 notebook cell 120; preserved mathematical operations
def stage18_run(steps=425, dt=p.dt, save_every=2, mode='conservative', limit_mm=50, verbose=True):
    if mode not in ('conservative', 'source_density'):
        raise ValueError(mode)
    fast_concrete.reset()
    fast_interface.reset()
    appendix_c_beam.reset()
    uc = np.zeros_like(xyz)
    vc = np.zeros_like(xyz)
    us = np.zeros_like(steel_xyz)
    vs = np.zeros_like(steel_xyz)
    zh = float(hammer_z_initial)
    vh = float(HAMMER_V0)
    mid = (np.abs(xyz[:, 0] - p.L / 2) <= dx / 2 + 1e-12) & (xyz[:, 2] >= p.H / 2)

    def forces(commit):
        fnc, nb = fast_concrete.evaluate(uc, dt, commit)
        fns = rod_forces(us)
        if mode == 'conservative':
            fic, fis, ni = conservative_interface(uc, us, commit)
            ac = (fnc + fic) / CONCRETE_NODE_MASS
            a_s = (fns + fis) / STEEL_MASS[:, None]
            mom = np.linalg.norm(fic.sum(axis=0) + fis.sum(axis=0))
        else:
            fic, fis, ni = fast_interface.evaluate(uc, us, commit)
            ac = fnc / CONCRETE_NODE_MASS + fic / p.rho_concrete
            a_s = fns / STEEL_MASS[:, None] + fis / p.rho_steel
            mom = np.linalg.norm(fic.sum(axis=0) * V + np.sum(fis * STEEL_CONTROL_VOLUME[:, None], axis=0))
        c, hammer_f, _, n_contact, penetration = sphere_contact(xyz + uc, zh)
        ac += c / CONCRETE_NODE_MASS
        ah = float(hammer_f[2]) / HAMMER_MASS - (9.81 if USE_GRAVITY else 0.0)
        return (ac, a_s, ah, float(hammer_f[2]), int(n_contact), float(penetration), nb, ni, mom)
    ac, ast, ah, fh, nc, pen, nb, ni, resid = forces(False)
    data = []
    error = None
    first = None
    start = time.perf_counter()
    for k in range(1, steps + 1):
        vc += 0.5 * dt * ac
        vs += 0.5 * dt * ast
        vh += 0.5 * dt * ah
        uc += dt * vc
        us += dt * vs
        zh += dt * vh
        apply_supports(uc, vc)
        ac, ast, ah, fh, nc, pen, nb, ni, resid = forces(True)
        vc += 0.5 * dt * ac
        vs += 0.5 * dt * ast
        vh += 0.5 * dt * ah
        apply_supports(uc, vc)
        if nc > 0 and first is None:
            first = k
        maxc = float(np.max(np.linalg.norm(uc, axis=1))) * 1000.0
        maxs = float(np.max(np.linalg.norm(us, axis=1))) * 1000.0
        if not np.isfinite(maxc + maxs + fh) or max(maxc, maxs) > limit_mm or pen > dx / 2:
            error = f'safety stop {k}: concrete {maxc:.4g} mm; steel {maxs:.4g} mm'
        if k == 1 or k % save_every == 0 or k == first or error:
            data.append(dict(step=k, time_ms=k * dt * 1000.0, impact_kN=max(fh, 0.0) / 1000.0, midspan_mm=-1000.0 * np.mean(uc[mid, 2]), max_concrete_mm=maxc, max_steel_mm=maxs, n_contact=nc, broken_concrete=int(nb), broken_interface=int(ni), hammer_velocity_mps=vh, interface_physical_force_residual_N=resid))
        if verbose and (k % 100 == 0 or error):
            print('step', k, 'force_kN', round(fh / 1000, 3), 'max_c_mm', round(maxc, 6), 'max_s_mm', round(maxs, 6), 'interface residual N', round(resid, 6))
        if error:
            break
    out = STAGE18_DIR / f'{mode}_history.csv'
    with out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    summ = dict(mode=mode, steps_completed=k, first_contact=first, safety_stop=error, max_recorded_force_kN=max((x['impact_kN'] for x in data)), max_concrete_mm=max((x['max_concrete_mm'] for x in data)), max_steel_mm=max((x['max_steel_mm'] for x in data)), elapsed_seconds=round(time.perf_counter() - start, 2), interpretation='diagnostic only; interface normalization and material constants not Bowen-verified')
    (STAGE18_DIR / f'{mode}_summary.json').write_text(json.dumps(summ, indent=2))
    print(summ)
    return (data, summ)
