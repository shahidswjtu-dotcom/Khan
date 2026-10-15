import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from numba import njit
import time, csv, json

# Stage 19 notebook cell 94; preserved mathematical operations
@njit(cache=False)
def _fast_concrete_kernel(X, U, I0, J0, r0, pair_weight, prev_stretch, broken, dt, E, ft, Gf, delta, fc, Vbond, commit):
    F = np.zeros((X.shape[0], 3), dtype=np.float64)
    failed_count = 0
    exp0 = 1.0 / (10.0 + 6.0 * fc / 10000000.0)
    theta = 10.0 ** (6.0 * exp0 - 2.0)
    reference_s0 = ft / E
    reference_sc = 5.0 * Gf / (6.0 * E * delta * reference_s0)
    base_c = 72.0 * E / (np.pi * delta ** 4)
    for k in range(I0.size):
        i, j = (I0[k], J0[k])
        ax = X[j, 0] + U[j, 0] - X[i, 0] - U[i, 0]
        ay = X[j, 1] + U[j, 1] - X[i, 1] - U[i, 1]
        az = X[j, 2] + U[j, 2] - X[i, 2] - U[i, 2]
        length = np.sqrt(ax * ax + ay * ay + az * az)
        if length <= 1e-15:
            raise ValueError('collapsed concrete bond')
        stretch = (length - r0[k]) / r0[k]
        rate = (stretch - prev_stretch[k]) / dt
        if rate <= 3e-06:
            dif_t = 1.0
            dif_e = 1.0
        else:
            ratio = rate / 3e-06
            dif_e = ratio ** 0.016
            if rate <= 1.0:
                dif_t = ratio ** exp0
            else:
                dif_t = theta * ratio ** (1.0 / 3.0)
        mult = dif_t / dif_e
        s0 = reference_s0 * mult
        sc = reference_sc * mult
        new_broken = broken[k] or stretch >= sc
        if new_broken:
            scalar = 0.0
            failed_count += 1
        else:
            scaled = r0[k] / delta
            c = base_c * (1.0 - scaled * scaled) ** 2 * dif_e
            if stretch <= s0:
                scalar = c * stretch
            else:
                scalar = c * s0 * (sc - stretch) / (sc - s0)
        amplitude = scalar * Vbond * pair_weight[k] / length
        fx = amplitude * ax
        fy = amplitude * ay
        fz = amplitude * az
        F[i, 0] += fx
        F[i, 1] += fy
        F[i, 2] += fz
        F[j, 0] -= fx
        F[j, 1] -= fy
        F[j, 2] -= fz
        if commit:
            prev_stretch[k] = stretch
            if new_broken:
                broken[k] = True
    return (F, failed_count)

