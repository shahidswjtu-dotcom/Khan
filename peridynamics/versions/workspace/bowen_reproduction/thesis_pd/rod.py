from .geometry import *
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from numba import njit
import time, csv, json

# Stage 19 notebook cell 44; preserved mathematical operations
def rod_forces(ust):
    us = np.asarray(ust, float)
    d0 = steel_xyz[rod_pairs[:, 1]] - steel_xyz[rod_pairs[:, 0]]
    unit = d0 / steel_lengths[:, None]
    ext = np.sum((us[rod_pairs[:, 1]] - us[rod_pairs[:, 0]]) * unit, axis=1)
    fs = (steel_E * steel_A / steel_lengths * ext)[:, None] * unit
    out = np.zeros_like(us)
    np.add.at(out, rod_pairs[:, 0], fs)
    np.add.at(out, rod_pairs[:, 1], -fs)
    return out

def coupling_forces(ucon, ust, k=COUPLING_STIFFNESS):
    uc = np.sum(np.asarray(ucon)[sample_idx] * sample_w[:, :, None], axis=1)
    rel = uc - np.asarray(ust)[sample_steel]
    fsteel_sample = k * rel
    fs = np.zeros_like(ust)
    fc = np.zeros_like(ucon)
    np.add.at(fs, sample_steel, fsteel_sample)
    for a in range(8):
        np.add.at(fc, sample_idx[:, a], -sample_w[:, a, None] * fsteel_sample)
    return (fc, fs)
