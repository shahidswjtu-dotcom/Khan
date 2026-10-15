from .geometry import *
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from numba import njit
import time, csv, json

# Stage 19 notebook cell 54; preserved mathematical operations
def sphere_contact(xyz_def, hammer_center_z, v_concrete=None, vz_hammer=0.0, k=CONTACT_STIFFNESS):
    center = np.array([hammer_xy[0], hammer_xy[1], hammer_center_z], dtype=float)
    d = xyz_def - center[None, :]
    distances = np.linalg.norm(d, axis=1)
    penetration = np.maximum(HAMMER_RADIUS - distances, 0.0)
    active = penetration > 0
    forces = np.zeros_like(xyz_def)
    if active.any():
        unit = d[active] / distances[active, None]
        magnitudes = k * penetration[active]
        forces[active] = magnitudes[:, None] * unit
    f_hammer = -np.sum(forces, axis=0)
    potential = 0.5 * k * np.sum(penetration[active] ** 2)
    return (forces, f_hammer, float(potential), int(active.sum()), float(penetration.max(initial=0.0)))
