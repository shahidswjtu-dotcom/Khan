import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from numba import njit
import time, csv, json

# Stage 19 notebook cell 3; preserved mathematical operations
@dataclass(frozen=True)
class Parameters:
    L: float = 0.6
    W: float = 0.07
    H: float = 0.1
    dx: float = 0.005
    horizon_ratio: float = 3.0
    support_span: float = 0.5
    E_concrete: float = 25000000000.0
    nu_concrete: float = 0.25
    rho_concrete: float = 2400.0
    ft_concrete: float = 4600000.0
    fc_concrete: float = 60000000.0
    Gf_concrete: float = 105.0
    E_steel: float = 200000000000.0
    rho_steel: float = 7750.0
    fy_steel: float = 482000000.0
    bar_diameter: float = 0.006
    n_bars: int = 2
    impactor_mass: float = 38.4
    impactor_speed: float = 10.8
    impactor_radius: float = 0.06
    dt: float = 3e-07
p = Parameters()
L, W, H, dx = (p.L, p.W, p.H, p.dx)
Nx, Ny, Nz = [round(t / dx) for t in (L, W, H)]
assert all(np.isclose(np.array([Nx, Ny, Nz]) * dx, [L, W, H]))
assert (Nx, Ny, Nz) == (120, 14, 20)

# Stage 19 notebook cell 5; preserved mathematical operations
x = (np.arange(Nx) + 0.5) * dx
y = (np.arange(Ny) + 0.5) * dx
z = (np.arange(Nz) + 0.5) * dx
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
xyz = np.column_stack((X.ravel(), Y.ravel(), Z.ravel())).astype(np.float64)
N = len(xyz)
Vpoint = dx ** 3
assert N == 33600
assert np.unique(xyz, axis=0).shape[0] == N
assert np.all(xyz.min(axis=0) > 0) and np.all(xyz.max(axis=0) < [L, W, H])
assert np.isclose(N * Vpoint, L * W * H)

# Stage 19 notebook cell 7; preserved mathematical operations
support_x = np.array([(p.L - p.support_span) / 2, (p.L + p.support_span) / 2])
bar_z = 0.01
bar_y = np.array([0.02, 0.05])
bar_x = np.array([0.02, p.L - 0.02])
bar_x_nodes = np.linspace(bar_x[0], bar_x[1], round((bar_x[1] - bar_x[0]) / dx) + 1)
steel_xyz = np.array([[xx, yy, bar_z] for yy in bar_y for xx in bar_x_nodes])
impact_xy = np.array([p.L / 2, p.W / 2])
illustrative_gap = 0.001
impact_center = np.array([impact_xy[0], impact_xy[1], p.H + p.impactor_radius + illustrative_gap])
assert np.isclose(np.diff(support_x)[0], 0.5)
assert len(bar_y) == p.n_bars
assert np.all((bar_y > 0) & (bar_y < W))
assert bar_z > p.bar_diameter / 2

# Stage 19 notebook cell 11; preserved mathematical operations
delta = p.horizon_ratio * dx
s0 = p.ft_concrete / p.E_concrete
sc = 5 * p.Gf_concrete / (6 * p.E_concrete * delta * s0)
c_constant = 12 * p.E_concrete / (np.pi * delta ** 4)

def c_weighted(r):
    r = np.asarray(r)
    return np.where(r <= delta, 72 * p.E_concrete / (np.pi * delta ** 4) * (1 - (r / delta) ** 2) ** 2, 0.0)

# Stage 19 notebook cell 29; preserved mathematical operations
from scipy.spatial import cKDTree
from time import perf_counter
try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda fun: fun
BUILD_FULL_NETWORK = True
if BUILD_FULL_NETWORK:
    t0 = perf_counter()
    tree = cKDTree(xyz)
    ij = tree.query_pairs(r=delta * (1.0 + 1e-10), output_type='ndarray').astype(np.int32)
    I, J = (ij[:, 0].copy(), ij[:, 1].copy())
    R0 = xyz[J] - xyz[I]
    LEN = np.linalg.norm(R0, axis=1)
    UNIT = R0 / LEN[:, None]
    counts = np.bincount(np.concatenate([I, J]), minlength=N)
    print(f'Network: N={N:,} particles, {len(I):,} undirected bonds; neighbors min/mean/max={counts.min()}/{counts.mean():.2f}/{counts.max()}')
    print(f'Build time: {perf_counter() - t0:.2f} s; bond length range {LEN.min() * 1000.0:.2f}-{LEN.max() * 1000.0:.2f} mm')
    assert N == 33600 and I.size > 500000
    assert np.all((LEN > 0) & (LEN <= delta * (1.0 + 1e-09)))
    assert np.all(I < J) and np.all(counts > 0)

# Stage 19 notebook cell 31; preserved mathematical operations
V = dx ** 3
USE_PARTIAL_VOLUME = True
USE_SURFACE_CORRECTION = False
volweight = np.clip((delta + dx / 2 - LEN) / dx, 0.0, 1.0) if USE_PARTIAL_VOLUME else np.ones_like(LEN)
occupancy = np.bincount(I, weights=volweight, minlength=N) + np.bincount(J, weights=volweight, minlength=N)
interior = (xyz[:, 0] > delta + dx) & (xyz[:, 0] < p.L - delta - dx) & (xyz[:, 1] > delta + dx) & (xyz[:, 1] < p.W - delta - dx) & (xyz[:, 2] > delta + dx) & (xyz[:, 2] < p.H - delta - dx)
reference_neighbors = float(np.median(occupancy[interior]))
node_correction = np.ones(N)
if USE_SURFACE_CORRECTION:
    node_correction = np.clip(reference_neighbors / np.maximum(occupancy, 1), 1.0, 3.0)
PAIR_WEIGHT = volweight * np.sqrt(node_correction[I] * node_correction[J])
assert interior.any()

# Stage 19 notebook cell 40; preserved mathematical operations
steel_xyz = np.asarray(steel_xyz, dtype=np.float64)
Nsteel = len(steel_xyz)
rod_pairs = np.concatenate([np.column_stack([np.arange(k * len(bar_x_nodes), (k + 1) * len(bar_x_nodes) - 1), np.arange(k * len(bar_x_nodes) + 1, (k + 1) * len(bar_x_nodes))]) for k in range(len(bar_y))]).astype(np.int32)
steel_E = 200000000000.0
steel_diameter = 0.006
steel_A = np.pi * steel_diameter ** 2 / 4
steel_lengths = np.linalg.norm(steel_xyz[rod_pairs[:, 1]] - steel_xyz[rod_pairs[:, 0]], axis=1)
assert np.all(steel_lengths > 0)

# Stage 19 notebook cell 42; preserved mathematical operations
def concrete_grid_eight(point):
    point = np.asarray(point, dtype=float)
    nx, ny, nz = (Nx, Ny, Nz)
    h = dx
    a = point / h - 0.5
    base = np.floor(a).astype(int)
    if np.any(base < 0) or np.any(base >= np.array([nx - 1, ny - 1, nz - 1])):
        raise ValueError(f'Coupling point {point} outside center-grid interpolation domain')
    frac = a - base
    inds = []
    ws = []
    for qx in (0, 1):
        for qy in (0, 1):
            for qz in (0, 1):
                inds.append(np.ravel_multi_index((base[0] + qx, base[1] + qy, base[2] + qz), (nx, ny, nz)))
                ws.append((frac[0] if qx else 1 - frac[0]) * (frac[1] if qy else 1 - frac[1]) * (frac[2] if qz else 1 - frac[2]))
    return (np.array(inds, dtype=np.int32), np.array(ws, dtype=float))
COUPLING_OFFSET = dx
COUPLING_STIFFNESS = 200000000.0
sample_steel = []
sample_pos = []
sample_idx = []
sample_w = []
for si, pos in enumerate(steel_xyz):
    for sign in (-1.0, 1.0):
        pt = pos + np.array([0.0, sign * COUPLING_OFFSET, 0.0])
        idx, w = concrete_grid_eight(pt)
        sample_steel.append(si)
        sample_pos.append(pt)
        sample_idx.append(idx)
        sample_w.append(w)
sample_steel = np.asarray(sample_steel, dtype=np.int32)
sample_pos = np.asarray(sample_pos, dtype=float)
sample_idx = np.asarray(sample_idx, dtype=np.int32)
sample_w = np.asarray(sample_w, dtype=float)
assert np.all(sample_w >= -1e-14)
np.testing.assert_allclose(sample_w.sum(axis=1), 1.0, atol=1e-13)
np.testing.assert_allclose((xyz[sample_idx] * sample_w[:, :, None]).sum(axis=1), sample_pos, atol=1e-13)

# Stage 19 notebook cell 56; preserved mathematical operations
CONCRETE_NODE_MASS = p.rho_concrete * dx ** 3
STEEL_MASS = np.zeros(Nsteel)
rod_mass_per_element = p.rho_steel * steel_A * steel_lengths
np.add.at(STEEL_MASS, rod_pairs[:, 0], 0.5 * rod_mass_per_element)
np.add.at(STEEL_MASS, rod_pairs[:, 1], 0.5 * rod_mass_per_element)
assert np.all(STEEL_MASS > 0)

# Stage 19 notebook cell 52; preserved mathematical operations
import time
HAMMER_MASS = p.impactor_mass
HAMMER_RADIUS = p.impactor_radius
HAMMER_V0 = -p.impactor_speed
INITIAL_GAP = 0.001
hammer_z_initial = float(np.max(xyz[:, 2]) + HAMMER_RADIUS + INITIAL_GAP)
hammer_xy = np.array([p.L / 2, p.W / 2])
CONTACT_STIFFNESS = 200000000.0
CONTACT_DAMPING = 0.0
USE_GRAVITY = False
support_at_x = np.min(np.abs(xyz[:, 0, None] - support_x[None, :]), axis=1) <= dx / 2 + 1e-12
support_at_z = xyz[:, 2] <= dx / 2 + 1e-12
support_mask = support_at_x & support_at_z
SUPPORT_DOF_Z_ONLY = True
anchor_candidates = np.flatnonzero(support_mask)
assert len(anchor_candidates) >= 2
anchor_x = anchor_candidates[0]
anchor_y = anchor_candidates[1]
assert abs(HAMMER_V0) == 10.8 and len(anchor_candidates) > 0
