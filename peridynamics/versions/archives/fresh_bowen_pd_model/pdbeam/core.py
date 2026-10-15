from dataclasses import dataclass, asdict
import json
import numpy as np
from .checkpoints import atomic_npz
from scipy.spatial import cKDTree

@dataclass(frozen=True)
class Benchmark:
    length: float = .600
    width: float = .070
    height: float = .100
    spacing: float = .005
    horizon: float = .015
    dt: float = 3e-7
    young: float = 25e9
    density: float = 2400.
    tensile_strength: float = 4.6e6
    fracture_energy: float = 105.
    steel_young: float = 200e9
    steel_density: float = 7750.
    steel_yield: float = 482e6
    bar_diameter: float = .006
    bar_height: float = .020
    hammer_mass: float = 38.4
    hammer_speed: float = 10.8
    hammer_radius: float = .060

    @property
    def onset(self):
        return self.tensile_strength / self.young

    @property
    def terminal(self):
        return 5*self.fracture_energy/(6*self.young*self.horizon*self.onset)


def particles(p):
    axes = []
    for size in (p.length, p.width, p.height):
        n = round(size/p.spacing)
        if not np.isclose(n*p.spacing, size):
            raise ValueError('Dimensions must be integral multiples of spacing')
        axes.append((np.arange(n)+.5)*p.spacing)
    return np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)

@dataclass
class Bonds:
    ends: np.ndarray
    reference: np.ndarray
    length: np.ndarray
    stiffness: np.ndarray
    volume_fraction: np.ndarray


def family(x, p):
    ends = cKDTree(x).query_pairs(p.horizon*(1+1e-12), output_type='ndarray')
    ends = ends[np.lexsort((ends[:,1], ends[:,0]))]
    reference = x[ends[:,1]]-x[ends[:,0]]
    length = np.linalg.norm(reference, axis=1)
    # Eq. 2-12: force-density micromodulus [N/m^6].
    stiffness = 72*p.young/(np.pi*p.horizon**4)*np.maximum(0, 1-(length/p.horizon)**2)**2
    fraction = np.clip((p.horizon+.5*p.spacing-length)/p.spacing, 0, 1)
    return Bonds(ends, reference, length, stiffness, fraction)


def concrete_force(x, u, bonds, p, alive, history='irreversible'):
    """Return physical nodal force [N] and updated flags, never force density.

    Static Eq. 2-5 only. Rate and boundary corrections must be assembled
    explicitly before this operator can represent the full impact experiment.
    """
    if history not in ('irreversible', 'instantaneous'):
        raise ValueError('Specify irreversible or instantaneous history')
    i, j = bonds.ends.T
    d = x[j]+u[j]-x[i]-u[i]
    r = np.linalg.norm(d, axis=1)
    if np.any(r == 0):
        raise ValueError('Collapsed bond requires short-range contact handling')
    s = (r-bonds.length)/bonds.length
    law = np.where(s <= p.onset, s,
                   p.onset*np.maximum(p.terminal-s, 0)/(p.terminal-p.onset))
    active = alive if history == 'irreversible' else np.ones_like(alive)
    magnitude = bonds.stiffness*law*p.spacing**6*bonds.volume_fraction*active
    pair = magnitude[:,None]*d/r[:,None]
    f = np.zeros_like(x)
    np.add.at(f, i, pair)
    np.add.at(f, j, -pair)
    return f, alive & (s <= p.terminal)


def interpolate_weights(point, vertices):
    """Eight-node axis-aligned trilinear stencil, Appendix C."""
    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    if np.any(hi <= lo) or np.any(point < lo-1e-12) or np.any(point > hi+1e-12):
        raise ValueError('Point must lie in a nondegenerate interpolation cell')
    q = np.clip((point-lo)/(hi-lo), 0, 1)
    bits = np.isclose(vertices, hi).astype(int)
    if len(set(map(tuple, bits))) != 8 or not np.all(np.isclose(vertices, lo)|np.isclose(vertices, hi)):
        raise ValueError('Require all eight distinct cell corners')
    return np.prod(np.where(bits, q, 1-q), axis=1)


def interface_pair(steel_position, concrete_positions, weights, reference_length,
                   stiffness, steel_volume, sample_volume):
    """Elastic PDROD sample force transfer in physical N.

    Caller must supply scientifically justified interface stiffness and volumes.
    Synthetic verification values are not Chapter 3 parameter definitions.
    """
    virtual = weights @ concrete_positions
    d = virtual-steel_position
    r = np.linalg.norm(d)
    if r == 0 or reference_length <= 0 or min(steel_volume, sample_volume) <= 0:
        raise ValueError('Invalid interface geometry or control volumes')
    force = stiffness*(r/reference_length-1)*steel_volume*sample_volume*d/r
    return force, -weights[:,None]*force


def project_contact(old_positions, predicted_positions, center, radius, masses, dt):
    """Eq. 2-45--47: sphere projection and nodal impulse.

    Center is the predicted center. Hammer recoil is applied by the caller.
    """
    d = predicted_positions-center
    distance = np.linalg.norm(d, axis=1)
    inside = distance < radius
    if np.any(inside & (distance == 0)):
        raise ValueError('Contact normal undefined at sphere center')
    corrected = predicted_positions.copy()
    corrected[inside] = center+radius*d[inside]/distance[inside,None]
    velocity = (corrected-old_positions)/dt
    impulse = masses[:,None]*(corrected-predicted_positions)/dt
    return corrected, velocity, -impulse.sum(axis=0)


def save_checkpoint(path, p, step, **state):
    atomic_npz(path, config=json.dumps(asdict(p), sort_keys=True),
                        step=np.int64(step), **state)


def load_checkpoint(path, p):
    with np.load(path, allow_pickle=False) as data:
        if str(data['config']) != json.dumps(asdict(p), sort_keys=True):
            raise ValueError('Checkpoint configuration differs')
        return int(data['step']), {k: data[k].copy() for k in data.files if k not in ('config','step')}
