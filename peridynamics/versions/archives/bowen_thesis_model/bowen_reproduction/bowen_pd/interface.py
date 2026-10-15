import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
from numba import njit
import time, csv, json

# Stage 19 notebook cell 74; preserved mathematical operations
@dataclass(frozen=True)
class AppendixCInterfaceParameters:
    bcc: float
    s01: float
    sc1: float
    horizon: float
    quadrature_volume: float

# Stage 19 notebook cell 75; preserved mathematical operations
class AppendixCPDROD:
    """Source-shaped translation of the Appendix C nonlinear steel-to-concrete force.

    Sample is a virtual material point whose displacement is trilinearly interpolated
    from eight concrete points. Forces on the real concrete nodes are -N_a*F_steel.
    Persistent broken flags prevent crack healing on subsequent evaluations.
    """

    def __init__(self, steel_reference, virtual_reference, indices, weights, steel_ids, params):
        self.xs = np.asarray(steel_reference, dtype=float)
        self.xv = np.asarray(virtual_reference, dtype=float)
        self.indices = np.asarray(indices, dtype=np.int32)
        self.weights = np.asarray(weights, dtype=float)
        self.steel_ids = np.asarray(steel_ids, dtype=np.int32)
        self.params = params
        self.reference = self.xv - self.xs[self.steel_ids]
        self.lengths = np.linalg.norm(self.reference, axis=1)
        if self.indices.ndim != 2 or self.indices.shape[1] != 8:
            raise ValueError('eight interpolation weights required')
        if np.any(self.lengths <= 0):
            raise ValueError('zero-length steel–concrete interface bond')
        if np.any(self.lengths > params.horizon*(1+1e-12)):
            raise ValueError('all coupling samples must lie inside horizon')
        np.testing.assert_allclose(self.weights.sum(axis=1), 1.0, atol=1e-12)
        if not 0 < params.s01 < params.sc1:
            raise ValueError('invalid interface stretches')
        self.broken = np.zeros(len(self.lengths), dtype=bool)
        self.prev_stretch = np.zeros(len(self.lengths))

    def evaluate(self, concrete_u, steel_u, commit=False):
        concrete_u = np.asarray(concrete_u, dtype=float)
        steel_u = np.asarray(steel_u, dtype=float)
        interp = np.sum(concrete_u[self.indices] * self.weights[:, :, None], axis=1)
        deformed = self.reference + interp - steel_u[self.steel_ids]
        new_lengths = np.linalg.norm(deformed, axis=1)
        if np.any(new_lengths <= 1e-16):
            raise FloatingPointError('collapsed interface bond')
        stretch = (new_lengths - self.lengths) / self.lengths
        prm = self.params
        c = prm.bcc * (1 - (self.lengths / prm.horizon) ** 2) ** 2
        scalar = np.where(stretch <= prm.s01, c * stretch, c * prm.s01 * np.clip((prm.sc1 - stretch) / (prm.sc1 - prm.s01), 0, 1))
        newly_broken = stretch > prm.sc1
        effective_broken = self.broken | newly_broken
        scalar = np.where(effective_broken, 0.0, scalar)
        sample_force = (scalar * prm.quadrature_volume)[:, None] * deformed / new_lengths[:, None]
        fsteel = np.zeros_like(steel_u)
        fconcrete = np.zeros_like(concrete_u)
        np.add.at(fsteel, self.steel_ids, sample_force)
        for a in range(8):
            np.add.at(fconcrete, self.indices[:, a], -self.weights[:, a, None] * sample_force)
        if commit:
            self.broken |= newly_broken
            self.prev_stretch[:] = stretch
        return (fconcrete, fsteel, stretch, effective_broken)

    def reset(self):
        self.broken.fill(False)
        self.prev_stretch.fill(0.0)

# Stage 19 notebook cell 96; preserved mathematical operations
@njit(cache=False)
def _fast_pdrod_kernel(uc, us, ref, ref_len, st_ids, indices, weights, broken, prev, bcc, s01, sc1, horizon, quadrature_volume, commit):
    Fcon = np.zeros((uc.shape[0], 3), dtype=np.float64)
    Fsteel = np.zeros((us.shape[0], 3), dtype=np.float64)
    failed_count = 0
    for k in range(ref_len.size):
        sid = st_ids[k]
        dxv = ref[k, 0] - us[sid, 0]
        dyv = ref[k, 1] - us[sid, 1]
        dzv = ref[k, 2] - us[sid, 2]
        for a in range(8):
            idx = indices[k, a]
            wt = weights[k, a]
            dxv += wt * uc[idx, 0]
            dyv += wt * uc[idx, 1]
            dzv += wt * uc[idx, 2]
        length = np.sqrt(dxv * dxv + dyv * dyv + dzv * dzv)
        if length <= 1e-15:
            raise ValueError('collapsed PDROD interface bond')
        stretch = (length - ref_len[k]) / ref_len[k]
        failed = broken[k] or stretch > sc1
        if failed:
            scalar = 0.0
            failed_count += 1
        else:
            c = bcc * (1.0 - (ref_len[k] / horizon) ** 2) ** 2
            if stretch <= s01:
                scalar = c * stretch
            else:
                scalar = c * s01 * (sc1 - stretch) / (sc1 - s01)
        amp = scalar * quadrature_volume / length
        fx = amp * dxv
        fy = amp * dyv
        fz = amp * dzv
        Fsteel[sid, 0] += fx
        Fsteel[sid, 1] += fy
        Fsteel[sid, 2] += fz
        for a in range(8):
            idx = indices[k, a]
            wt = weights[k, a]
            Fcon[idx, 0] -= wt * fx
            Fcon[idx, 1] -= wt * fy
            Fcon[idx, 2] -= wt * fz
        if commit:
            broken[k] = failed
            prev[k] = stretch
    return (Fcon, Fsteel, failed_count)

class FastInterface:

    def __init__(self, original):
        self.original = original
        self.broken = np.zeros_like(original.broken)
        self.prev_stretch = np.zeros_like(original.prev_stretch)

    def reset(self):
        self.broken[:] = False
        self.prev_stretch[:] = 0.0

    def evaluate(self, uc, us, commit=False):
        it = self.original
        p0 = it.params
        return _fast_pdrod_kernel(uc, us, it.reference, it.lengths, it.steel_ids, it.indices, it.weights, self.broken, self.prev_stretch, p0.bcc, p0.s01, p0.sc1, p0.horizon, p0.quadrature_volume, commit)

# Stage 19 notebook cell 118; preserved mathematical operations
@njit(cache=False)
def _stage18_pdrod_conservative(uc, us, ref, ref_len, sids, ids, w, broken, prev, bcc, s01, sc1, horizon, vol, steel_cv, commit):
    fc = np.zeros((uc.shape[0], 3), dtype=np.float64)
    fs = np.zeros((us.shape[0], 3), dtype=np.float64)
    nf = 0
    for k in range(ref_len.size):
        sid = sids[k]
        ax = ref[k, 0] - us[sid, 0]
        ay = ref[k, 1] - us[sid, 1]
        az = ref[k, 2] - us[sid, 2]
        for a in range(8):
            j = ids[k, a]
            ww = w[k, a]
            ax += ww * uc[j, 0]
            ay += ww * uc[j, 1]
            az += ww * uc[j, 2]
        curr = np.sqrt(ax * ax + ay * ay + az * az)
        if curr < 1e-15:
            raise ValueError('collapsed PDROD bond')
        st = (curr - ref_len[k]) / ref_len[k]
        fail = broken[k] or st > sc1
        if fail:
            scalar = 0.0
            nf += 1
        else:
            c = bcc * (1.0 - (ref_len[k] / horizon) ** 2) ** 2
            if st <= s01:
                scalar = c * st
            else:
                scalar = c * s01 * (sc1 - st) / (sc1 - s01)
        amp = scalar * vol * steel_cv[sid] / curr
        fx = amp * ax
        fy = amp * ay
        fz = amp * az
        fs[sid, 0] += fx
        fs[sid, 1] += fy
        fs[sid, 2] += fz
        for a in range(8):
            j = ids[k, a]
            ww = w[k, a]
            fc[j, 0] -= ww * fx
            fc[j, 1] -= ww * fy
            fc[j, 2] -= ww * fz
        if commit:
            broken[k] = fail
            prev[k] = st
    return (fc, fs, nf)

def conservative_interface(uc, us, commit=False):
    obj = fast_interface.original
    pp = obj.params
    return _stage18_pdrod_conservative(uc, us, obj.reference, obj.lengths, obj.steel_ids, obj.indices, obj.weights, fast_interface.broken, fast_interface.prev_stretch, pp.bcc, pp.s01, pp.sc1, pp.horizon, pp.quadrature_volume, STEEL_CONTROL_VOLUME, commit)

# Stage 19 notebook cell 104; preserved mathematical operations
@dataclass(frozen=True)
class AppendixCParameters:
    bcc: float
    s01: float
    sc1: float
    delta: float
    vol: float

class AppendixCPublishedKernel:
    """Vectorized form of Appendix C using explicit 0-based steel/concrete indexing."""

    def __init__(self, steel_xyz, concrete_xyz, steel_id, concrete8_id, weights8, virtual_xyz, params: AppendixCParameters):
        self.xs = np.asarray(steel_xyz, float)
        self.xc = np.asarray(concrete_xyz, float)
        self.st = np.asarray(steel_id, np.int64)
        self.idx = np.asarray(concrete8_id, np.int64)
        self.w = np.asarray(weights8, float)
        self.xv = np.asarray(virtual_xyz, float)
        self.p = params
        if self.idx.shape != (len(self.st), 8) or self.w.shape != self.idx.shape:
            raise ValueError('Each coupling entry must have exactly 8 concrete indices and weights')
        if self.xv.shape != (len(self.st), 3):
            raise ValueError('virtual_xyz must have shape (K,3)')
        if np.any(self.st < 0) or np.any(self.st >= len(self.xs)) or np.any(self.idx < 0) or np.any(self.idx >= len(self.xc)):
            raise IndexError('Input connectivity is outside the respective steel or concrete array')
        if np.max(np.abs(self.w.sum(axis=1) - 1)) > 1e-10:
            raise ValueError('Shape weights must sum to one')
        self.reference = self.xv - self.xs[self.st]
        self.idist = np.linalg.norm(self.reference, axis=1)
        if np.any(self.idist <= 0):
            raise ValueError('PDROD bonds cannot have zero initial length')
        if not (params.bcc > 0 and params.delta > 0 and (params.vol > 0) and (0 < params.s01 < params.sc1)):
            raise ValueError('Invalid PDROD constants')
        if np.any(self.idist > params.delta):
            raise ValueError('Neighbors must be within delta')
        self.failfail = np.ones(len(self.st), dtype=np.int8)

    def evaluate(self, uc, us, old_uc=None, old_us=None, dt=None, mode='literal', commit=False):
        if mode not in ('literal', 'irreversible_corrected'):
            raise ValueError('mode must be literal or irreversible_corrected')
        uc = np.asarray(uc, float)
        us = np.asarray(us, float)
        if uc.shape != self.xc.shape or us.shape != self.xs.shape:
            raise ValueError('Displacements must match geometry shape')
        dispdisp = np.einsum('ka,kad->kd', self.w, uc[self.idx])
        direction = self.reference + dispdisp - us[self.st]
        nlength = np.linalg.norm(direction, axis=1)
        if np.any(nlength < 1e-15):
            raise FloatingPointError('Current interface length is zero')
        sssss = (nlength - self.idist) / self.idist
        if old_uc is not None and old_us is not None and (dt is not None):
            if dt <= 0:
                raise ValueError('dt must be positive')
            old_interp = np.einsum('ka,kad->kd', self.w, np.asarray(old_uc)[self.idx])
            oldnlength = np.linalg.norm(self.reference + old_interp - np.asarray(old_us)[self.st], axis=1)
            s1 = np.abs((nlength - oldnlength) / self.idist / dt)
        else:
            oldnlength = np.full_like(nlength, np.nan)
            s1 = np.full_like(nlength, np.nan)
        bc = self.p.bcc * (1 - (self.idist / self.p.delta) ** 2) ** 2
        scalar = np.where(sssss <= self.p.s01, bc * sssss, np.where(sssss <= self.p.sc1, bc * self.p.s01 * (self.p.sc1 - sssss) / (self.p.sc1 - self.p.s01), 0.0))
        new_failed = sssss > self.p.sc1
        if mode == 'irreversible_corrected':
            scalar = np.where((self.failfail == 0) | new_failed, 0.0, scalar)
        dforce = (scalar * self.p.vol)[:, None] * (direction / nlength[:, None])
        fs = np.zeros_like(us)
        fc = np.zeros_like(uc)
        np.add.at(fs, self.st, dforce)
        for a in range(8):
            np.add.at(fc, self.idx[:, a], -self.w[:, a, None] * dforce)
        if commit:
            self.failfail[new_failed] = 0
        dmgpar1 = float(np.sum(self.failfail * self.p.vol))
        dmgpar2 = float(len(self.failfail) * self.p.vol)
        return dict(concrete_force=fc, steel_force=fs, stretch=sssss, stretch_rate=s1, oldnlength=oldnlength, current_length=nlength, dforce=dforce, new_failed=new_failed, failfail=self.failfail.copy(), sample_damage=1 - dmgpar1 / dmgpar2 if dmgpar2 else np.nan)

    def reset(self):
        self.failfail[:] = 1
