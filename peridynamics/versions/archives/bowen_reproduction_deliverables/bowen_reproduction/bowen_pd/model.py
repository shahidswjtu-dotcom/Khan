"""Explicit legacy reconstruction; every demo coefficient remains labelled."""
from types import SimpleNamespace
from . import geometry as g, concrete as c, interface as interface, rod, contact
import numpy as np

def create(corrected=False, family=False):
    # geometry module is immutable baseline. Copy corrected arrays, never mutate it.
    ns=SimpleNamespace(**{k:v for k,v in vars(g).items() if not k.startswith('__')})
    if corrected:
        ns.steel_xyz=g.steel_xyz.copy(); ns.steel_xyz[:,2]=.020
        ns.sample_pos=g.sample_pos.copy(); ns.sample_pos[:,2]=.020
        pairs=[g.concrete_grid_eight(q) for q in ns.sample_pos]
        ns.sample_idx=np.array([q[0] for q in pairs]);ns.sample_w=np.array([q[1] for q in pairs])
    if family:
        # Console2.f90 lines240--278: integer offsets over a 7^3 lattice,
        # nonzero vectors inside horizon, then 8-node interpolation.
        # Geometry/coefficients remain Chapter3 reconstruction, not author code.
        samples=[];sid=[];indices=[];weights=[]
        from itertools import product
        lower=g.xyz.min(0);upper=g.xyz.max(0)
        for i,pos in enumerate(ns.steel_xyz):
            for offset in product(range(-3,4),repeat=3):
                offset=np.asarray(offset)*g.dx
                if np.linalg.norm(offset)<1e-12 or np.linalg.norm(offset)>g.delta*(1+1e-12):continue
                pt=pos+offset
                if np.any(pt<lower-1e-12) or np.any(pt>=upper-1e-12):continue
                idx,w=g.concrete_grid_eight(pt)
                samples.append(pt);sid.append(i);indices.append(idx);weights.append(w)
        ns.sample_pos=np.array(samples);ns.sample_steel=np.array(sid,dtype=np.int32)
        ns.sample_idx=np.array(indices,dtype=np.int32);ns.sample_w=np.array(weights)
    # extracted functions use module globals; bind independent copies by function cloning
    import types
    scope=dict(vars(ns))
    scope.update(CONCRETE_NODE_MASS=g.CONCRETE_NODE_MASS,STEEL_CONTROL_VOLUME=g.STEEL_MASS/g.p.rho_steel)
    scope['rod_forces']=types.FunctionType(rod.rod_forces.__code__,scope)
    scope['sphere_contact']=types.FunctionType(contact.sphere_contact.__code__,scope,argdefs=contact.sphere_contact.__defaults__)
    params=interface.AppendixCInterfaceParameters(72*g.p.E_concrete/(np.pi*g.delta**4),g.p.ft_concrete/g.p.E_concrete,5*g.p.Gf_concrete/(6*g.p.E_concrete*g.delta*(g.p.ft_concrete/g.p.E_concrete)),g.delta,g.V)
    obj=interface.AppendixCPDROD(ns.steel_xyz,ns.sample_pos,ns.sample_idx,ns.sample_w,ns.sample_steel,params)
    ns.scope=scope;ns.params=params;ns.interface=obj
    ns.broken=np.zeros(len(g.I),bool);ns.prev=np.zeros(len(g.I))
    ns.interface_broken=np.zeros(len(ns.sample_steel),bool);ns.interface_prev=np.zeros(len(ns.sample_steel))
    return ns

def evaluate(m,uc,us,zh,dt,commit=False,penalty=True):
    g=m.scope;p=g['p'];it=m.interface;pp=m.params
    from .energies import concrete_energy,interface_energy,rod_energy
    ce,old_ce,cr=concrete_energy(g['xyz'],uc,g['I'],g['J'],g['LEN'],g['PAIR_WEIGHT'],m.prev,m.broken,dt,p.E_concrete,p.ft_concrete,p.Gf_concrete,g['delta'],p.fc_concrete,g['V'])
    ie,ir=interface_energy(uc,us,it.reference,it.lengths,it.steel_ids,it.indices,it.weights,m.interface_prev,m.interface_broken,pp.bcc,pp.s01,pp.sc1,pp.horizon,pp.quadrature_volume,g['STEEL_CONTROL_VOLUME'])
    fc,nb=c._fast_concrete_kernel(g['xyz'],uc,g['I'],g['J'],g['LEN'],g['PAIR_WEIGHT'],m.prev,m.broken,dt,p.E_concrete,p.ft_concrete,p.Gf_concrete,g['delta'],p.fc_concrete,g['V']**2,commit)
    fs=g['rod_forces'](us)
    ic,ist,ni=interface._stage18_pdrod_conservative(uc,us,it.reference,it.lengths,it.steel_ids,it.indices,it.weights,m.interface_broken,m.interface_prev,pp.bcc,pp.s01,pp.sc1,pp.horizon,pp.quadrature_volume,g['STEEL_CONTROL_VOLUME'],commit)
    cont,hf,ep,nc,pen=g['sphere_contact'](g['xyz']+uc,zh) if penalty else (np.zeros_like(uc),np.zeros(3),0.,0,0.)
    return dict(concrete_envelope_potential_J=ce,old_potential_current_rate_J=old_ce,concrete_rupture_release_J=cr,interface_envelope_potential_J=ie,interface_rupture_release_J=ir,rod_potential_J=rod_energy(m,us),fc=fc+ic+cont,fs=fs+ist,hf=hf,ep=ep,nc=nc,pen=pen,nb=nb,ni=ni,internal_c=fc+ic,internal_s=fs+ist,interface_residual=float(np.linalg.norm(ic.sum(0)+ist.sum(0))))
