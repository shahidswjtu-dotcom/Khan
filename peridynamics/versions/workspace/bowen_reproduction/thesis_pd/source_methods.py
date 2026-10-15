"""Source-supported missing methods, independently verifiable before integration.

No material fitting. These functions are outside the verified engine package:
its checkpoint-compatible source hash and legacy runs remain unchanged.
"""
import numpy as np
from scipy.spatial import cKDTree

def surface_energy_correction(xyz,I,J,E,nu,delta,spacing,strain=1e-3,weighted=True):
    """Eq2-40 and Console2 lines459--529, 3 reference strain patch tests.

    Return per-node directional correction and symmetric per-bond ellipsoidal
    interpolation as used in Console2 lines628--653. Weighted kernel is Eq2-12;
    constant kernel reproduces the supplied Fortran's local micro-modulus.
    """
    if not (-1<nu<.5) or E<=0 or spacing<=0 or delta<=0 or strain<=0:
        raise ValueError('Invalid material/geometric inputs')
    d=xyz[J]-xyz[I];r=np.linalg.norm(d,axis=1)
    if np.any(r<=0) or np.any(r>delta*(1+1e-10)):raise ValueError('Invalid horizon pairs')
    volume=spacing**3;vc=np.clip((delta+spacing/2-r)/spacing,0,1)
    c=(72*E/(np.pi*delta**4)*(1-(r/delta)**2)**2 if weighted else np.full_like(r,12*E/(np.pi*delta**4)))
    lam=E*nu/((1+nu)*(1-2*nu));mu=E/(2*(1+nu))
    reference=.5*(lam+2*mu)*strain**2
    correction=np.empty((len(xyz),3));energies=np.empty_like(correction)
    for axis in range(3):
        dd=d.copy();dd[:,axis]*=1+strain
        stretch=(np.linalg.norm(dd,axis=1)-r)/r
        # A point receives half of the physical unique pair potential / volume.
        per_endpoint=.25*c*r*stretch**2*volume*vc
        W=np.bincount(I,weights=per_endpoint,minlength=len(xyz))+np.bincount(J,weights=per_endpoint,minlength=len(xyz))
        if np.any(W<=0):raise ValueError('No strain energy in a particle family')
        energies[:,axis]=W;correction[:,axis]=reference/W
    midpoint=.5*(correction[I]+correction[J]);unit=d/r[:,None]
    factor=1/np.sqrt(np.sum((unit/midpoint)**2,axis=1))
    return correction,factor,energies,reference

def short_range_repulsion(xref,u,E,delta,spacing,search='all_current',micromodulus='constant'):
    """Eqs2-48--49, current-neighbour search independent of broken flags.

    N=law*ViVj. Constant c=12E/(pi delta^4) is explicit for all-current pairs.
    A reference-distance weighted c for new pairs outside delta is undefined
    by the supplied text: reject that combination rather than extrapolate.
    """
    if min(E,delta,spacing)<=0:raise ValueError('Positive material/geometric values required')
    x=np.asarray(xref)+np.asarray(u)
    if search!='all_current' or micromodulus!='constant':
        raise ValueError('Only explicit all-current/constant convention supported')
    pairs=cKDTree(x).query_pairs(1.35*spacing,output_type='ndarray')
    F=np.zeros_like(x,dtype=float)
    if len(pairs)==0:return F,0.,0
    I,J=pairs.T;d=x[J]-x[I];r=np.linalg.norm(d,axis=1)
    if np.any(r<=1e-15):raise FloatingPointError('Coincident particles: repulsion direction undefined')
    reference=np.linalg.norm(xref[J]-xref[I],axis=1)
    cutoff=np.minimum(.9*reference,1.35*spacing);penetration=np.maximum(cutoff-r,0)
    cs=1.5*12*E/(np.pi*delta**4);k=cs/delta*spacing**6
    force=(-k*penetration/r)[:,None]*d
    np.add.at(F,I,force);np.add.at(F,J,-force)
    return F,float(.5*k*np.sum(penetration**2)),int(np.sum(penetration>0))
