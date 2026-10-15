"""Frozen-parameter energy accounting, separate from force integration."""
import numpy as np
from numba import njit
from .operators import law_energy,dynamic_thresholds

@njit(cache=True)
def threshold_ledger(x,u,previous,ends,r0,coeff,alive,a_stored,b_stored,a0,b0,dt,rate,irreversible):
    exchange=0.; rupture=0.
    for k in range(len(ends)):
        i,j=ends[k]
        dx=x[j,0]+u[j,0]-x[i,0]-u[i,0]
        dy=x[j,1]+u[j,1]-x[i,1]-u[i,1]
        dz=x[j,2]+u[j,2]-x[i,2]-u[i,2]
        r=np.sqrt(dx*dx+dy*dy+dz*dz); s=r/r0[k]-1
        a,b=a0,b0
        if rate:
            ox=x[j,0]+previous[j,0]-x[i,0]-previous[i,0]
            oy=x[j,1]+previous[j,1]-x[i,1]-previous[i,1]
            oz=x[j,2]+previous[j,2]-x[i,2]-previous[i,2]
            oldr=np.sqrt(ox*ox+oy*oy+oz*oz)
            a,b=dynamic_thresholds((r-oldr)/(r0[k]*dt),a0,b0)
        if alive[k] or not irreversible:
            oldpsi=law_energy(s,a_stored[k],b_stored[k])[1]
            psi=law_energy(s,a,b)[1]
            exchange+=coeff[k]*r0[k]*(psi-oldpsi)
            if irreversible and s>b:
                rupture+=coeff[k]*r0[k]*psi
        a_stored[k]=a; b_stored[k]=b
    return exchange,rupture

@njit(cache=True)
def concrete_potential(x,u,ends,r0,coeff,alive,a,b,irreversible):
    energy=0.
    for k in range(len(ends)):
        if irreversible and not alive[k]: continue
        i,j=ends[k]; d=x[j]+u[j]-x[i]-u[i]
        s=np.sqrt(np.dot(d,d))/r0[k]-1
        energy+=coeff[k]*r0[k]*law_energy(s,a[k],b[k])[1]
    return energy

@njit(cache=True)
def interface_rupture(x,u,sx,su,owners,stencils,weights,r0,coeff,alive,a,b,irreversible):
    loss=0.
    if not irreversible: return loss
    for k in range(len(owners)):
        if not alive[k]: continue
        virtual=np.zeros(3)
        for q in range(8):
            j=stencils[k,q]; virtual+=weights[k,q]*(x[j]+u[j])
        d=virtual-sx[owners[k]]-su[owners[k]]
        s=np.sqrt(np.dot(d,d))/r0[k]-1
        if s>b: loss+=coeff[k]*r0[k]*law_energy(s,a,b)[1]
    return loss

@njit(cache=True)
def rod_potential(x,u,ends,r0,area,young,plastic):
    total=0.
    for k in range(len(ends)):
        i,j=ends[k]; d=x[j]+u[j]-x[i]-u[i]
        elastic=np.sqrt(np.dot(d,d))/r0[k]-1-plastic[k]
        total+=.5*young*elastic**2*area*r0[k]
    return total

@njit(cache=True)
def interface_potential(x,u,sx,su,owners,stencils,weights,r0,coeff,alive,a,b,irreversible):
    total=0.
    for k in range(len(owners)):
        if irreversible and not alive[k]: continue
        virtual=np.zeros(3)
        for q in range(8):
            j=stencils[k,q]; virtual+=weights[k,q]*(x[j]+u[j])
        d=virtual-sx[owners[k]]-su[owners[k]]
        s=np.sqrt(np.dot(d,d))/r0[k]-1
        total+=coeff[k]*r0[k]*law_energy(s,a,b)[1]
    return total
