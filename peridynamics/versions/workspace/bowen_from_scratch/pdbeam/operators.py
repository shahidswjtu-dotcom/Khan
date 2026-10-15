"""Independent compiled operators. All returned forces are physical newtons."""
import numpy as np
from numba import njit

@njit(cache=True)
def dynamic_thresholds(rate, onset, terminal):
    # Eqs. 2-13--15, fc=60 MPa; signed tensile bond rate.
    if rate <= 3e-6:
        return onset, terminal
    alpha=1/(10+6*60/10)
    dif=(rate/3e-6)**alpha if rate<=1 else 10**(6*alpha-2)*(rate/3e-6)**(1/3)
    modulus=(rate/3e-6)**.016
    return onset*dif/modulus, terminal*dif/modulus

@njit(cache=True)
def law_energy(s,a,b):
    if s<=a:
        return s,.5*s*s
    if s<=b:
        q=s-a
        return a*(b-s)/(b-a), .5*a*a+a*q-.5*a*q*q/(b-a)
    return 0.,.5*a*b

@njit(cache=True)
def concrete(x,u,previous,ends,r0,coefficient,alive,onset,terminal,dt,irreversible,rate_enabled):
    f=np.zeros_like(x); potential=0.
    for k in range(len(ends)):
        if irreversible and not alive[k]: continue
        i,j=ends[k]
        dx=x[j,0]+u[j,0]-x[i,0]-u[i,0]
        dy=x[j,1]+u[j,1]-x[i,1]-u[i,1]
        dz=x[j,2]+u[j,2]-x[i,2]-u[i,2]
        r=np.sqrt(dx*dx+dy*dy+dz*dz)
        if r<1e-15: continue
        s=r/r0[k]-1
        a,b=onset,terminal
        if rate_enabled:
            ox=x[j,0]+previous[j,0]-x[i,0]-previous[i,0]
            oy=x[j,1]+previous[j,1]-x[i,1]-previous[i,1]
            oz=x[j,2]+previous[j,2]-x[i,2]-previous[i,2]
            oldr=np.sqrt(ox*ox+oy*oy+oz*oz)
            a,b=dynamic_thresholds((r-oldr)/(r0[k]*dt),a,b)
        value,psi=law_energy(s,a,b)
        alive[k]=alive[k] and s<=b
        if irreversible and not alive[k]: continue
        scale=coefficient[k]*value/r
        f[i,0]+=scale*dx; f[j,0]-=scale*dx
        f[i,1]+=scale*dy; f[j,1]-=scale*dy
        f[i,2]+=scale*dz; f[j,2]-=scale*dz
        potential+=coefficient[k]*r0[k]*psi
    return f,potential

@njit(cache=True)
def coupling(x,u,sx,su,owners,stencils,weights,r0,coefficients,alive,onset,terminal,irreversible):
    fc=np.zeros_like(x); fs=np.zeros_like(sx); potential=0.
    for k in range(len(owners)):
        if irreversible and not alive[k]: continue
        a=owners[k]; virtual=np.zeros(3)
        for q in range(8):
            j=stencils[k,q]
            virtual+=weights[k,q]*(x[j]+u[j])
        d=virtual-sx[a]-su[a]; r=np.sqrt(np.dot(d,d))
        if r<1e-15:
            continue
        value,psi=law_energy(r/r0[k]-1,onset,terminal)
        alive[k]=alive[k] and r/r0[k]-1<=terminal
        if irreversible and not alive[k]:
            continue
        pair=coefficients[k]*value*d/r
        fs[a]+=pair
        for q in range(8):
            fc[stencils[k,q]]-=weights[k,q]*pair
        potential+=coefficients[k]*r0[k]*psi
    return fc,fs,potential

@njit(cache=True)
def rods(x,u,ends,r0,area,young,yield_stress,plastic):
    f=np.zeros_like(x); potential=0.; dissipation=0.
    for k in range(len(ends)):
        i,j=ends[k]; d=x[j]+u[j]-x[i]-u[i]; r=np.sqrt(np.dot(d,d))
        if r<1e-15:
            continue
        strain=r/r0[k]-1
        trial=young*(strain-plastic[k])
        stress=min(yield_stress,max(-yield_stress,trial))
        increment=(trial-stress)/young
        plastic[k]+=increment
        dissipation+=yield_stress*abs(increment)*area*r0[k]
        pair=area*stress*d/r
        f[i]+=pair; f[j]-=pair
        potential+=.5*stress**2/young*area*r0[k]
    return f,potential,dissipation

@njit(cache=True)
def repulsion(reference,current,ends,spacing,horizon,young):
    f=np.zeros_like(current); potential=0.
    stiffness=1.5*12*young/(np.pi*horizon**4)
    coefficient=stiffness/horizon*spacing**6
    for k in range(len(ends)):
        i,j=ends[k]; d=current[j]-current[i]; r=np.sqrt(np.dot(d,d))
        old=reference[j]-reference[i]
        exclusion=min(.9*np.sqrt(np.dot(old,old)),1.35*spacing)
        if 1e-15<r<exclusion:
            pair=coefficient*(r-exclusion)*d/r
            f[i]+=pair; f[j]-=pair
            potential+=.5*coefficient*(r-exclusion)**2
    return f,potential
