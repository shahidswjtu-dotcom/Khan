"""Exact instantaneous potential of the implemented bilinear envelope.

Rupture release is the potential removed when a bond is deleted. Rate-induced
potential changes are logged separately. These are algorithmic energies, not a
validated unloading/dissipation model or proof of discrete fracture Gf.
"""
import numpy as np
from numba import njit

@njit
def primitive(s,a,b):
    if s<=a:return .5*s*s
    t=min(s,b)-a
    return .5*a*a+a*t-.5*a*t*t/(b-a)

@njit
def concrete_energy(X,U,I,J,r0,weight,previous,broken,dt,E,ft,Gf,delta,fc,V):
    stored=0.;old_current=0.;release=0.
    exponent=1/(10+6*fc/1e7);theta=10**(6*exponent-2)
    a0=ft/E;b0=5*Gf/(6*E*delta*a0)
    for k in range(len(I)):
        if broken[k]:continue
        i=I[k];j=J[k]
        ax=X[j,0]+U[j,0]-X[i,0]-U[i,0]
        ay=X[j,1]+U[j,1]-X[i,1]-U[i,1]
        az=X[j,2]+U[j,2]-X[i,2]-U[i,2]
        length=np.sqrt(ax*ax+ay*ay+az*az);s=(length-r0[k])/r0[k]
        rate=(s-previous[k])/dt
        if rate<=3e-6:t=1.;e=1.
        else:
            ratio=rate/3e-6;e=ratio**.016
            t=ratio**exponent if rate<=1 else theta*ratio**(1/3)
        a=a0*t/e;b=b0*t/e
        coeff=72*E/(np.pi*delta**4)*(1-(r0[k]/delta)**2)**2*e*r0[k]*V**2*weight[k]
        old_current+=coeff*primitive(previous[k],a,b)
        en=coeff*primitive(s,a,b)
        if s>=b:release+=en
        else:stored+=en
    return stored,old_current,release

@njit
def interface_energy(uc,us,ref,r0,sids,ids,w,previous,broken,bcc,a,b,horizon,V,steel_cv):
    stored=0.;release=0.
    for k in range(len(r0)):
        if broken[k]:continue
        d=ref[k]-us[sids[k]]
        for j in range(8):d=d+w[k,j]*uc[ids[k,j]]
        s=(np.sqrt(np.dot(d,d))-r0[k])/r0[k]
        coeff=bcc*(1-(r0[k]/horizon)**2)**2*V*steel_cv[sids[k]]*r0[k]
        en=coeff*primitive(s,a,b)
        if s>b:release+=en
        else:stored+=en
    return stored,release

def rod_energy(m,us):
    unit=(m.steel_xyz[m.rod_pairs[:,1]]-m.steel_xyz[m.rod_pairs[:,0]])/m.steel_lengths[:,None]
    ext=np.sum((us[m.rod_pairs[:,1]]-us[m.rod_pairs[:,0]])*unit,axis=1)
    return float(np.sum(.5*m.steel_E*m.steel_A/m.steel_lengths*ext**2))
