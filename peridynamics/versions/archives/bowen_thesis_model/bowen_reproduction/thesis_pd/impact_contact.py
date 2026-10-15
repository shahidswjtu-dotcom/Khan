"""Thesis Eqs 2-45--47, explicit projection, axial rigid hammer.

Force is impulse/dt in N. Hammer position uses pre-contact predicted motion;
recoil changes next step's velocity. This split and axial guide are assumptions.
"""
import numpy as np

def project_contact(xref,u,old_u,v,zh,vh,dt,mass,hammer_mass,xy,radius):
    center=np.array([*xy,zh]);d=xref+u-center
    dist=np.linalg.norm(d,axis=1);active=dist<radius
    if np.any(dist[active]<=1e-15): raise FloatingPointError('Particle at hammer centre')
    correction=np.zeros_like(u)
    correction[active]=d[active]*((radius/dist[active])-1)[:,None]
    before=v[active].copy()
    u[active]+=correction[active]
    v[active]=(u[active]-old_u[active])/dt
    impulse=np.sum(mass*(v[active]-before),axis=0)
    old_vh=vh;vh-=impulse[2]/hammer_mass
    delta_ke=.5*mass*np.sum(v[active]**2-before**2)+.5*hammer_mass*(vh**2-old_vh**2)
    return vh,-impulse/dt,int(active.sum()),float(np.max(radius-dist[active],initial=0)),float(delta_ke)
