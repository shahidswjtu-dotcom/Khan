import numpy as np
import pytest
from source_methods import surface_energy_correction,short_range_repulsion
from bowen_pd.geometry import xyz,I,J,p,delta,dx

def test_directional_energy_patch_correction():
 node,pair,energy,reference=surface_energy_correction(xyz,I,J,p.E_concrete,p.nu_concrete,delta,dx)
 np.testing.assert_allclose(node*energy,reference,rtol=1e-13)
 assert np.all(pair>0) and np.all(np.isfinite(pair))
 # All permutations of endpoint IDs give exactly the same symmetric correction.
 _,reverse,_,_=surface_energy_correction(xyz,J,I,p.E_concrete,p.nu_concrete,delta,dx)
 np.testing.assert_allclose(pair,reverse,rtol=1e-12)

def test_repulsion_force_potential_momentum():
 x=np.array([[0.,0,0],[.005,0,0]])
 u=np.zeros_like(x);u[1,0]=-.001
 force,energy,count=short_range_repulsion(x,u,p.E_concrete,delta,dx)
 assert count==1 and force[0,0]<0 and force[1,0]>0 and energy>0
 np.testing.assert_allclose(force.sum(0),0)
 eps=1e-8;up=u.copy();um=u.copy();up[1,0]+=eps;um[1,0]-=eps
 derivative=(short_range_repulsion(x,up,p.E_concrete,delta,dx)[1]-short_range_repulsion(x,um,p.E_concrete,delta,dx)[1])/(2*eps)
 assert force[1,0]==pytest.approx(-derivative,rel=1e-8)
 np.testing.assert_allclose(np.cross(x+u,force).sum(0),0)

def test_new_contact_outside_original_horizon():
 x=np.array([[0.,0,0],[.03,0,0]])
 u=np.zeros_like(x);u[1,0]=-.026
 force,energy,count=short_range_repulsion(x,u,p.E_concrete,delta,dx)
 assert count==1 and force[1,0]>0
 # New pair contact must not depend on membership or survival of original PD bonds.

def test_unloaded_reference_no_repulsion():
 x=np.array([[0.,0,0],[dx,0,0],[0,dx,0]])
 force,energy,count=short_range_repulsion(x,np.zeros_like(x),p.E_concrete,delta,dx)
 assert count==0 and energy==0 and np.max(abs(force))==0
