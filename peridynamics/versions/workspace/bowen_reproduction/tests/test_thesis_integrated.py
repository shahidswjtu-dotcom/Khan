import numpy as np
from thesis_pd.model import create,evaluate
from thesis_pd.integrator import Simulation

def test_integrated_surface_and_self_contact():
 m=create(corrected=True,family=True)
 from bowen_pd.geometry import PAIR_WEIGHT
 assert np.max(abs(m.PAIR_WEIGHT-PAIR_WEIGHT))>0
 np.testing.assert_array_equal(m.PAIR_WEIGHT,m.scope['PAIR_WEIGHT'])
 assert len(m.interface.lengths)==27572
 u=np.zeros_like(m.xyz);us=np.zeros_like(m.steel_xyz)
 r=evaluate(m,u,us,m.hammer_z_initial,m.p.dt)
 assert r['repulsion_pairs']==0 and r['repulsion_potential_J']==0
 assert np.max(abs(r['fc']))<1e-7
 # Move an interior point toward its neighbour; short-range force must be active even without broken bonds.
 k=17000;u[k,0]=-.0015
 r=evaluate(m,u,us,m.hammer_z_initial,m.p.dt)
 assert r['repulsion_pairs']>0 and r['repulsion_potential_J']>0
 np.testing.assert_allclose(r['internal_c'].sum(0)+r['internal_s'].sum(0),0,atol=2e-7)

def test_projection_restart_including_contact_observation(tmp_path):
 a=Simulation('family_projection');b=Simulation('family_projection')
 for _ in range(350):a.advance()
 a.checkpoint(tmp_path/'c.npz');b.restore(tmp_path/'c.npz')
 assert a.record()==b.record()
 for _ in range(3):a.advance();b.advance()
 np.testing.assert_array_equal(a.u,b.u);np.testing.assert_array_equal(a.v,b.v)
 assert a.record()==b.record()
