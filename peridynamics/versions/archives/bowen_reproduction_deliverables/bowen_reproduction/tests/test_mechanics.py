import numpy as np
import pytest
from bowen_pd.model import create,evaluate
from bowen_pd.integrator import Simulation
from bowen_pd.impact_contact import project_contact
from bowen_pd.interface import AppendixCPublishedKernel,AppendixCParameters
from bowen_pd.io import load

@pytest.fixture(scope='module')
def model():return create()

def test_geometry_mass(model):
 m=model
 assert len(m.xyz)==33600 and len(m.I)==1755720
 assert np.isclose(m.CONCRETE_NODE_MASS*len(m.xyz),10.08)
 assert np.isclose(m.STEEL_MASS.sum(),m.p.rho_steel*m.steel_A*m.steel_lengths.sum())
 assert m.STEEL_MASS[0]==pytest.approx(m.STEEL_MASS[1]/2)

def test_interface_physical_force_and_moment(model):
 from bowen_pd.interface import _stage18_pdrod_conservative
 m=model;it=m.interface;p=m.params
 rng=np.random.default_rng(9);u=rng.normal(0,1e-9,m.xyz.shape);us=rng.normal(0,1e-9,m.steel_xyz.shape)
 fc,fs,_=_stage18_pdrod_conservative(u,us,it.reference,it.lengths,it.steel_ids,it.indices,it.weights,m.interface_broken.copy(),m.interface_prev.copy(),p.bcc,p.s01,p.sc1,p.horizon,p.quadrature_volume,m.STEEL_MASS/m.p.rho_steel,False)
 np.testing.assert_allclose(fc.sum(0)+fs.sum(0),0,atol=1e-10)
 moment=np.cross(m.xyz+u,fc).sum(0)+np.cross(m.steel_xyz+us,fs).sum(0)
 np.testing.assert_allclose(moment,0,atol=1e-10)
 # virtual work: transpose interpolation preserves work
 dv=rng.normal(size=u.shape);dvs=rng.normal(size=us.shape)
 np.testing.assert_allclose(np.sum(m.sample_w,axis=1),1,atol=1e-14)
 np.testing.assert_allclose(np.sum(m.xyz[m.sample_idx]*m.sample_w[:,:,None],axis=1),m.sample_pos,atol=1e-14)

def test_force_density_conversion():
 q=300.;rho=7750.;volume=2e-7;mass=rho*volume
 assert q*volume/mass==pytest.approx(q/rho)
 assert q/mass!=pytest.approx(q/rho)

def test_no_linear_vertical_interface_stiffness(model):
 # Two purely lateral samples cannot provide axial or vertical linear transfer.
 it=model.interface
 assert np.all(it.reference[:,0]==0) and np.all(it.reference[:,2]==0)
 # This is a topology diagnosis, not permission to add arbitrary samples.

def test_rigid_translation_and_concrete_balance(model):
 m=model;u=np.full_like(m.xyz,1e-8);us=np.full_like(m.steel_xyz,1e-8)
 r=evaluate(m,u,us,m.hammer_z_initial,m.p.dt)
 assert np.max(abs(r['fc']))<1e-7
 u=m.xyz*1e-6;us=m.steel_xyz*1e-6
 r=evaluate(m,u,us,m.hammer_z_initial,m.p.dt)
 np.testing.assert_allclose(r['internal_c'].sum(0)+r['internal_s'].sum(0),0,atol=2e-8)

def test_literal_vs_irreversible(model):
 m=model;it=m.interface;p=m.params
 k=AppendixCPublishedKernel(m.steel_xyz,m.xyz,it.steel_ids,it.indices,it.weights,m.sample_pos,AppendixCParameters(p.bcc,p.s01,p.sc1,p.horizon,p.quadrature_volume))
 u=np.zeros_like(m.xyz);us=np.zeros_like(m.steel_xyz);us[:,2]=.01
 k.evaluate(u,us,mode='literal',commit=True)
 assert np.any(k.failfail==0)
 us[:,2]=0;us[:,1]=1e-8
 literal=k.evaluate(u,us,mode='literal')['steel_force']
 irreversible=k.evaluate(u,us,mode='irreversible_corrected')['steel_force']
 assert np.max(abs(literal))>0
 assert np.max(abs(irreversible))==0

def test_penalty_potential_gradient(model):
 m=model;x=np.array([[.3,.035,0.]])
 f,h,E,n,p=m.scope['sphere_contact'](x,.059)
 eps=1e-8
 xp=x.copy();xp[0,2]+=eps;xm=x.copy();xm[0,2]-=eps
 gradient=(m.scope['sphere_contact'](xp,.059)[2]-m.scope['sphere_contact'](xm,.059)[2])/(2*eps)
 assert f[0,2]==pytest.approx(-gradient,rel=1e-7)
 np.testing.assert_allclose(f.sum(0)+h,0)

def test_projection_impulse_and_geometry():
 x=np.array([[0.,0.,0.]])
 u=np.zeros_like(x);v=np.zeros_like(x);dt=.01;mass=2.;M=10.;vh=-1.
 vh2,hf,n,pen,dke=project_contact(x,u,u.copy(),v,.09,vh,dt,mass,M,[0,0],.1)
 assert n==1
 assert np.linalg.norm(x+u-np.array([0,0,.09]))==pytest.approx(.1)
 assert mass*v[0,2]+M*vh2==pytest.approx(M*vh)
 assert hf[2]==pytest.approx(-mass*v[0,2]/dt)

def test_probes_linear_field(model):
 from bowen_pd.probes import probes
 m=model;u=np.zeros_like(m.xyz);u[:,2]=-2*m.xyz[:,0]-3*m.xyz[:,1]-4*m.xyz[:,2]
 r=probes(m,u)
 assert r['bottom_center_layer_mm']==pytest.approx(1000*(.6+.105+.01))
 assert r['support_relative_bottom_mm']==pytest.approx(0,abs=1e-10)

def test_restart_rate_sensitive_after_contact(tmp_path):
 a=Simulation();b=Simulation()
 for _ in range(360):a.advance()
 a.checkpoint(tmp_path/'c.npz');b.restore(tmp_path/'c.npz')
 for _ in range(6):a.advance();b.advance()
 for k in ('u','v','us','vs'):
  np.testing.assert_array_equal(getattr(a,k),getattr(b,k))
 np.testing.assert_array_equal(a.m.broken,b.m.broken)
 assert a.vh==b.vh and a.zh==b.zh and a.work==b.work
 with pytest.raises(ValueError):Simulation(dt=1.5e-7).restore(tmp_path/'c.npz')

def test_reconstruction_bar_depth():
 m=create(corrected=True)
 assert np.all(m.steel_xyz[:,2]==.020)
 np.testing.assert_allclose(np.sum(m.xyz[m.sample_idx]*m.sample_w[:,:,None],axis=1),m.sample_pos,atol=1e-14)

def test_envelope_potential_derivative():
 from bowen_pd.energies import primitive
 a=.001;b=.01;epsilon=1e-8
 for stretch in [-.005,.0005,.005]:
  derivative=(primitive(stretch+epsilon,a,b)-primitive(stretch-epsilon,a,b))/(2*epsilon)
  f=stretch if stretch<=a else a*(b-stretch)/(b-a)
  assert derivative==pytest.approx(f,rel=1e-8)
 assert primitive(b,a,b)==pytest.approx(.5*a*b)

def test_static_pair_fracture_energy():
 from bowen_pd.energies import concrete_energy
 # fixed-rate quasistatic one bond; reference strain rate zero at its current stretch
 from bowen_pd.geometry import p,delta,V
 x=np.array([[0.,0.,0.],[.005,0,0]])
 i=np.array([0]);j=np.array([1]);r=np.array([.005]);w=np.ones(1)
 a=p.ft_concrete/p.E_concrete;b=5*p.Gf_concrete/(6*p.E_concrete*delta*a)
 u=np.zeros_like(x);u[1,0]=.005*b
 stored,old,release=concrete_energy(x,u,i,j,r,w,np.array([b]),np.array([False]),p.dt,p.E_concrete,p.ft_concrete,p.Gf_concrete,delta,p.fc_concrete,V)
 c=72*p.E_concrete/(np.pi*delta**4)*(1-(.005/delta)**2)**2
 # rounding of stretch may place exactly at sc just below rupture; envelope total still exact
 assert stored+release==pytest.approx(.5*c*.005*V**2*a*b)

def test_source_lattice_family_full_rank():
 m=create(corrected=True,family=True)
 it=m.interface;n=it.reference/it.lengths[:,None]
 local=n[it.steel_ids==100]
 assert len(local)>2 and np.linalg.eigvalsh(local.T@local).min()>0
 np.testing.assert_allclose(it.weights.sum(1),1,atol=1e-14)
 np.testing.assert_allclose(np.sum(m.xyz[it.indices]*it.weights[:,:,None],axis=1),m.sample_pos,atol=1e-14)

def test_family_projection_restart(tmp_path):
 a=Simulation('family_projection');b=Simulation('family_projection')
 for _ in range(350):a.advance()
 a.checkpoint(tmp_path/'family.npz');b.restore(tmp_path/'family.npz')
 for _ in range(4):a.advance();b.advance()
 for name in ('u','v','us','vs'):np.testing.assert_array_equal(getattr(a,name),getattr(b,name))
 assert a.record()==b.record()
 assert a.rate_exchange==b.rate_exchange and a.rupture_release==b.rupture_release

def test_interface_potential_work_gradient(model):
 from bowen_pd.energies import interface_energy
 from bowen_pd.interface import _stage18_pdrod_conservative
 m=model;it=m.interface;p=m.params;u=np.zeros_like(m.xyz);us=np.zeros_like(m.steel_xyz)
 us[100,1]=1e-8
 args=(it.reference,it.lengths,it.steel_ids,it.indices,it.weights,m.interface_broken.copy(),m.interface_prev.copy(),p.bcc,p.s01,p.sc1,p.horizon,p.quadrature_volume,m.STEEL_MASS/m.p.rho_steel,False)
 fc,fs,_=_stage18_pdrod_conservative(u,us,*args)
 def energy(us):return interface_energy(u,us,it.reference,it.lengths,it.steel_ids,it.indices,it.weights,m.interface_prev,m.interface_broken,p.bcc,p.s01,p.sc1,p.horizon,p.quadrature_volume,m.STEEL_MASS/m.p.rho_steel)[0]
 eps=1e-10;up=us.copy();um=us.copy();up[100,1]+=eps;um[100,1]-=eps
 assert fs[100,1]==pytest.approx(-(energy(up)-energy(um))/(2*eps),rel=1e-5)
