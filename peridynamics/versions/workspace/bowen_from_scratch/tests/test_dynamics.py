import numpy as np
import pytest
from pdbeam.operators import rods,dynamic_thresholds,law_energy,repulsion
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction
from pdbeam.solver import Simulation

def test_rod_objectivity_and_plastic_return():
    x=np.array([[0.,0.,0.],[1.,0.,0.]])
    ends=np.array([[0,1]]); lengths=np.array([1.]); plastic=np.zeros(1)
    rotated=np.array([[0.,0.,0.],[0.,1.,0.]])
    force,energy,diss=rods(x,rotated-x,ends,lengths,2.,10.,1.,plastic)
    np.testing.assert_allclose(force,0,atol=1e-12)
    assert energy==0 and diss==0
    u=np.array([[0.,0.,0.],[.2,0.,0.]])
    force,energy,diss=rods(x,u,ends,lengths,2.,10.,1.,plastic)
    assert force[0,0]==pytest.approx(2)
    assert plastic[0]==pytest.approx(.1)
    assert diss==pytest.approx(.2)
    force,_,diss=rods(x,u,ends,lengths,2.,10.,1.,plastic)
    assert diss==pytest.approx(0,abs=1e-14)


def test_bilinear_potential_gradient():
    for s in (-.001,.00005,.0003):
        force,energy=law_energy(s,.0001,.001)
        h=1e-9
        gradient=(law_energy(s+h,.0001,.001)[1]-law_energy(s-h,.0001,.001)[1])/(2*h)
        assert gradient==pytest.approx(force,rel=1e-7)


def test_signed_rate_quasistatic_limit():
    assert dynamic_thresholds(-100.,.1,.3)==(.1,.3)
    assert dynamic_thresholds(3e-6,.1,.3)==(.1,.3)
    a,b=dynamic_thresholds(10.,.1,.3)
    assert b/a==pytest.approx(3)


def test_self_contact_potential_and_action_reaction():
    ref=np.array([[0.,0.,0.],[.005,0.,0.]])
    pos=np.array([[0.,0.,0.],[.003,0.,0.]])
    ends=np.array([[0,1]])
    force,energy=repulsion(ref,pos,ends,.005,.015,25e9)
    assert force[0,0]<0 and energy>0
    np.testing.assert_allclose(force.sum(axis=0),0)
    h=1e-9; plus=pos.copy(); minus=pos.copy()
    plus[1,0]+=h; minus[1,0]-=h
    gradient=(repulsion(ref,plus,ends,.005,.015,25e9)[1]-repulsion(ref,minus,ends,.005,.015,25e9)[1])/(2*h)
    assert -gradient==pytest.approx(force[1,0],rel=1e-7)


def test_unresolved_interface_cannot_run():
    with pytest.raises(ValueError,match='interface inputs'): Assembly()


def test_assembled_contact_restart_and_momentum(tmp_path):
    p=Benchmark(); model=Assembly(p,Reconstruction.bulk_interface_hypothesis(p))
    direct=Simulation(model)
    direct.hz-=.00265 # synthetic initial penetration to exercise contact and constraints immediately
    row=direct.advance()
    assert row['force_N']>0
    assert row['momentum_residual_kg_m_s']<1e-9
    direct.save(tmp_path/'checkpoint.npz')
    resumed=Simulation(model); resumed.restore(tmp_path/'checkpoint.npz')
    for _ in range(2):
        a=direct.advance(); b=resumed.advance()
        assert a.keys()==b.keys()
        np.testing.assert_equal(list(a.values()),list(b.values()))
    for name in ('u','v','su','sv','alive','interface_alive','plastic'):
        np.testing.assert_array_equal(getattr(direct,name),getattr(resumed,name))
