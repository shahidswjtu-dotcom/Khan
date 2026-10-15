import numpy as np
import pytest
from pdbeam.core import *

def test_chapter3_geometry_and_mass():
    p=Benchmark(); x=particles(p)
    assert x.shape==(33600,3)
    assert len(x)*p.spacing**3*p.density == pytest.approx(10.08)
    assert p.onset == pytest.approx(.000184)
    assert p.terminal > p.onset


def test_affine_force_units_and_momentum():
    p=Benchmark(); x=np.array([[0.,0.,0.],[p.spacing,0.,0.]])
    b=family(x,p); u=.00001*x
    f,alive=concrete_force(x,u,b,p,np.ones(1,dtype=bool))
    expected=b.stiffness[0]*.00001*p.spacing**6*b.volume_fraction[0]
    assert f[0,0]==pytest.approx(expected)
    assert np.linalg.norm(f.sum(axis=0)) < 1e-12
    assert f[0,0]/(p.density*p.spacing**3)==pytest.approx(expected/(p.density*p.spacing**3))


def test_rigid_motion_has_no_internal_force():
    p=Benchmark(); x=np.array([[0.,0.,0.],[p.spacing,0.,0.]])
    b=family(x,p); angle=.73
    rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    u=x@rotation.T+np.array([1.,2.,3.])-x
    f,_=concrete_force(x,u,b,p,np.ones(1,dtype=bool))
    assert np.max(abs(f)) < 1e-7


def test_failure_history_is_explicit():
    p=Benchmark(); x=np.array([[0.,0.,0.],[p.spacing,0.,0.]])
    b=family(x,p)
    _,flags=concrete_force(x,2*p.terminal*x,b,p,np.ones(1,dtype=bool))
    f,_=concrete_force(x,.5*p.onset*x,b,p,flags)
    assert not flags[0] and np.all(f==0)
    f,_=concrete_force(x,.5*p.onset*x,b,p,flags,history='instantaneous')
    assert f[0,0]>0


def test_interface_reproduces_affine_fields_and_balances_torque():
    vertices=np.stack(np.meshgrid([0.,1.],[0.,1.],[0.,1.],indexing='ij'),axis=-1).reshape(-1,3)
    point=np.array([.2,.3,.7]); w=interpolate_weights(point,vertices)
    assert w.sum()==pytest.approx(1)
    np.testing.assert_allclose(w@vertices,point)
    A=np.array([[1.1,.2,0],[0,.9,.1],[.1,0,1.2]]); shift=np.array([.3,.4,.5])
    np.testing.assert_allclose(w@(vertices@A.T+shift), point@A.T+shift)
    steel=np.array([-.2,.3,.7])
    fs,fc=interface_pair(steel,vertices,w,.3,2.,.4,.5)
    np.testing.assert_allclose(fs+fc.sum(axis=0),0,atol=1e-15)
    np.testing.assert_allclose(np.cross(steel,fs)+np.cross(vertices,fc).sum(axis=0),0,atol=1e-15)


def test_contact_impulse_and_free_flight():
    old=np.array([[0.,0.,-1.1],[2.,0.,0.]])
    pred=np.array([[0.,0.,-.9],[2.,0.,0.]])
    mass=np.array([2.,3.]); dt=.1
    corrected,v,reaction=project_contact(old,pred,np.zeros(3),1.,mass,dt)
    np.testing.assert_allclose(corrected[0],[0,0,-1])
    np.testing.assert_allclose(corrected[1],pred[1])
    np.testing.assert_allclose((mass[:,None]*(v-(pred-old)/dt)).sum(axis=0)+reaction,0,atol=1e-14)


def test_restart_roundtrip_and_config_guard(tmp_path):
    p=Benchmark(); path=tmp_path/'state.npz'; u=np.arange(12).reshape(4,3)
    save_checkpoint(path,p,7,u=u,alive=np.array([True,False]))
    step,state=load_checkpoint(path,p)
    assert step==7
    np.testing.assert_array_equal(state['u'],u)
    with pytest.raises(ValueError): load_checkpoint(path,Benchmark(dt=1e-7))
