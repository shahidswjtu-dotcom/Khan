import numpy as np
import pytest
from pdbeam.core import Benchmark,family,concrete_force
from pdbeam.operators import concrete
from pdbeam.energy import concrete_potential,threshold_ledger

def pair():
    p=Benchmark(); x=np.array([[0.,0.,0.],[p.spacing,0.,0.]])
    b=family(x,p); coefficient=b.stiffness*p.spacing**6*b.volume_fraction
    return p,x,b,coefficient


def test_compiled_static_force_matches_independent_operator():
    p,x,b,c=pair(); u=.5*p.onset*x
    f,_=concrete_force(x,u,b,p,np.ones(1,dtype=bool))
    compiled,_=concrete(x,u,u,b.ends,b.length,c,np.ones(1,dtype=bool),p.onset,p.terminal,p.dt,True,False)
    np.testing.assert_allclose(f,compiled,rtol=1e-10)


def test_concrete_force_is_negative_energy_gradient():
    p,x,b,c=pair(); u=.5*p.onset*x
    alive=np.ones(1,dtype=bool); a=np.array([p.onset]); terminal=np.array([p.terminal])
    force,_=concrete(x,u,u,b.ends,b.length,c,alive.copy(),p.onset,p.terminal,p.dt,True,False)
    h=1e-10; plus=u.copy(); minus=u.copy(); plus[1,0]+=h; minus[1,0]-=h
    derivative=(concrete_potential(x,plus,b.ends,b.length,c,alive,a,terminal,True)-concrete_potential(x,minus,b.ends,b.length,c,alive,a,terminal,True))/(2*h)
    assert -derivative==pytest.approx(force[1,0],rel=1e-7)


def test_rate_exchange_ledger_at_fixed_geometry():
    p,x,b,c=pair(); u=.5*p.onset*x; previous=np.zeros_like(x)
    alive=np.ones(1,dtype=bool); a=np.array([p.onset]); terminal=np.array([p.terminal])
    before=concrete_potential(x,u,b.ends,b.length,c,alive,a,terminal,True)
    exchange,loss=threshold_ledger(x,u,previous,b.ends,b.length,c,alive,a,terminal,p.onset,p.terminal,p.dt,True,True)
    after=concrete_potential(x,u,b.ends,b.length,c,alive,a,terminal,True)
    assert after-before==pytest.approx(exchange,abs=1e-15)
    assert loss==0


def test_deleted_rupture_energy_ledger():
    p,x,b,c=pair(); u=2*p.terminal*x
    alive=np.ones(1,dtype=bool); a=np.array([p.onset]); terminal=np.array([p.terminal])
    before=concrete_potential(x,u,b.ends,b.length,c,alive,a,terminal,True)
    exchange,loss=threshold_ledger(x,u,u,b.ends,b.length,c,alive,a,terminal,p.onset,p.terminal,p.dt,False,True)
    concrete(x,u,u,b.ends,b.length,c,alive,p.onset,p.terminal,p.dt,True,False)
    after=concrete_potential(x,u,b.ends,b.length,c,alive,a,terminal,True)
    assert after+loss-before-exchange==pytest.approx(0,abs=1e-15)


def test_compiled_literal_history_transmits_force_after_prior_failure():
    p,x,b,c=pair(); u=.5*p.onset*x
    flags=np.zeros(1,dtype=bool)
    f,energy=concrete(x,u,u,b.ends,b.length,c,flags.copy(),p.onset,p.terminal,p.dt,False,False)
    assert f[0,0]>0 and energy>0
    f,energy=concrete(x,u,u,b.ends,b.length,c,flags.copy(),p.onset,p.terminal,p.dt,True,False)
    assert np.all(f==0) and energy==0
