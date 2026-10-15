import numpy as np
from printed_pd.concrete import _fast_concrete_kernel
from printed_pd.interface import _stage18_pdrod_conservative

def test_printed_concrete_prior_flag_does_not_gate_force():
 x=np.array([[0.,0,0],[.005,0,0]]);u=np.zeros_like(x)
 broken=np.array([True]);prev=np.array([0.])
 u[1,0]=1e-9
 F,_=_fast_concrete_kernel(x,u,np.array([0]),np.array([1]),np.array([.005]),np.ones(1),prev,broken,3e-7,25e9,4.6e6,105.,.015,60e6,(.005**3)**2,True)
 assert np.max(abs(F))>0 and broken[0]
 np.testing.assert_allclose(F.sum(0),0)

def test_printed_interface_prior_flag_persists_but_force_reloads():
 uc=np.zeros((8,3));us=np.zeros((1,3));us[0,1]=1e-8
 ref=np.array([[0.,.005,0.]])
 broken=np.array([True]);previous=np.zeros(1)
 fc,fs,_=_stage18_pdrod_conservative(uc,us,ref,np.array([.005]),np.array([0]),np.arange(8)[None,:],np.full((1,8),.125),broken,previous,1e19,.000184,.001268,.015,.005**3,np.array([1.4e-7]),True)
 assert np.max(abs(fs))>0 and broken[0]
 np.testing.assert_allclose(fc.sum(0)+fs.sum(0),0,atol=1e-12)
