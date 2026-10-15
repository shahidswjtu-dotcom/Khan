# Stage 7 independent verification of SOURCE-SHAPED interface law.
import numpy.testing as npt
rng7=np.random.default_rng(77)
uc0=np.zeros_like(xyz);us0=np.zeros_like(steel_xyz)
fc,fs,s,damage=pdrod_c.evaluate(uc0,us0)
assert np.max(np.abs(fc))==0 and np.max(np.abs(fs))==0
npt.assert_allclose(s,0,atol=1e-12)

trans=np.array([1e-5,2e-5,-3e-5])
fc,fs,_,_=pdrod_c.evaluate(uc0+trans,us0+trans)
assert np.linalg.norm(fc)+np.linalg.norm(fs)<1e-3

uc=rng7.normal(0,1e-8,xyz.shape); us=rng7.normal(0,1e-8,steel_xyz.shape)
fc,fs,_,_=pdrod_c.evaluate(uc,us)
net=fc.sum(axis=0)+fs.sum(axis=0)
scale=np.abs(fc).sum()+np.abs(fs).sum()
assert np.linalg.norm(net)<2e-12*max(scale,1.0)

# The forces are central to the interpolated location, so the total moment must also cancel.
mom=np.cross(xyz+uc,fc).sum(axis=0)+np.cross(steel_xyz+us,fs).sum(axis=0)
assert np.linalg.norm(mom)<1e-10*max(scale,1.0)

# Stretch beyond sc1: the bond breaks irreversibly; on unloading it cannot heal.
u_break=np.zeros_like(xyz);u_break[sample_idx[0]]=0.0
# Isolated one-sample reproduction removes coupling overlap while retaining the 8-point map.
unit=AppendixCPDROD(steel_xyz,sample_pos[:1],sample_idx[:1],sample_w[:1],sample_steel[:1],DEMO_INTERFACE)
u_steel=np.zeros_like(steel_xyz)
u_steel[unit.steel_ids[0]]=-2*DEMO_INTERFACE.sc1*unit.reference[0]
_,_,s_test,fail=unit.evaluate(uc0,u_steel,commit=True)
assert fail[0] and unit.broken[0]
_,fs_unload,_,fail_again=unit.evaluate(uc0,us0,commit=False)
assert fail_again[0] and np.linalg.norm(fs_unload)==0
print('PASS: zero state; translation; equal/opposite force; zero total moment; irreversible rupture')
print('PASS: Appendix C source-shaped nonlinear interface unit tests')
