# Verification 1: source-to-physical steel force equivalence.
fast_interface.reset()
rng=np.random.default_rng(1818)
u_test=rng.normal(0,1.e-9,xyz.shape)
s_test=rng.normal(0,1.e-9,steel_xyz.shape)
fc_phys,fs_phys,nf=conservative_interface(u_test,s_test)
fc_source,fs_source,nf2=fast_interface.evaluate(u_test,s_test,commit=False)
np.testing.assert_allclose(fs_phys,fs_source*STEEL_CONTROL_VOLUME[:,None],rtol=5.e-7,atol=1.e-7)
# Verification 2: physical force and moment balance.
fr=fc_phys.sum(axis=0)+fs_phys.sum(axis=0)
force_scale=max(np.sum(np.linalg.norm(fs_phys,axis=1)),1e-30)
assert np.linalg.norm(fr)/force_scale <1.e-10
# Moment must be computed about current DEFORMED steel & interpolated virtual nodes;
# net torque of a central bond is zero if virtual point weights reproduce position.
cur_c=xyz+u_test;cur_s=steel_xyz+s_test
M=np.cross(cur_c,fc_phys).sum(axis=0)+np.cross(cur_s,fs_phys).sum(axis=0)
mom_scale=max(force_scale*delta,1e-30)
print('Normalized physical force residual:',np.linalg.norm(fr)/force_scale)
print('Normalized physical moment residual:',np.linalg.norm(M)/mom_scale)
assert np.linalg.norm(M)/mom_scale<1.e-8
# Verification 3: undeformed and uniform translation states.
z_c=np.zeros_like(xyz);z_s=np.zeros_like(steel_xyz)
f0,s0,_=conservative_interface(z_c,z_s)
assert np.max(np.abs(f0))<1e-9 and np.max(np.abs(s0))<1e-9
trans=np.array([3.e-5,-1.e-5,2.e-5])
ft,st,_=conservative_interface(z_c+trans,z_s+trans)
assert np.max(np.abs(ft))<1e-4 and np.max(np.abs(st))<1e-4
print('PASS: source scalar-force regression, physical force/moment conservation, rigid translation')