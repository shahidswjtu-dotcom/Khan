# Direct comparison: source-literal interface vs earlier irreversible interface
# in elastic range, both should agree numerically (same demo inputs and mesh).
uf=np.zeros_like(xyz); sf=np.zeros_like(steel_xyz)
r0=appendix_c_beam.evaluate(uf,sf,mode='literal')
assert np.max(np.abs(r0['concrete_force']))==0.
assert np.max(np.abs(r0['steel_force']))==0.

rng14=np.random.default_rng(147)
uf=rng14.normal(0,1.e-8,xyz.shape);sf=rng14.normal(0,1.e-8,steel_xyz.shape)
r1=appendix_c_beam.evaluate(uf,sf,mode='literal')
fast_interface.reset()
fc_old,fs_old,_=fast_interface.evaluate(uf,sf,False)
np.testing.assert_allclose(r1['concrete_force'],fc_old,rtol=1e-8,atol=2e-4)
np.testing.assert_allclose(r1['steel_force'],fs_old,rtol=1e-8,atol=2e-4)
assert np.linalg.norm(r1['concrete_force'].sum(axis=0)+r1['steel_force'].sum(axis=0))<1e-4
print('PASS zero state, source-vs-previous elastic force agreement, and total force equilibrium')
