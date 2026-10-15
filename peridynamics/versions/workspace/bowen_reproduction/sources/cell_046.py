rng=np.random.default_rng(2026)
uc=np.zeros_like(xyz);us=np.zeros_like(steel_xyz)
fc,fs=coupling_forces(uc,us)
assert np.max(np.abs(fc))==0 and np.max(np.abs(fs))==0
print('PASS 1: undeformed interface has zero force')

# Common translation of steel and concrete: no artificial coupling traction.
translation=np.array([1.1e-5,-1.5e-5,2.0e-5])
fc,fs=coupling_forces(np.tile(translation,(N,1)),np.tile(translation,(Nsteel,1)))
assert np.linalg.norm(fc)+np.linalg.norm(fs)<1e-6
print('PASS 2: rigid translation gives zero coupling force')

# Affine displacement reproduced exactly by interpolation.
A=np.array([[1e-5,2e-6,0],[3e-6,-1e-5,1e-6],[0,0,2e-5]])
ucon=xyz@A.T; usteel=steel_xyz@A.T
interpolated=np.sum(ucon[sample_idx]*sample_w[:,:,None],axis=1)
np.testing.assert_allclose(interpolated,sample_pos@A.T,rtol=0,atol=5e-17)
print('PASS 3: trilinear interpolation exactly reproduces affine displacement')

# Net force, moment and discrete virtual work for arbitrary (small) displacements.
ucon=rng.normal(0,1e-7,xyz.shape); usteel=rng.normal(0,1e-7,steel_xyz.shape)
fc,fs=coupling_forces(ucon,usteel)
net=fc.sum(axis=0)+fs.sum(axis=0)
scale=np.sum(np.abs(fc))+np.sum(np.abs(fs))
assert np.linalg.norm(net)<1e-9*max(scale,1)
# Noncoincident coupling bonds do exert a physical couple if the chosen force is not central.
# The identity checked is the work conjugacy of interpolation/force transfer, not zero net torque.
ducon=rng.normal(0,1e-8,xyz.shape); dusteel=rng.normal(0,1e-8,steel_xyz.shape)
lhs=np.sum(fc*ducon)+np.sum(fs*dusteel)
rel=np.sum(ucon[sample_idx]*sample_w[:,:,None],axis=1)-usteel[sample_steel]
drel=np.sum(ducon[sample_idx]*sample_w[:,:,None],axis=1)-dusteel[sample_steel]
rhs=np.sum(COUPLING_STIFFNESS*rel*drel)
np.testing.assert_allclose(lhs,-rhs,rtol=1e-10,atol=1e-12)
print('PASS 4: coupling action/reaction and virtual work:',np.linalg.norm(net),'N')

fr=rod_forces(usteel)
assert np.linalg.norm(fr.sum(axis=0))<1e-8
print('PASS 5: rod axial forces are equal/opposite; net force:',np.linalg.norm(fr.sum(axis=0)),'N')
print('STAGE 4 UNIT TESTS PASSED (not impact validation).')
