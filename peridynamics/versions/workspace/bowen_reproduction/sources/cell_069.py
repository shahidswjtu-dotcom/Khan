rng=np.random.default_rng(641)
uc_a=rng.normal(0,1e-7,size=xyz.shape)
us_a=rng.normal(0,1e-7,size=steel_xyz.shape)

def coupling_potential(uc,us):
    interp=np.sum(uc[sample_idx]*sample_w[:,:,None],axis=1)
    rel=interp-us[sample_steel]
    return .5*COUPLING_STIFFNESS*np.sum(rel*rel)

fc_a,fs_a=coupling_forces(uc_a,us_a)
np.testing.assert_allclose(fc_a.sum(axis=0)+fs_a.sum(axis=0),np.zeros(3),atol=5e-10)
# select an actually mapped concrete material point
ii=int(sample_idx[len(sample_idx)//2,3]);jj=int(sample_steel[len(sample_steel)//2]);eps=1e-8
up=uc_a.copy();um=uc_a.copy();up[ii,2]+=eps;um[ii,2]-=eps
fdc=-(coupling_potential(up,us_a)-coupling_potential(um,us_a))/(2*eps)
np.testing.assert_allclose(fdc,fc_a[ii,2],atol=1e-3,rtol=1e-5)
up=us_a.copy();um=us_a.copy();up[jj,1]+=eps;um[jj,1]-=eps
fds=-(coupling_potential(uc_a,up)-coupling_potential(uc_a,um))/(2*eps)
np.testing.assert_allclose(fds,fs_a[jj,1],atol=1e-3,rtol=1e-5)
print('PASS: interface momentum balance and penalty-potential gradients')
