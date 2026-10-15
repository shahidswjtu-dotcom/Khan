def concrete_grid_eight(point):
    point=np.asarray(point,dtype=float)
    nx,ny,nz=Nx,Ny,Nz
    h=dx
    a=point/h-0.5
    base=np.floor(a).astype(int)
    if np.any(base<0) or np.any(base>=np.array([nx-1,ny-1,nz-1])):
        raise ValueError(f'Coupling point {point} outside center-grid interpolation domain')
    frac=a-base
    inds=[]; ws=[]
    for qx in (0,1):
        for qy in (0,1):
            for qz in (0,1):
                inds.append(np.ravel_multi_index((base[0]+qx,base[1]+qy,base[2]+qz),(nx,ny,nz)))
                ws.append((frac[0] if qx else 1-frac[0])*(frac[1] if qy else 1-frac[1])*(frac[2] if qz else 1-frac[2]))
    return np.array(inds,dtype=np.int32),np.array(ws,dtype=float)

# Couple steel nodes at their current locations to interpolated concrete kinematic points,
# with an optionally offset sampling location to avoid a zero-length bond.
# A transverse offset on each rod side allows a finite interface bond direction.
COUPLING_OFFSET=dx      # a provisional geometric choice, not established from Bowen's source
COUPLING_STIFFNESS=2.0e8 # N/m per coupling location: provisional, not calibrated
sample_steel=[];sample_pos=[];sample_idx=[];sample_w=[]
for si, pos in enumerate(steel_xyz):
    for sign in (-1.,1.):
        pt=pos+np.array([0.,sign*COUPLING_OFFSET,0.])
        idx,w=concrete_grid_eight(pt)
        sample_steel.append(si);sample_pos.append(pt);sample_idx.append(idx);sample_w.append(w)
sample_steel=np.asarray(sample_steel,dtype=np.int32)
sample_pos=np.asarray(sample_pos,dtype=float)
sample_idx=np.asarray(sample_idx,dtype=np.int32)
sample_w=np.asarray(sample_w,dtype=float)
assert np.all(sample_w>=-1e-14)
np.testing.assert_allclose(sample_w.sum(axis=1),1.,atol=1e-13)
np.testing.assert_allclose((xyz[sample_idx]*sample_w[:,:,None]).sum(axis=1),sample_pos,atol=1e-13)
print('Coupling samples:',len(sample_steel),'eight-node maps:',sample_idx.shape)
