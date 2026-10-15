import numpy.testing as npt
reset_stage8()
zuc=np.zeros_like(xyz); zus=np.zeros_like(steel_xyz)
f0,s0,h0,E0,n0,pen0,nb0,ni0=source_shaped_forces(zuc,zus,hammer_z_initial,p.dt,commit=False)
assert np.max(np.abs(f0))==0 and np.max(np.abs(s0))==0
assert n0==0 and nb0==0 and ni0==0
# Rigid translation together of concrete/steel, without contact: no internal forces.
rig=np.array([2e-7,-1e-7,3e-7])
fr,sr,*_=source_shaped_forces(zuc+rig,zus+rig,hammer_z_initial+rig[2],p.dt,commit=False)
assert np.max(np.abs(fr))<1e-4 and np.max(np.abs(sr))<1e-4
# Apply a tiny relative steel movement; interface reactions must balance.
pert=zus.copy();pert[:,0]=1e-8
fc,fs,h,E,ns,mp,_,_=source_shaped_forces(zuc,pert,hammer_z_initial,p.dt,commit=False)
residual=np.linalg.norm(fc.sum(axis=0)+fs.sum(axis=0)+np.array([0.,0.,h]))
scale=max(np.linalg.norm(fc.sum(axis=0)),np.linalg.norm(fs.sum(axis=0)),1.)
assert residual/scale<1e-9, (residual,scale)
assert np.linalg.norm(fs)>0
assert np.count_nonzero(interface8.broken)==0  # read-only calls did not commit
print('PASS: zero-state and common rigid translation')
print('PASS: Stage 7 PDROD participates in integrated forces')
print(f'PASS: integrated pairwise action–reaction, relative residual={residual/scale:.3e}')
reset_stage8()