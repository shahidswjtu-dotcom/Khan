x = (np.arange(Nx)+0.5)*dx
y = (np.arange(Ny)+0.5)*dx
z = (np.arange(Nz)+0.5)*dx
X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
xyz = np.column_stack((X.ravel(),Y.ravel(),Z.ravel())).astype(np.float64)
N = len(xyz); Vpoint=dx**3
assert N==33600
assert np.unique(xyz,axis=0).shape[0]==N
assert np.all(xyz.min(axis=0)>0) and np.all(xyz.max(axis=0)<[L,W,H])
assert np.isclose(N*Vpoint,L*W*H)
print(f'PASS: {N:,} unique points; volume {N*Vpoint:.7f} m³ (target {L*W*H:.7f} m³)')
print('Coordinate min:',xyz.min(axis=0),'max:',xyz.max(axis=0))
