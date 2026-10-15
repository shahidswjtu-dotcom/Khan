V=dx**3
# Independent toggles for sensitivity; defaults are conservative documented choices.
USE_PARTIAL_VOLUME=True
USE_SURFACE_CORRECTION=False
volweight=np.clip((delta+dx/2-LEN)/dx,0.0,1.0) if USE_PARTIAL_VOLUME else np.ones_like(LEN)

# Relative neighborhood occupancy, based on weighted discrete interaction counts.
occupancy=np.bincount(I,weights=volweight,minlength=N)+np.bincount(J,weights=volweight,minlength=N)
interior=(xyz[:,0]>delta+dx)&(xyz[:,0]<p.L-delta-dx)&(xyz[:,1]>delta+dx)&(xyz[:,1]<p.W-delta-dx)&(xyz[:,2]>delta+dx)&(xyz[:,2]<p.H-delta-dx)
reference_neighbors=float(np.median(occupancy[interior]))
node_correction=np.ones(N)
if USE_SURFACE_CORRECTION:
    node_correction=np.clip(reference_neighbors/np.maximum(occupancy,1),1.0,3.0)
# Geometric symmetric mean retains equal/opposite bond forces.
PAIR_WEIGHT=volweight*np.sqrt(node_correction[I]*node_correction[J])
print('Interior weighted neighbor reference:',reference_neighbors)
print('Surface correction enabled:',USE_SURFACE_CORRECTION,'; pair quadrature weights:',PAIR_WEIGHT.min(),PAIR_WEIGHT.max())
assert interior.any()
