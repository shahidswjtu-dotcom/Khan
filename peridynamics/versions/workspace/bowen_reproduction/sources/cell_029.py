from scipy.spatial import cKDTree
from time import perf_counter
try:
    from numba import njit
    HAVE_NUMBA=True
except ImportError:
    HAVE_NUMBA=False
    def njit(*args,**kwargs):
        if args and callable(args[0]): return args[0]
        return lambda fun: fun

# Full Bowen discretization: 120 x 14 x 20 = 33,600 concrete particles.
# The 1e-10 relative tolerance retains exact-horizon lattice bonds.
BUILD_FULL_NETWORK=True
if BUILD_FULL_NETWORK:
    t0=perf_counter()
    tree=cKDTree(xyz)
    ij=tree.query_pairs(r=delta*(1.0+1e-10),output_type='ndarray').astype(np.int32)
    I,J=ij[:,0].copy(),ij[:,1].copy()
    R0=xyz[J]-xyz[I]
    LEN=np.linalg.norm(R0,axis=1)
    UNIT=R0/LEN[:,None]
    counts=np.bincount(np.concatenate([I,J]),minlength=N)
    print(f'Network: N={N:,} particles, {len(I):,} undirected bonds; neighbors min/mean/max={counts.min()}/{counts.mean():.2f}/{counts.max()}')
    print(f'Build time: {perf_counter()-t0:.2f} s; bond length range {LEN.min()*1e3:.2f}-{LEN.max()*1e3:.2f} mm')
    assert N==33600 and I.size>500000
    assert np.all((LEN>0)&(LEN<=delta*(1.0+1e-9)))
    assert np.all(I<J) and np.all(counts>0)
