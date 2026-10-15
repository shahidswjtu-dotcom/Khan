BUILD_NEIGHBORS = False
if BUILD_NEIGHBORS:
    from scipy.spatial import cKDTree
    tree=cKDTree(xyz)
    pairs=tree.query_pairs(r=delta*(1+1e-12), output_type='ndarray')
    distances=np.linalg.norm(xyz[pairs[:,0]]-xyz[pairs[:,1]],axis=1)
    counts=np.bincount(pairs.reshape(-1),minlength=N)
    print('Undirected pairs:',len(pairs),'max / min / mean neighbors:',counts.max(),counts.min(),counts.mean())
    assert np.max(distances) <= delta*(1+1e-10)
else:
    print('Topology construction skipped. Set BUILD_NEIGHBORS=True after verifying geometry.')
