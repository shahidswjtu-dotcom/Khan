# Standalone Stage 4 inputs from earlier cells: xyz, dx, p, steel_xyz, bar_x_nodes, bar_y, V.
steel_xyz=np.asarray(steel_xyz,dtype=np.float64)
Nsteel=len(steel_xyz)
rod_pairs=np.concatenate([np.column_stack([np.arange(k*len(bar_x_nodes),(k+1)*len(bar_x_nodes)-1),
                                            np.arange(k*len(bar_x_nodes)+1,(k+1)*len(bar_x_nodes))])
                          for k in range(len(bar_y))]).astype(np.int32)
steel_E=200e9
steel_diameter=0.006  # nominal Chapter 3 reinforcement; confirm against original figure
steel_A=np.pi*steel_diameter**2/4
steel_lengths=np.linalg.norm(steel_xyz[rod_pairs[:,1]]-steel_xyz[rod_pairs[:,0]],axis=1)
assert np.all(steel_lengths>0)
print('Steel nodes:',Nsteel,'rod elements:',len(rod_pairs),'nominal area [mm²]:',steel_A*1e6)
