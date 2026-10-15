support_x = np.array([(p.L-p.support_span)/2,(p.L+p.support_span)/2])  # DERIVED
bar_z = 0.010                          # ASSUMED from illustrative 80 mm dimension, verify
bar_y = np.array([0.020,0.050])        # ASSUMED (not dimensioned in source)
bar_x = np.array([0.020,p.L-0.020])     # ASSUMED: rod termination / anchorage length
bar_x_nodes = np.linspace(bar_x[0],bar_x[1],round((bar_x[1]-bar_x[0])/dx)+1)
steel_xyz=np.array([[xx,yy,bar_z] for yy in bar_y for xx in bar_x_nodes])
impact_xy=np.array([p.L/2,p.W/2])
illustrative_gap=0.001                 # ASSUMED only for illustration
impact_center=np.array([impact_xy[0],impact_xy[1],p.H+p.impactor_radius+illustrative_gap])
print('Support x positions:',support_x,'m; spacing:',support_x[1]-support_x[0],'m')
print('ASSUMED transverse bar centers:',bar_y,'m; ASSUMED bar depth z:',bar_z,'m')
print('ASSUMED rod node count:',len(steel_xyz))
print('Illustrative initial sphere center:',impact_center,'m')
assert np.isclose(np.diff(support_x)[0],0.5)
assert len(bar_y)==p.n_bars
assert np.all((bar_y>0)&(bar_y<W))
assert bar_z>p.bar_diameter/2
