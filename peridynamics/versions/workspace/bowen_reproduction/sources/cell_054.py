def sphere_contact(xyz_def, hammer_center_z, v_concrete=None, vz_hammer=0.0, k=CONTACT_STIFFNESS):
    center=np.array([hammer_xy[0],hammer_xy[1],hammer_center_z],dtype=float)
    d=xyz_def-center[None,:]
    distances=np.linalg.norm(d,axis=1)
    penetration=np.maximum(HAMMER_RADIUS-distances,0.)
    active=penetration>0
    forces=np.zeros_like(xyz_def)
    if active.any():
        unit=d[active]/distances[active,None]
        magnitudes=k*penetration[active]
        forces[active]=magnitudes[:,None]*unit
    f_hammer = -np.sum(forces,axis=0)
    potential=0.5*k*np.sum(penetration[active]**2)
    return forces, f_hammer, float(potential),int(active.sum()),float(penetration.max(initial=0.0))

# Standalone test: one point 0.25 mm inside the sphere surface, and out-of-contact configuration.
point=np.array([[p.L/2,p.W/2,0.0]])
f0,h0,E0,n0,pen0=sphere_contact(point,100.0)
assert n0==0 and np.all(f0==0) and np.all(h0==0)
cz=HAMMER_RADIUS - 0.00025
f1,h1,E1,n1,pen1=sphere_contact(point,cz)
assert n1==1 and f1[0,2]<0 and h1[2]>0
np.testing.assert_allclose(f1.sum(axis=0)+h1,0.0,atol=1e-10)
np.testing.assert_allclose(E1,0.5*CONTACT_STIFFNESS*pen1**2,rtol=1e-9)
print('PASS: hammer/contact forces are equal and opposite; positive repulsion and potential verified')
print('Single-point penetration (mm):',pen1*1e3,'; force (kN):',abs(h1[2])/1e3)
