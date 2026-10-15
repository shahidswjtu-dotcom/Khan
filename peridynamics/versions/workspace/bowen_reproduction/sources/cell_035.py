import numpy.testing as npt
t0=perf_counter()
zero=np.zeros_like(xyz)
Fzero,s,broken=model.evaluate(zero,p.dt)
assert np.max(np.abs(Fzero))<1e-7
print('PASS 1: undeformed state has zero internal force')

# Nonuniform, smooth displacement: small relative to stretching threshold.
ux=1.e-6*np.sin(np.pi*xyz[:,0]/p.L)*np.sin(np.pi*xyz[:,2]/p.H)
u=np.column_stack([ux,0.4*ux,-0.2*ux])
F,s,broken=model.evaluate(u,p.dt,override_rate=0.0)
net=F.sum(axis=0)
scale=np.sum(np.linalg.norm(F,axis=1))
assert np.linalg.norm(net) <= max(1.e-8, 1.e-10*scale)
# Central forces imply vanishing internal torque about arbitrary fixed origin.
torque=np.sum(np.cross(xyz+u,F),axis=0)
assert np.linalg.norm(torque)<max(1.e-8, 1.e-9*scale*p.L)
print(f'PASS 2: net internal force {np.linalg.norm(net):.3e} N, torque {np.linalg.norm(torque):.3e} N m')

# Uniform infinitesimal isotropic dilation, below first damage threshold.
eps=1.e-5
u_aff=eps*(xyz-np.array([p.L/2,p.W/2,p.H/2]))
Faff,saff,baff=model.evaluate(u_aff,p.dt,override_rate=0.0)
force_mag=np.linalg.norm(Faff,axis=1)
print('Interior affine-force max [N]:',force_mag[interior].max(),'boundary max [N]:',force_mag[~interior].max())
assert force_mag[interior].max()<1.e-7
assert not np.any(baff)
print('PASS 3: affine-strain interior force cancellation')

# Fracture test: force remains zero even after unloading to intact stretch.
# Shift one corner particle far away, mark all its connected bonds failed.
corner=0
u_break=np.zeros_like(xyz); u_break[corner]=[0.08,0.0,0.0]
Fb,sb,bb=model.evaluate(u_break,p.dt,commit=True,override_rate=0.0)
assert np.any(model.broken)
broken_before=model.broken.copy()
Fu,su,bu=model.evaluate(zero,p.dt,commit=True,override_rate=0.0)
assert np.all(model.broken[broken_before])
assert np.max(np.abs(Fu))<1.e-7
print('PASS 4: damage remains after unloading;',broken_before.sum(),'bonds broken')
# Reset history before any physical simulation.
model.broken[:]=False;model.prev_stretch[:]=0.0
print(f'All Stage 3 assembly tests passed. Elapsed {perf_counter()-t0:.2f} s')
