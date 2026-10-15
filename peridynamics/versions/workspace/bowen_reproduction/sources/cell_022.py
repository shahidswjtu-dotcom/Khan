delta=p.horizon_ratio*p.dx
r=delta*0.5
s0,sc=map(float,thresholds_at_rate(0.0))
c=float(micromodulus(r))

f,fail=bond_scalar_force(np.array([s0/2,s0,s0+(sc-s0)/2,sc+1.e-7]),0.0,r)
assert np.isclose(f[0],c*s0/2,rtol=1e-12)
assert np.isclose(f[1],c*s0,rtol=1e-12)
assert np.isclose(f[2],c*s0/2,rtol=1e-12)
assert f[3]==0 and fail[-1]
assert not np.any(fail[:3])

h=BondHistory(1)
assert h.update([sc*1.1],[0.0],[r])[0]==0
assert h.update([s0*0.1],[0.0],[r])[0]==0  # no healing
assert h.broken[0]

dft,dfe=normal_concrete_difs(1e-2)
a,b=thresholds_at_rate(1e-2)
assert np.isclose(a/s0,dft/dfe) and np.isclose(b/sc,dft/dfe)
assert np.all(np.asarray(normal_concrete_difs(np.array([0.0,1e-3,10.0]))[0])>=1.0)
assert float(micromodulus(delta*1.1))==0.0
print('PASS: bond elastic slope, softening, irreversible failure, rate factors, vectorized shape.')
