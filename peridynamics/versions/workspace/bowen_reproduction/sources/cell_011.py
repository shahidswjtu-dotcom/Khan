delta=p.horizon_ratio*dx
s0=p.ft_concrete/p.E_concrete
sc=5*p.Gf_concrete/(6*p.E_concrete*delta*s0)
c_constant=12*p.E_concrete/(np.pi*delta**4)
def c_weighted(r):
    r=np.asarray(r)
    return np.where(r<=delta,72*p.E_concrete/(np.pi*delta**4)*(1-(r/delta)**2)**2,0.0)
print(f'Initial stretch s0 = {s0:.8e}')
print(f'Critical stretch sc = {sc:.8e}')
print(f'Constant micromodulus = {c_constant:.8e} [model units]')
print(f'Weighted c(1 dx) = {float(c_weighted(dx)):.8e} [model units]')
assert 0<s0<sc
