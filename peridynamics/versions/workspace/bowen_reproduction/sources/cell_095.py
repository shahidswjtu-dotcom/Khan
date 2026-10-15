# Compare the SAME states in NumPy and compiled implementations, with no history commits.
from time import perf_counter
rng10=np.random.default_rng(110)
for label, amplitude in [('zero',0.0),('tiny',1.e-8),('elastic',1.e-7)]:
    U=rng10.normal(size=xyz.shape)*amplitude
    concrete8.broken.fill(False); concrete8.prev_stretch.fill(0.0)
    fast_concrete.reset()
    baseline=concrete8.evaluate(U,p.dt,commit=False)[0]
    optimized,nf=fast_concrete.evaluate(U,p.dt,commit=False)
    difference=np.max(np.abs(baseline-optimized))
    scale=max(np.max(np.abs(baseline)),1.)
    np.testing.assert_allclose(optimized,baseline,rtol=2.e-10,atol=1.e-7)
    print(f'{label}: max abs error={difference:.4g} N, normalized error={difference/scale:.3g}')

# Separate failure/history consistency test at high stretch, then unload.
U=np.zeros_like(xyz)
U[:,0]=5.e-3*(xyz[:,0]/p.L)**2
concrete8.broken.fill(False); concrete8.prev_stretch.fill(0.0)
fast_concrete.reset()
f_np,_,fail_np=concrete8.evaluate(U,p.dt,commit=True)
f_opt,nf=fast_concrete.evaluate(U,p.dt,commit=True)
assert np.array_equal(concrete8.broken,fast_concrete.broken)
np.testing.assert_allclose(f_np,f_opt,rtol=2.e-10,atol=1.e-7)
f_np2=concrete8.evaluate(np.zeros_like(U),p.dt,commit=True)[0]
f_opt2,nf2=fast_concrete.evaluate(np.zeros_like(U),p.dt,commit=True)
assert np.array_equal(concrete8.broken,fast_concrete.broken)
np.testing.assert_allclose(f_np2,f_opt2,rtol=2.e-10,atol=1.e-7)
print('PASS: irreversible history and unloading equivalence; broken bonds:',nf)
# Return original states to zeros for repeatable studies.
reset_stage8();fast_concrete.reset()