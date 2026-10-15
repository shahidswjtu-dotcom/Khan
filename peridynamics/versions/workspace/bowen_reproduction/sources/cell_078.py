def bowen_concrete_force_with_modulus_DIF(stretch,stretch_rate,reference_distance,broken=False):
    """Thesis-informed correction candidate: apply DIF_E to elastic micromodulus AND DIF_t/DIF_E to stretches."""
    s,r,l,b=np.broadcast_arrays(np.asarray(stretch,float),np.asarray(stretch_rate,float),
                                 np.asarray(reference_distance,float),np.asarray(broken,bool))
    dif_t,dif_e=normal_concrete_difs(r)
    start=(p.ft_concrete/p.E_concrete)*dif_t/dif_e
    end=(5*p.Gf_concrete/(6*p.E_concrete*(p.horizon_ratio*p.dx)*(p.ft_concrete/p.E_concrete)))*dif_t/dif_e
    coeff=micromodulus(l)*dif_e
    strength=np.where(s<=start,coeff*s,coeff*start*np.clip((end-s)/(end-start),0,1))
    failed=b | (s>=end)
    return np.where(failed,0.,strength),failed

# Source-independent numeric consistency: elastic slope increases with DIF_E.
probe_s=1e-8
probe_l=0.5*p.horizon_ratio*p.dx
for rate in (0.,1e-3,1.):
    strength,_=bowen_concrete_force_with_modulus_DIF(probe_s,rate,probe_l)
    expected=micromodulus(probe_l)*normal_concrete_difs(rate)[1]*probe_s
    npt.assert_allclose(strength,expected,rtol=1e-12)
print('PASS: rate-dependent elastic-slope test')
