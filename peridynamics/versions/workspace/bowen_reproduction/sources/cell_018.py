import numpy as np
import matplotlib.pyplot as plt

EPS_QS=3.0e-6  # thesis Eq. (2-13/14): s^-1
F0=10.0e6     # thesis parameter: Pa

def normal_concrete_difs(stretch_rate, fc=p.fc_concrete):
    """Return DIF_Nt and DIF_E from thesis Eqs. (2-13)-(2-14).

    Positive tensile stretch rate is used; compressive/nonpositive rate gives 1.
    Note: the thesis does not completely specify the numerical treatment of
    negative stretch rates; this is a stated implementation convention.
    """
    rate=np.maximum(np.asarray(stretch_rate,dtype=float),0.0)
    exponent=1.0/(10.0+6.0*fc/F0)
    theta=10.0**(6.0*exponent-2.0)
    ratio=np.maximum(rate/EPS_QS,1.0)
    dif_t=np.where(rate<=EPS_QS,1.0,np.where(rate<=1.0,ratio**exponent,theta*ratio**(1.0/3.0)))
    dif_e=np.where(rate<=EPS_QS,1.0,ratio**0.016)
    return dif_t,dif_e

def thresholds_at_rate(stretch_rate, E=p.E_concrete, ft=p.ft_concrete,
                       Gf=p.Gf_concrete, delta=p.horizon_ratio*p.dx):
    """Thesis Eqs. (2-6), (2-8), (2-15), evaluated from *original* constants."""
    s0_ref=ft/E
    sc_ref=5.0*Gf/(6.0*E*delta*s0_ref)
    dif_t,dif_e=normal_concrete_difs(stretch_rate)
    factor=dif_t/dif_e
    return np.asarray(s0_ref*factor),np.asarray(sc_ref*factor)

def micromodulus(distance, delta=p.horizon_ratio*p.dx, E=p.E_concrete, improved=True):
    """Thesis Eq. (2-12) (improved), or Eq. (2-4) (constant)."""
    distance=np.asarray(distance,dtype=float)
    if improved:
        c=72*E/(np.pi*delta**4)*np.maximum(0.0,1.0-(distance/delta)**2)**2
    else:
        c=np.full_like(distance,12*E/(np.pi*delta**4))
    return np.where((distance>0)&(distance<=delta),c,0.0)

print('Quasistatic s0, sc:',thresholds_at_rate(0.0))
for rate in [0.0,3e-6,1e-3,1.0,10.0,100.0]:
    dft,dfe=normal_concrete_difs(rate)
    a,b=thresholds_at_rate(rate)
    print(f'rate={rate:9.3g}/s: DIF_t={float(dft):7.4f}, DIF_E={float(dfe):7.4f}, s0={float(a):.6e}, sc={float(b):.6e}')
