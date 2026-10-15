from dataclasses import dataclass
@dataclass(frozen=True)
class AppendixCInterfaceParameters:
    bcc: float    # reference micromodulus (N/m^6 in a force-density interpretation)
    s01: float    # interface tensile-onset stretch
    sc1: float    # interface critical failure stretch
    horizon: float
    quadrature_volume: float   # Appendix C variable `vol`, NOT yet identified as exact interface quadrature

# Demonstration inputs ONLY, derived from concrete constants; do NOT claim Bowen's interface calibration.
DEMO_INTERFACE=AppendixCInterfaceParameters(
    bcc=float(72*p.E_concrete/(np.pi*(p.horizon_ratio*dx)**4)),
    s01=float(p.ft_concrete/p.E_concrete),
    sc1=float(5*p.Gf_concrete/(6*p.E_concrete*(p.horizon_ratio*dx)*(p.ft_concrete/p.E_concrete))),
    horizon=float(p.horizon_ratio*dx),
    quadrature_volume=float(V))
print('Demo interface parameters (NOT confirmed Bowen):', DEMO_INTERFACE)
