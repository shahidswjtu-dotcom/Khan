import numpy as np
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction
from pdbeam.solver import Simulation
from audit_projection_energy import AuditedSimulation,potential_change

def test_projection_audit_does_not_modify_mechanics():
    p=Benchmark(); m=Assembly(p,Reconstruction.bulk_interface_hypothesis(p))
    reference=Simulation(m); audit=AuditedSimulation(m)
    reference.hz-=.00265; audit.hz-=.00265
    a=reference.advance(); b=audit.advance()
    for k in a:
        np.testing.assert_equal(a[k],b[k])
    for name in ('u','v','su','sv','alive','interface_alive','plastic'):
        np.testing.assert_array_equal(getattr(reference,name),getattr(audit,name))
    assert potential_change(m,audit,m.x+audit.u,m.x+audit.u)==0
