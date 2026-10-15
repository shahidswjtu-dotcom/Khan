"""Sufficient linear timestep bound for the initial elastic assembly.

Rayleigh/Jensen inequality: lambda_max <= 2 max_i(sum spring shares / mass_i).
Self-contact and subsequent nonlinear geometry are outside this initial bound.
"""
import json
from pathlib import Path
import numpy as np
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction
p=Benchmark(); m=Assembly(p,Reconstruction.bulk_interface_hypothesis(p))
n=len(m.x); total=n+len(m.sx)
stiffness=np.zeros(total)
for side in (0,1): np.add.at(stiffness,m.b.ends[:,side],m.coefficient/m.b.length)
for side in (0,1): np.add.at(stiffness,n+m.rod_ends[:,side],p.steel_young*m.area/m.rod_length)
k=m.interface_coefficient/m.interface_length
np.add.at(stiffness,n+m.owners,k)
for corner in range(8): np.add.at(stiffness,m.stencils[:,corner],k*m.weights[:,corner])
mass=np.concatenate((m.mass,m.steel_mass))
bound=2*np.max(stiffness/mass)
sufficient=2/np.sqrt(bound)
report={'initial_lambda_upper_bound_per_s2':float(bound),
        'sufficient_linear_dt_s':float(sufficient),'used_dt_s':p.dt,
        'used_dt_over_bound':float(p.dt/sufficient),
        'passed_initial_linear_bound':bool(p.dt<sufficient),
        'limitation':'Initial elastic concrete/interface/rod assembly only; excludes evolving nonlinear contact and force-peak convergence.'}
Path('outputs/stability.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
