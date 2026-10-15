"""Assess the assembled correction, rather than just its calibration fields."""
import json
from pathlib import Path
import numpy as np
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction
p=Benchmark(); m=Assembly(p,Reconstruction.bulk_interface_hypothesis(p)); eps=1e-6
rows=[]
for axis in range(3):
    vector=m.b.reference.copy(); vector[:,axis]*=1+eps
    stretch=np.linalg.norm(vector,axis=1)/m.b.length-1
    pair_energy=.5*m.coefficient*m.b.length*stretch**2
    density=np.zeros(len(m.x))
    for side in (0,1): np.add.at(density,m.b.ends[:,side],.5*pair_energy/p.spacing**3)
    continuum=.6*p.young*eps**2
    ratio=density/continuum
    rows.append({'axis':'xyz'[axis],'global_energy_ratio':float(pair_energy.sum()/(continuum*len(m.x)*p.spacing**3)),
                 'node_ratio_min':float(ratio.min()),'node_ratio_max':float(ratio.max()),
                 'node_ratio_mean':float(ratio.mean())})
report={'fields':rows,'note':'Ellipsoidal averaging is an approximate surface correction. Calibration equality by construction does not imply exact assembled energy at every particle.'}
Path('outputs/surface_patch.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
