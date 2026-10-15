"""Build documented concrete discretization; does not claim a full RC solver."""
import json
from pathlib import Path
import numpy as np
from pdbeam.core import Benchmark, particles, family, concrete_force

p=Benchmark()
x=particles(p)
b=family(x,p)
f,flags=concrete_force(x,np.zeros_like(x),b,p,np.ones(len(b.ends),dtype=bool))
out=Path(__file__).parent/'outputs'; out.mkdir(exist_ok=True)
np.savez_compressed(out/'reference_mesh.npz',positions=x,ends=b.ends,
                    reference_length=b.length,stiffness=b.stiffness,
                    volume_fraction=b.volume_fraction)
report=dict(particles=len(x),unique_bonds=len(b.ends),
            concrete_mass_kg=len(x)*p.spacing**3*p.density,
            max_initial_force_N=float(np.max(abs(f))),
            bar_height_from_bottom_m=p.bar_height,
            timestep_s=p.dt,stage19_dependencies=[],
            full_impact_ready=False)
(out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
