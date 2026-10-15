"""Exercise the installed wheel without importing either original project."""
from pathlib import Path
import json
import sys
import numpy as np

root=Path('/workspace/work/python-build')
sys.path.insert(0,str(root/'wheel-site'))
import pdbeam
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly, Reconstruction
from pdbeam.solver import Simulation

assert Path(pdbeam.__file__).is_relative_to(root/'wheel-site')
p=Benchmark()
model=Assembly(p,Reconstruction.bulk_interface_hypothesis(p))
sim=Simulation(model)
row=sim.advance()
checkpoint=root/'wheel-smoke-checkpoint.npz'
sim.save(checkpoint)
restored=Simulation(model)
restored.restore(checkpoint)
assert sim.step==restored.step==1
for name in ('u','v','su','sv','alive','interface_alive','plastic'):
    np.testing.assert_array_equal(getattr(sim,name),getattr(restored,name))
    assert np.isfinite(getattr(sim,name)).all()
assert row['momentum_residual_kg_m_s']<1e-9
result=dict(imported_package=str(pdbeam.__file__),particles=len(model.x),
            steel_nodes=len(model.sx),unique_concrete_bonds=len(model.b.ends),
            steps=sim.step,checkpoint_roundtrip_exact=True,
            momentum_residual_kg_m_s=row['momentum_residual_kg_m_s'],
            scientific_benchmark_validated=False)
(root/'results'/'wheel-smoke-details.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
