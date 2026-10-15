from time import perf_counter
from pathlib import Path
import csv, json

# Benchmark force evaluations independent of time integration.
# Note: repeated force evaluations at the undeformed state are not a stability test.
reset_stage8()
u0=np.zeros_like(xyz); s0=np.zeros_like(steel_xyz)
t0=perf_counter()
for _ in range(2):
    _=source_shaped_forces(u0,s0,float(hammer_z_initial),p.dt,commit=False)
force_eval_seconds=(perf_counter()-t0)/2
print(f'Full-mesh force evaluation: {force_eval_seconds:.3f} s')
print(f'Estimated 350-step minimum at two evaluations/step: {700*force_eval_seconds/60:.1f} min (ignores overhead)')
print('This timing depends on your CPU and is only a rough lower-bound estimate.')