import argparse,csv,json,time
from pathlib import Path
import numpy as np
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction
from pdbeam.solver import Simulation

parser=argparse.ArgumentParser()
parser.add_argument('--steps',type=int,default=5000)
parser.add_argument('--dt',type=float,default=3e-7)
parser.add_argument('--out',default='outputs/fresh_bulk_hypothesis')
parser.add_argument('--restart')
parser.add_argument('--accept-bulk-interface-hypothesis',action='store_true',required=True,
                    help='Explicitly selects unverified interface reconstruction, not author inputs')
a=parser.parse_args(); p=Benchmark(dt=a.dt)
out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
model=Assembly(p,Reconstruction.bulk_interface_hypothesis(p)); sim=Simulation(model)
if a.restart: sim.restore(a.restart)
(out/'manifest.json').write_text(model.manifest+'\n')
start=time.monotonic(); peak=-float('inf'); displacement=-float('inf'); last=None
with (out/'history.csv').open('w',newline='') as stream:
    writer=None
    while sim.step<a.steps:
        last=sim.advance()
        if writer is None:
            writer=csv.DictWriter(stream,fieldnames=list(last)); writer.writeheader()
        writer.writerow(last)
        peak=max(peak,last['force_N']); displacement=max(displacement,last['bottom_midspan_down_m'])
        if sim.step%100==0:
            stream.flush()
            print(json.dumps({'step':sim.step,'peak_kN':peak/1000,'bottom_mm':displacement*1000,'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
        if sim.step%1000==0: sim.save(out/'checkpoint.npz')
sim.save(out/'checkpoint.npz')
np.savez_compressed(out/'final_state.npz',positions=model.x,displacement=sim.u,steel_positions=model.sx,steel_displacement=sim.su,ends=model.b.ends,alive=sim.alive)
summary=dict(final=last,segment_peak_force_kN=peak/1000,segment_max_bottom_midspan_mm=displacement*1000,
             elapsed_s=time.monotonic()-start,original_author_model_reproduced=False,
             interface_definition=model.choices.interface_definition)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
