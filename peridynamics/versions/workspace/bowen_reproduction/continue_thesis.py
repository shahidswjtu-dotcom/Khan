"""Continue a diagnosed detached-fragment trigger without changing mechanics.

The 50mm limit remains for material with surviving nonzero-stiffness bonds and
steel. Unbound fragment flight is logged separately. Core source hash stays
identical; external driver hash and changed diagnostic policy are recorded.
"""
import argparse,csv,json,time,hashlib
from pathlib import Path
import numpy as np
from thesis_pd.integrator import Simulation
from thesis_pd.io import source_hash

p=argparse.ArgumentParser();p.add_argument('--resume',required=True);p.add_argument('--out',required=True);p.add_argument('--steps',type=int,default=5000);p.add_argument('--dt',type=float,default=3e-7);a=p.parse_args()
out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
s=Simulation('family_projection',a.dt);s.restore(a.resume);start_step=s.step;m=s.m
(out/'manifest.json').write_text(json.dumps(dict(source_hash=source_hash(),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),mode=s.mode,dt=s.dt,resume=a.resume,start_step=start_step,requested_step=a.steps,stop_policy='50mm attached concrete or steel; detached flight tracked separately; stop at nonfinite diagnostics',justification='Prior50mm trigger had0 surviving nonzero-stiffness bonds,0.94J particle kinetic energy; attached maximum28.56mm. No mechanics, material or integration changes.',model_hash=s.model_hash()),indent=2))
start=time.perf_counter();stop=None;last=None
with (out/'history.csv').open('w',newline='') as f:
 writer=None
 while s.step<a.steps:
  r=s.advance()
  alive=(~m.broken)&(m.LEN<m.delta*(1-1e-12))
  degree=np.bincount(m.I[alive],minlength=len(m.xyz))+np.bincount(m.J[alive],minlength=len(m.xyz))
  displacement=np.linalg.norm(s.u,axis=1)*1000
  r['max_attached_concrete_mm']=float(displacement[degree>0].max(initial=0))
  r['detached_particle_count']=int(np.sum(degree==0))
  if writer is None:writer=csv.DictWriter(f,fieldnames=list(r));writer.writeheader()
  writer.writerow(r);last=r
  if not np.all(np.isfinite(list(r.values()))):stop='Nonfinite diagnostics'
  elif r['max_attached_concrete_mm']>50 or r['max_steel_mm']>50:stop='50mm attached-material/steel diagnostic stop'
  if s.step%250==0 or stop:s.checkpoint(out/'checkpoint.npz');f.flush()
  if s.step%500==0 or stop:print('continued',s.step,'bottom mm',r['support_relative_bottom_mm'],'attached max mm',r['max_attached_concrete_mm'],flush=True)
  if stop:break
s.checkpoint(out/'final.npz')
damage=(np.bincount(m.I,weights=m.broken,minlength=len(m.xyz))+np.bincount(m.J,weights=m.broken,minlength=len(m.xyz)))/m.counts
np.savez(out/'damage.npz',xyz=m.xyz,u=s.u,damage=damage,steel_xyz=m.steel_xyz,us=s.us)
summary=dict(start_step=start_step,steps_completed=s.step,complete_requested_interval=s.step==a.steps,stop=stop,last=last,elapsed_s=time.perf_counter()-start)
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(summary)
