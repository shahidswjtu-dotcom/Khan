"""Separate supported-material and detached-particle displacement diagnostics."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction

root=Path(__file__).parent
parser=argparse.ArgumentParser()
parser.add_argument('--checkpoint',default='outputs/hooked_dt3e7/checkpoint.npz')
parser.add_argument('--out',default='outputs/hooked_report/connectivity.json')
a=parser.parse_args()
with np.load(root/a.checkpoint,allow_pickle=False) as d:
    manifest=json.loads(str(d['manifest']))
p=Benchmark(**manifest['benchmark'])
choices=dict(manifest['reconstruction']); choices['bar_y']=tuple(choices['bar_y'])
m=Assembly(p,Reconstruction(**choices))
with np.load(root/a.checkpoint,allow_pickle=False) as d:
    alive=d['alive']; ia=d['interface_alive']; u=d['u']; su=d['su']; step=int(d['step']); terminal=d['threshold_b']
n=len(m.x); ns=len(m.sx)
current=m.x+u
length=np.linalg.norm(current[m.b.ends[:,1]]-current[m.b.ends[:,0]],axis=1)
active_concrete=(length/m.b.length-1<=terminal)&(m.coefficient>0)
virtual=np.sum(m.weights[:,:,None]*current[m.stencils],axis=1)
interface_stretch=np.linalg.norm(virtual-m.sx[m.owners]-su[m.owners],axis=1)/m.interface_length-1
active=(interface_stretch<=m.choices.interface_terminal)&(m.interface_coefficient>0)
if m.choices.irreversible:
    active_concrete&=alive; active&=ia
parts=[m.b.ends[active_concrete],m.rod_ends+n]
for corner in range(8):
    selected=active&(m.weights[:,corner]>1e-12)
    parts.append(np.column_stack((m.stencils[selected,corner],n+m.owners[selected])))
ends=np.concatenate(parts)
graph=coo_matrix((np.ones(len(ends)),(ends[:,0],ends[:,1])),shape=(n+ns,n+ns)).tocsr()
count,labels=connected_components(graph,directed=False)
supported_labels=np.unique(labels[:n][m.fixed[:,2]])
supported=np.isin(labels,supported_labels)
probe_supported=m.probe&supported[:n]
nsingle=np.bincount(labels)[labels]==1
report={'step':step,'components':count,'supported_concrete_nodes':int(supported[:n].sum()),
    'detached_concrete_nodes':int((~supported[:n]).sum()),
    'isolated_concrete_nodes':int(nsingle[:n].sum()),
    'probe_nodes':int(m.probe.sum()),'supported_probe_nodes':int(probe_supported.sum()),
    'mean_bottom_probe_down_mm':-float(u[m.probe,2].mean())*1000,
    'mean_supported_bottom_probe_down_mm':-float(u[probe_supported,2].mean())*1000 if probe_supported.any() else None,
    'maximum_supported_concrete_displacement_mm':float(np.linalg.norm(u[supported[:n]],axis=1).max())*1000,
    'maximum_supported_downward_displacement_mm':float((-u[supported[:n],2]).max())*1000,
    'maximum_all_concrete_displacement_mm':float(np.linalg.norm(u,axis=1).max())*1000,
    'steel_midspan_down_mm':[-float(su[np.argmin(np.linalg.norm(m.sx-np.array([.3,y,.02]),axis=1)),2])*1000 for y in m.choices.bar_y],
    'definition':'Graph includes nonzero-stiffness bonds whose current stretch is below frozen terminal thresholds, positive-weight active interface links, and rods. Prior failure flags gate links only for irreversible history. Supported components intersect vertical support particles.'}
destination=root/a.out
destination.parent.mkdir(parents=True,exist_ok=True)
destination.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
