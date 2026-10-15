from dataclasses import dataclass, asdict
from itertools import product
import json
import hashlib
from pathlib import Path
import numpy as np
from .core import Benchmark, particles, family, interpolate_weights

@dataclass(frozen=True)
class Reconstruction:
    # These are declared hypotheses, not recovered Chapter 3 input values.
    interface_definition: str = 'UNRESOLVED'
    interface_stiffness: float = 0.
    interface_onset: float = 0.
    interface_terminal: float = 0.
    bar_y: tuple = (.020,.050)
    bar_x_start: float = .025
    bar_x_end: float = .575
    hook_rise: float = 0.
    support_width: float = .010
    initial_gap: float = 0.
    irreversible: bool = True
    rate_enabled: bool = True

    @classmethod
    def bulk_interface_hypothesis(cls,p):
        return cls(hook_rise=.060,interface_definition='hypothesis: concrete weighted micromodulus and thresholds; rod lumped control volume',
                   interface_stiffness=72*p.young/(np.pi*p.horizon**4),
                   interface_onset=p.onset,interface_terminal=p.terminal)

class Assembly:
    def __init__(self,p=Benchmark(),choices=Reconstruction()):
        if choices.interface_definition=='UNRESOLVED':
            raise ValueError('No recovered Chapter 3 interface inputs: explicitly select and label a hypothesis')
        if not (choices.interface_stiffness>0 and 0<choices.interface_onset<choices.interface_terminal):
            raise ValueError('Invalid interface parameters')
        self.p=p; self.choices=choices; self.x=particles(p); self.b=family(self.x,p)
        # Source-energy calibration: 0.1% uniaxial strain without lateral strain.
        energy=np.zeros((len(self.x),3)); eps=.001
        for axis in range(3):
            deformed=self.b.reference.copy(); deformed[:,axis]*=1+eps
            s=np.linalg.norm(deformed,axis=1)/self.b.length-1
            e=.25*self.b.stiffness*s*s*self.b.length*p.spacing**3*self.b.volume_fraction
            np.add.at(energy[:,axis],self.b.ends[:,0],e)
            np.add.at(energy[:,axis],self.b.ends[:,1],e)
        self.node_surface=.6*p.young*eps**2/energy
        midpoint=.5*(self.node_surface[self.b.ends[:,0]]+self.node_surface[self.b.ends[:,1]])
        direction=self.b.reference/self.b.length[:,None]
        surface=1/np.sqrt(np.sum((direction/midpoint)**2,axis=1))
        self.coefficient=self.b.stiffness*p.spacing**6*self.b.volume_fraction*surface
        n=round((choices.bar_x_end-choices.bar_x_start)/p.spacing)+1
        xs=np.linspace(choices.bar_x_start,choices.bar_x_end,n)
        if choices.hook_rise<0 or p.bar_height+choices.hook_rise>=p.height:
            raise ValueError('End hooks must remain inside the beam')
        count=round(choices.hook_rise/p.spacing)
        if not np.isclose(count*p.spacing,choices.hook_rise):
            raise ValueError('Hook rise must be an integral number of rod spacings')
        paths=[]
        for y in choices.bar_y:
            path=[(choices.bar_x_start,y,p.bar_height+k*p.spacing) for k in range(count,0,-1)]
            path.extend((x,y,p.bar_height) for x in xs)
            path.extend((choices.bar_x_end,y,p.bar_height+k*p.spacing) for k in range(1,count+1))
            paths.append(path)
        n=len(paths[0])
        self.sx=np.array(paths).reshape(-1,3)
        self.rod_ends=np.array([(bar*n+k,bar*n+k+1) for bar in range(2) for k in range(n-1)])
        self.rod_length=np.linalg.norm(self.sx[self.rod_ends[:,1]]-self.sx[self.rod_ends[:,0]],axis=1)
        self.area=np.pi*p.bar_diameter**2/4
        self.steel_volume=np.zeros(len(self.sx))
        for side in range(2):
            np.add.at(self.steel_volume,self.rod_ends[:,side],.5*self.area*self.rod_length)
        self.steel_mass=p.steel_density*self.steel_volume
        owners=[]; stencils=[]; weights=[]; lengths=[]; coefficients=[]
        shape=np.rint(np.array([p.length,p.width,p.height])/p.spacing).astype(int)
        offsets=np.array(list(product(range(-3,4),repeat=3)))*p.spacing
        offsets=offsets[(np.linalg.norm(offsets,axis=1)>0)&(np.linalg.norm(offsets,axis=1)<=p.horizon*(1+1e-12))]
        bits=np.array(list(product(range(2),repeat=3)))
        for a,position in enumerate(self.sx):
            for delta in offsets:
                point=position+delta
                q=(point-.5*p.spacing)/p.spacing
                low=np.floor(q+1e-10).astype(int)
                if np.any(low<0) or np.any(low>=shape-1):
                    continue
                ijk=low+bits
                ids=np.ravel_multi_index(ijk.T,tuple(shape))
                w=interpolate_weights(point,self.x[ids])
                r=np.linalg.norm(delta)
                vc=np.clip((p.horizon+.5*p.spacing-r)/p.spacing,0,1)
                owners.append(a); stencils.append(ids); weights.append(w); lengths.append(r)
                coefficients.append(choices.interface_stiffness*(1-(r/p.horizon)**2)**2*self.steel_volume[a]*p.spacing**3*vc)
        self.owners=np.array(owners); self.stencils=np.array(stencils); self.weights=np.array(weights)
        self.interface_length=np.array(lengths); self.interface_coefficient=np.array(coefficients)
        self.fixed=np.zeros_like(self.x,dtype=bool)
        bottom=self.x[:,2]<p.spacing
        for support in (.05,.55):
            self.fixed[:,2]|=bottom&(abs(self.x[:,0]-support)<=choices.support_width/2+1e-12)
        # Minimal in-plane rigid-body anchors, explicitly reconstructed.
        left=np.argmin(np.linalg.norm(self.x-np.array([.05,.035,.0025]),axis=1))
        right=np.argmin(np.linalg.norm(self.x-np.array([.55,.035,.0025]),axis=1))
        self.fixed[left,:2]=True; self.fixed[right,1]=True
        self.mass=np.full(len(self.x),p.density*p.spacing**3)
        self.probe=(abs(self.x[:,0]-.3)<=p.spacing/2+1e-12)&bottom
        self.manifest=json.dumps({'benchmark':asdict(p),'reconstruction':asdict(choices),
                                  'solver_version':'fresh-1','source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(Path(__file__).parent.glob('*.py'))}},sort_keys=True)
