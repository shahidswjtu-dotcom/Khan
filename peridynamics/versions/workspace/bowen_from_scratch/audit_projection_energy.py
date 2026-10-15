"""Audit a checkpoint continuation without altering the mechanical operator.

Projection changes configuration and hence stored energy. Report its change
separately from kinetic contact/constraint changes recorded by the base solver.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from pdbeam.core import Benchmark,project_contact
from pdbeam.assembly import Assembly,Reconstruction
from pdbeam.solver import Simulation
from pdbeam.energy import concrete_potential,interface_potential
from pdbeam.operators import repulsion


def potential_change(m,sim,before,after):
    changed=np.any(before!=after,axis=1)
    ids=np.flatnonzero(changed)
    if not len(ids): return 0.
    selected=changed[m.b.ends[:,0]]|changed[m.b.ends[:,1]]
    pairs=m.b.ends[selected]; r0=m.b.length[selected]; coeff=m.coefficient[selected]
    alive=sim.alive[selected]; a=sim.threshold_a[selected]; b=sim.threshold_b[selected]
    before_u=before-m.x; after_u=after-m.x
    energy=concrete_potential(m.x,after_u,pairs,r0,coeff,alive,a,b,m.choices.irreversible)-concrete_potential(m.x,before_u,pairs,r0,coeff,alive,a,b,m.choices.irreversible)
    selected=np.any(changed[m.stencils],axis=1)
    params=(m.sx,sim.su,m.owners[selected],m.stencils[selected],m.weights[selected],m.interface_length[selected],m.interface_coefficient[selected],sim.interface_alive[selected],m.choices.interface_onset,m.choices.interface_terminal,m.choices.irreversible)
    energy+=interface_potential(m.x,after_u,*params)-interface_potential(m.x,before_u,*params)
    # Union of neighbours near moved nodes before/after; unchanged pairs cancel.
    keys=[]; n=len(m.x)
    for coordinates in (before,after):
        neighbours=cKDTree(coordinates).query_ball_point(coordinates[ids],1.35*m.p.spacing)
        for i,adjacent in zip(ids,neighbours):
            for j in adjacent:
                if i!=j: keys.append(min(i,j)*n+max(i,j))
    keys=np.unique(np.asarray(keys,dtype=np.int64))
    near=np.column_stack((keys//n,keys%n))
    args=(near,m.p.spacing,m.p.horizon,m.p.young)
    energy+=repulsion(m.x,after,*args)[1]-repulsion(m.x,before,*args)[1]
    return float(energy)

class AuditedSimulation(Simulation):
    def __init__(self,m):
        super().__init__(m)
        self.projection_potential=0.; self.constraint_potential=0.

    def advance(self):
        m=self.m; p=m.p; old=self.u.copy(); velocity=self.v.copy()
        row=super().advance()
        predicted=m.x+old+p.dt*(velocity+p.dt*self.old_force/m.mass[:,None])
        corrected,_,_=project_contact(m.x+old,predicted,np.array([p.length/2,p.width/2,self.hz]),p.hammer_radius,m.mass,p.dt)
        self.projection_potential+=potential_change(m,self,predicted,corrected)
        self.constraint_potential+=potential_change(m,self,corrected,m.x+self.u)
        row['contact_projection_potential_change_J']=self.projection_potential
        row['constraint_projection_potential_change_J']=self.constraint_potential
        row['augmented_energy_residual_J']=row['energy_balance_residual_J']-self.projection_potential-self.constraint_potential
        return row

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--checkpoint',required=True)
    parser.add_argument('--steps',type=int,default=100)
    parser.add_argument('--dt',type=float,default=3e-7)
    parser.add_argument('--out',default='outputs/projection_energy_audit')
    a=parser.parse_args(); p=Benchmark(dt=a.dt)
    m=Assembly(p,Reconstruction.bulk_interface_hypothesis(p))
    sim=AuditedSimulation(m); sim.restore(a.checkpoint)
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    start=sim.step
    with (out/'history.csv').open('w',newline='') as f:
        writer=None
        for _ in range(a.steps):
            row=sim.advance()
            if writer is None:
                writer=csv.DictWriter(f,fieldnames=row.keys()); writer.writeheader()
            writer.writerow(row)
    result=dict(start_step=start,end_step=sim.step,final=row,
                note='Projection potential sums cover this segment only; augmented residual retains the pre-segment residual offset. Compare changes, not its absolute value.')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
