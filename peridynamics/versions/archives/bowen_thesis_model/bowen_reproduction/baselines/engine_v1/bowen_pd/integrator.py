"""Restartable baseline Verlet and equation-based projection alternatives."""
import csv,json,time,platform,sys
from dataclasses import asdict
from pathlib import Path
import numpy as np
from .model import create,evaluate
from .probes import probes
from .impact_contact import project_contact
from .io import save,load,source_hash

class Simulation:
    def __init__(self,mode='legacy',dt=3e-7):
        if mode not in ('legacy','projection'):raise ValueError(mode)
        if dt<=0:raise ValueError('dt must be positive')
        self.mode=mode;self.dt=dt;self.m=create(corrected=mode=='projection');m=self.m
        self.u=np.zeros_like(m.xyz);self.v=self.u.copy();self.us=np.zeros_like(m.steel_xyz);self.vs=self.us.copy()
        self.zh=m.hammer_z_initial;self.vh=m.HAMMER_V0;self.step=0;self.first=None
        self.impulse=np.zeros(3);self.guide_impulse=np.zeros(3)
        self.constraint_ke=0.;self.contact_ke=0.;self.work=0.
        self.f=evaluate(m,self.u,self.us,self.zh,dt,penalty=mode=='legacy')
        self.initial_ke=.5*m.HAMMER_MASS*self.vh**2
    def constrain(self):
        m=self.m;before=self.v.copy();u_before=self.u.copy()
        self.u[m.support_mask,2]=0;self.v[m.support_mask,2]=0
        self.u[m.anchor_x,0]=0;self.v[m.anchor_x,0]=0
        self.u[m.anchor_y,1]=0;self.v[m.anchor_y,1]=0
        self.impulse+=m.CONCRETE_NODE_MASS*np.sum(self.v-before,axis=0)
        self.constraint_ke+=.5*m.CONCRETE_NODE_MASS*np.sum(self.v**2-before**2)
        return float(np.max(abs(self.u-u_before)))
    def advance(self):
        m=self.m;dt=self.dt;f=self.f
        old_u=self.u.copy();old_us=self.us.copy();old_zh=self.zh
        if self.mode=='legacy':
            self.v+=.5*dt*f['fc']/m.CONCRETE_NODE_MASS
            self.vs+=.5*dt*f['fs']/m.STEEL_MASS[:,None]
            self.vh+=.5*dt*f['hf'][2]/m.HAMMER_MASS
        else:
            self.v+=dt*f['fc']/m.CONCRETE_NODE_MASS
            self.vs+=dt*f['fs']/m.STEEL_MASS[:,None]
        self.u+=dt*self.v;self.us+=dt*self.vs;self.zh+=dt*self.vh
        self.constrain()
        if self.mode=='projection':
            self.vh,hf,nc,pen,cke=project_contact(m.xyz,self.u,old_u,self.v,self.zh,self.vh,dt,m.CONCRETE_NODE_MASS,m.HAMMER_MASS,m.hammer_xy,m.HAMMER_RADIUS)
            self.contact_ke+=cke;self.guide_impulse[:2]+=dt*hf[:2]*(-1)
            self.constrain()
        new=evaluate(m,self.u,self.us,self.zh,dt,commit=True,penalty=self.mode=='legacy')
        if self.mode=='legacy':
            self.v+=.5*dt*new['fc']/m.CONCRETE_NODE_MASS
            self.vs+=.5*dt*new['fs']/m.STEEL_MASS[:,None]
            self.vh+=.5*dt*new['hf'][2]/m.HAMMER_MASS
            self.guide_impulse[:2]-=.5*dt*(f['hf'][:2]+new['hf'][:2])
            self.constrain()
        # Actual internal-force work; this is NOT labelled fracture dissipation.
        self.work+=.5*np.sum((f['internal_c']+new['internal_c'])*(self.u-old_u))+.5*np.sum((f['internal_s']+new['internal_s'])*(self.us-old_us))
        self.f=new;self.step+=1
        if self.mode=='projection':self.contact=(hf,nc,pen)
        else:self.contact=(new['hf'],new['nc'],new['pen'])
        if self.contact[1] and self.first is None:self.first=self.step
        return self.record()
    def record(self):
        m=self.m;f=self.f
        hf,nc,pen=getattr(self,'contact',(f['hf'],f['nc'],f['pen']))
        kc=.5*m.CONCRETE_NODE_MASS*np.sum(self.v**2);ks=.5*np.sum(m.STEEL_MASS[:,None]*self.vs**2);kh=.5*m.HAMMER_MASS*self.vh**2
        mom=m.CONCRETE_NODE_MASS*self.v.sum(0)+np.sum(m.STEEL_MASS[:,None]*self.vs,axis=0)+np.array([0,0,m.HAMMER_MASS*self.vh])
        residual=mom-np.array([0,0,m.HAMMER_MASS*m.HAMMER_V0])-self.impulse-self.guide_impulse
        # Baseline uses continuous penalty; projection exchanges impulse explicitly.
        budget=kc+ks+kh+f['ep']-self.initial_ke-self.work-self.constraint_ke-self.contact_ke
        r=dict(step=self.step,time_ms=self.step*self.dt*1000,impact_kN=float(hf[2]/1000),
            max_concrete_mm=float(np.max(np.linalg.norm(self.u,axis=1))*1000),max_steel_mm=float(np.max(np.linalg.norm(self.us,axis=1))*1000),
            n_contact=nc,penetration_mm=pen*1000,broken_concrete=int(m.broken.sum()),broken_interface=int(m.interface_broken.sum()),
            hammer_velocity_mps=self.vh,hammer_ke_J=kh,concrete_ke_J=kc,steel_ke_J=ks,contact_potential_J=f['ep'],
            cumulative_internal_work_J=self.work,constraint_ke_change_J=self.constraint_ke,projection_ke_change_J=self.contact_ke,
            numerical_work_budget_residual_J=float(budget),momentum_residual_kg_mps=float(np.linalg.norm(residual)),
            interface_physical_force_residual_N=f['interface_residual'])
        r.update(probes(m,self.u));return r
    def checkpoint(self,path):
        arrays={k:getattr(self,k) for k in ('u','v','us','vs','impulse','guide_impulse')}
        arrays.update(broken=self.m.broken,prev=self.m.prev,interface_broken=self.m.interface_broken,interface_prev=self.m.interface_prev)
        for k,v in self.f.items():arrays['force_'+k]=np.asarray(v)
        meta=dict(mode=self.mode,dt=self.dt,source_hash=source_hash(),step=self.step,first=self.first,zh=self.zh,vh=self.vh,work=self.work,constraint_ke=self.constraint_ke,contact_ke=self.contact_ke,initial_ke=self.initial_ke)
        save(path,arrays,meta)
    def restore(self,path):
        arrays,meta=load(path,dict(mode=self.mode,dt=self.dt,source_hash=source_hash()))
        for k in ('u','v','us','vs','impulse','guide_impulse'):setattr(self,k,arrays[k])
        for k in ('broken','prev','interface_broken','interface_prev'):setattr(self.m,k,arrays[k])
        for k in ('step','first','zh','vh','work','constraint_ke','contact_ke','initial_ke'):setattr(self,k,meta[k])
        # Do NOT re-evaluate a rate-sensitive force after loading committed stretch:
        # it would use zero rate instead of the previous step's rate.
        self.f={k[6:]:v for k,v in arrays.items() if k.startswith('force_')}

def run(out,steps=550,mode='legacy',dt=3e-7,resume=None,checkpoint_every=250):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);s=Simulation(mode,dt)
    if resume:s.restore(resume)
    manifest=dict(mode=mode,dt=dt,requested_steps=steps,source_hash=source_hash(),parameters=asdict(s.m.p),python=sys.version,platform=platform.platform(),provenance='provisional coupling quadrature and inputs; no Bowen reproduction claim',resume=str(resume) if resume else None)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    start=time.perf_counter();stop=None;peak=-np.inf;deflection=0.;last=None
    with (out/('history_resume.csv' if resume else 'history.csv')).open('w',newline='') as file:
        writer=None
        while s.step<steps:
            try:r=s.advance()
            except (ValueError,FloatingPointError) as e:stop=f'{type(e).__name__}: {e}';break
            if writer is None:writer=csv.DictWriter(file,fieldnames=list(r));writer.writeheader()
            writer.writerow(r);last=r;peak=max(peak,r['impact_kN']);deflection=max(deflection,r['support_relative_bottom_mm'])
            if not np.all(np.isfinite(list(r.values()))):stop='nonfinite state'
            elif r['max_concrete_mm']>50 or r['max_steel_mm']>50:stop='displacement safety limit 50 mm (diagnosis trigger, not physical failure criterion)'
            elif mode=='legacy' and r['penetration_mm']>s.m.dx*500:stop='penalty penetration > dx/2'
            if s.step%checkpoint_every==0 or stop:
                s.checkpoint(out/'checkpoint.npz');file.flush()
            if s.step%500==0 or stop:print(mode,s.step,'force kN',round(r['impact_kN'],3),'bottom mm',round(r['support_relative_bottom_mm'],5),'budget J',round(r['numerical_work_budget_residual_J'],4),flush=True)
            if stop:break
    s.checkpoint(out/'final.npz')
    damage=(np.bincount(s.m.I,weights=s.m.broken,minlength=len(s.m.xyz))+np.bincount(s.m.J,weights=s.m.broken,minlength=len(s.m.xyz)))/s.m.counts
    np.savez(out/'damage.npz',xyz=s.m.xyz,u=s.u,damage=damage,steel_xyz=s.m.steel_xyz,us=s.us)
    summ=dict(steps_completed=s.step,first_contact=s.first,stop=stop,elapsed_s=time.perf_counter()-start,peak_force_kN=float(peak),max_support_relative_bottom_mm=deflection,last=last,complete_requested_interval=s.step==steps)
    (out/'summary.json').write_text(json.dumps(summ,indent=2));return s,summ
