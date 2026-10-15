"""Restartable explicit solver for a declared Chapter 3 reconstruction."""
import numpy as np
from scipy.spatial import cKDTree
from .core import project_contact
from .operators import concrete, coupling, rods, repulsion
from .energy import threshold_ledger, concrete_potential, interface_rupture, rod_potential, interface_potential

class Simulation:
    def __init__(self,model):
        self.m=model; self.step=0
        self.u=np.zeros_like(model.x); self.previous=self.u.copy(); self.v=self.u.copy()
        self.su=np.zeros_like(model.sx); self.sv=self.su.copy()
        self.alive=np.ones(len(model.b.ends),dtype=bool)
        self.interface_alive=np.ones(len(model.owners),dtype=bool)
        self.plastic=np.zeros(len(model.rod_ends))
        p=model.p
        self.hz=p.height+p.hammer_radius+model.choices.initial_gap
        self.hv=-p.hammer_speed
        self.external_impulse=np.zeros(3)
        self.contact_energy=0.; self.support_energy=0.; self.internal_work=0.; self.plastic_dissipation=0.
        self.old_force=np.zeros_like(self.u); self.old_steel_force=np.zeros_like(self.su)
        self.rate_exchange=0.; self.rupture_loss=0.
        self.threshold_a=np.full(len(model.b.ends),p.onset)
        self.threshold_b=np.full(len(model.b.ends),p.terminal)
        self.initial_energy=.5*p.hammer_mass*p.hammer_speed**2
        self.initial_momentum=np.array([0.,0.,-p.hammer_mass*p.hammer_speed])

    def kinetic(self):
        m=self.m
        return .5*np.sum(m.mass[:,None]*self.v**2)+.5*np.sum(m.steel_mass[:,None]*self.sv**2)+.5*m.p.hammer_mass*self.hv**2

    def advance(self):
        m=self.m; p=m.p; c=m.choices
        exchange,loss=threshold_ledger(m.x,self.u,self.previous,m.b.ends,m.b.length,m.coefficient,self.alive,self.threshold_a,self.threshold_b,p.onset,p.terminal,p.dt,c.rate_enabled,c.irreversible)
        self.rate_exchange+=exchange; self.rupture_loss+=loss
        self.rupture_loss+=interface_rupture(m.x,self.u,m.sx,self.su,m.owners,m.stencils,m.weights,m.interface_length,m.interface_coefficient,self.interface_alive,c.interface_onset,c.interface_terminal,c.irreversible)
        fc,ec=concrete(m.x,self.u,self.previous,m.b.ends,m.b.length,m.coefficient,self.alive,p.onset,p.terminal,p.dt,c.irreversible,c.rate_enabled)
        fic,fis,ei=coupling(m.x,self.u,m.sx,self.su,m.owners,m.stencils,m.weights,m.interface_length,m.interface_coefficient,self.interface_alive,c.interface_onset,c.interface_terminal,c.irreversible)
        fr,er,diss=rods(m.sx,self.su,m.rod_ends,m.rod_length,m.area,p.steel_young,p.steel_yield,self.plastic)
        self.plastic_dissipation+=diss
        current=m.x+self.u
        close=cKDTree(current).query_pairs(1.35*p.spacing,output_type='ndarray')
        fsr,esr=repulsion(m.x,current,close,p.spacing,p.horizon,p.young)
        force=fc+fic+fsr; steel_force=fis+fr
        old=self.u.copy(); olds=self.su.copy()
        self.v+=p.dt*force/m.mass[:,None]; self.sv+=p.dt*steel_force/m.steel_mass[:,None]
        predicted=old+p.dt*self.v
        self.su+=p.dt*self.sv
        self.hz+=p.dt*self.hv
        before_contact=self.kinetic()
        corrected,velocity,reaction=project_contact(m.x+old,m.x+predicted,np.array([p.length/2,p.width/2,self.hz]),p.hammer_radius,m.mass,p.dt)
        self.u=corrected-m.x; self.v=velocity
        self.hv+=reaction[2]/p.hammer_mass
        self.external_impulse[:2]-=reaction[:2] # axial hammer guide
        self.contact_energy+=self.kinetic()-before_contact
        before_support=self.kinetic()
        support_impulse=np.where(m.fixed,-m.mass[:,None]*self.v,0)
        self.external_impulse+=support_impulse.sum(axis=0)
        self.u[m.fixed]=0.; self.v[m.fixed]=0.
        self.support_energy+=self.kinetic()-before_support
        # Path work diagnostic. Not a substitute for fracture/rate energy accounting.
        self.internal_work+=np.sum(.5*(self.old_force+force)*(self.u-old))+np.sum(.5*(self.old_steel_force+steel_force)*(self.su-olds))
        self.old_force=force; self.old_steel_force=steel_force
        self.previous=old; self.step+=1
        momentum=np.sum(m.mass[:,None]*self.v,axis=0)+np.sum(m.steel_mass[:,None]*self.sv,axis=0)+np.array([0,0,p.hammer_mass*self.hv])
        kinetic=self.kinetic()
        # End-state stored energy uses the force step's frozen thresholds.
        stored=float('nan'); residual=float('nan')
        if self.step%100==0:
            stored=concrete_potential(m.x,self.u,m.b.ends,m.b.length,m.coefficient,self.alive,self.threshold_a,self.threshold_b,c.irreversible)
            ip=interface_potential(m.x,self.u,m.sx,self.su,m.owners,m.stencils,m.weights,m.interface_length,m.interface_coefficient,self.interface_alive,c.interface_onset,c.interface_terminal,c.irreversible)
            rp=rod_potential(m.sx,self.su,m.rod_ends,m.rod_length,m.area,p.steel_young,self.plastic)
            updated=m.x+self.u
            near=cKDTree(updated).query_pairs(1.35*p.spacing,output_type='ndarray')
            _,sp=repulsion(m.x,updated,near,p.spacing,p.horizon,p.young)
            stored+=ip+rp+sp
            residual=kinetic+stored+self.rupture_loss+self.plastic_dissipation-self.initial_energy-self.rate_exchange-self.contact_energy-self.support_energy
        row=dict(step=self.step,time_s=self.step*p.dt,force_N=reaction[2]/p.dt,
                 bottom_midspan_down_m=-float(self.u[m.probe,2].mean()),
                 max_particle_displacement_m=float(np.linalg.norm(self.u,axis=1).max()),
                 kinetic_J=float(kinetic),stored_energy_J=float(stored),energy_balance_residual_J=float(residual),rate_potential_exchange_J=float(self.rate_exchange),rupture_deleted_energy_J=float(self.rupture_loss),potential_pre_step_J=float(ec+ei+er+esr),
                 plastic_dissipation_J=float(self.plastic_dissipation),
                 contact_energy_change_J=float(self.contact_energy),support_energy_change_J=float(self.support_energy),
                 internal_path_work_J=float(self.internal_work),
                 work_balance_residual_J=float(kinetic-self.initial_energy-self.internal_work-self.contact_energy-self.support_energy),
                 momentum_residual_kg_m_s=float(np.linalg.norm(momentum-self.initial_momentum-self.external_impulse)),
                 broken_bonds=int((~self.alive).sum()),broken_interface=int((~self.interface_alive).sum()),hammer_velocity_m_s=float(self.hv))
        if not np.isfinite(kinetic) or np.max(np.linalg.norm(self.v,axis=1))>1000:
            raise FloatingPointError('Nonfinite state or particle speed above diagnostic 1000 m/s limit')
        return row

    def save(self,path):
        arrays=['u','previous','v','su','sv','alive','interface_alive','plastic','external_impulse','old_force','old_steel_force','threshold_a','threshold_b']
        scalars=['step','hz','hv','contact_energy','support_energy','internal_work','plastic_dissipation','rate_exchange','rupture_loss']
        np.savez_compressed(path,manifest=self.m.manifest,**{k:getattr(self,k) for k in arrays+scalars})

    def restore(self,path):
        with np.load(path,allow_pickle=False) as data:
            if str(data['manifest'])!=self.m.manifest:
                raise ValueError('Solver/reconstruction manifest differs')
            for k in data.files:
                if k=='manifest': continue
                value=data[k]
                setattr(self,k,value.copy() if value.ndim else value.item())
