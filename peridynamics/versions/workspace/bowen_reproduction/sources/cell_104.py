import numpy as np
from dataclasses import dataclass

@dataclass(frozen=True)
class AppendixCParameters:
    bcc: float
    s01: float
    sc1: float
    delta: float
    vol: float

class AppendixCPublishedKernel:
    """Vectorized form of Appendix C using explicit 0-based steel/concrete indexing."""
    def __init__(self, steel_xyz, concrete_xyz, steel_id, concrete8_id,
                 weights8, virtual_xyz, params: AppendixCParameters):
        self.xs = np.asarray(steel_xyz, float)
        self.xc = np.asarray(concrete_xyz, float)
        self.st = np.asarray(steel_id, np.int64)
        self.idx = np.asarray(concrete8_id, np.int64)
        self.w = np.asarray(weights8, float)
        self.xv = np.asarray(virtual_xyz, float)
        self.p = params
        if self.idx.shape != (len(self.st), 8) or self.w.shape != self.idx.shape:
            raise ValueError('Each coupling entry must have exactly 8 concrete indices and weights')
        if self.xv.shape != (len(self.st), 3): raise ValueError('virtual_xyz must have shape (K,3)')
        if np.any(self.st<0) or np.any(self.st>=len(self.xs)) or np.any(self.idx<0) or np.any(self.idx>=len(self.xc)):
            raise IndexError('Input connectivity is outside the respective steel or concrete array')
        if np.max(np.abs(self.w.sum(axis=1)-1))>1e-10: raise ValueError('Shape weights must sum to one')
        self.reference = self.xv-self.xs[self.st]
        self.idist = np.linalg.norm(self.reference, axis=1)
        if np.any(self.idist<=0): raise ValueError('PDROD bonds cannot have zero initial length')
        if not (params.bcc>0 and params.delta>0 and params.vol>0 and 0<params.s01<params.sc1):
            raise ValueError('Invalid PDROD constants')
        if np.any(self.idist>params.delta): raise ValueError('Neighbors must be within delta')
        self.failfail = np.ones(len(self.st), dtype=np.int8)

    def evaluate(self, uc, us, old_uc=None, old_us=None, dt=None,
                 mode='literal', commit=False):
        if mode not in ('literal','irreversible_corrected'):
            raise ValueError('mode must be literal or irreversible_corrected')
        uc=np.asarray(uc,float); us=np.asarray(us,float)
        if uc.shape != self.xc.shape or us.shape != self.xs.shape:
            raise ValueError('Displacements must match geometry shape')
        dispdisp=np.einsum('ka,kad->kd',self.w, uc[self.idx])
        direction=self.reference+dispdisp-us[self.st]
        nlength=np.linalg.norm(direction,axis=1)
        if np.any(nlength<1e-15): raise FloatingPointError('Current interface length is zero')
        sssss=(nlength-self.idist)/self.idist
        if old_uc is not None and old_us is not None and dt is not None:
            if dt<=0: raise ValueError('dt must be positive')
            old_interp=np.einsum('ka,kad->kd',self.w,np.asarray(old_uc)[self.idx])
            oldnlength=np.linalg.norm(self.reference+old_interp-np.asarray(old_us)[self.st],axis=1)
            s1=np.abs((nlength-oldnlength)/self.idist/dt)
        else:
            oldnlength=np.full_like(nlength,np.nan); s1=np.full_like(nlength,np.nan)
        bc=self.p.bcc*(1-(self.idist/self.p.delta)**2)**2
        scalar=np.where(sssss<=self.p.s01,bc*sssss,
                    np.where(sssss<=self.p.sc1,
                        bc*self.p.s01*(self.p.sc1-sssss)/(self.p.sc1-self.p.s01),0.0))
        new_failed=sssss>self.p.sc1
        if mode=='irreversible_corrected':
            scalar=np.where((self.failfail==0)|new_failed,0.,scalar)
        # Appendix C: dforce = bc*(stretch law)*vol*(deformed vector)/nlength
        dforce=(scalar*self.p.vol)[:,None]*(direction/nlength[:,None])
        fs=np.zeros_like(us); fc=np.zeros_like(uc)
        np.add.at(fs,self.st,dforce)
        for a in range(8): np.add.at(fc,self.idx[:,a],-self.w[:,a,None]*dforce)
        # Printed code computes failure flag after dforce and does not use it above.
        if commit: self.failfail[new_failed]=0
        dmgpar1=float(np.sum(self.failfail*self.p.vol))
        dmgpar2=float(len(self.failfail)*self.p.vol)
        return dict(concrete_force=fc,steel_force=fs,stretch=sssss,stretch_rate=s1,
                    oldnlength=oldnlength,current_length=nlength,dforce=dforce,
                    new_failed=new_failed,failfail=self.failfail.copy(),
                    sample_damage=(1-dmgpar1/dmgpar2 if dmgpar2 else np.nan))
    def reset(self): self.failfail[:]=1
