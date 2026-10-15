def bond_scalar_force(stretch, stretch_rate, distance, broken=False, improved=True):
    """Return scalar bond force density and irreversible failure flag.

    The demonstration assumes instantaneous elastic/softening envelope forces.
    It does not implement an unstated unloading path for partially softened bonds.
    """
    s,sd,r,d=np.broadcast_arrays(np.asarray(stretch,float),np.asarray(stretch_rate,float),
                                np.asarray(distance,float),np.asarray(broken,bool))
    s0,sc=thresholds_at_rate(sd)
    c=micromodulus(r,improved=improved)
    fail=d | (s>=sc)
    elastic=c*s
    # For positive stretches above s0: linear decrease to zero at sc.
    soft=c*s0*np.clip((sc-s)/np.maximum(sc-s0,1.e-30),0.0,1.0)
    f=np.where(s<=s0,elastic,soft)
    return np.where(fail,0.0,f),fail

# A lightweight state holder so persistent failure cannot heal after unloading.
class BondHistory:
    def __init__(self,n):
        self.broken=np.zeros(n,dtype=bool)
        self.max_stretch=np.full(n,-np.inf)
    def update(self,stretch,rate,length):
        stretch=np.asarray(stretch,dtype=float)
        self.max_stretch=np.maximum(self.max_stretch,stretch)
        f,new_broken=bond_scalar_force(stretch,rate,length,self.broken)
        self.broken |= new_broken
        return f
