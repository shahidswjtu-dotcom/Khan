@njit
def accumulate_unique_pairs(iidx,jidx, vec, amplitudes, n_particles):
    f=np.zeros((n_particles,3),dtype=np.float64)
    for k in range(len(iidx)):
        i=iidx[k];j=jidx[k]
        fx=amplitudes[k]*vec[k,0];fy=amplitudes[k]*vec[k,1];fz=amplitudes[k]*vec[k,2]
        f[i,0]+=fx;f[i,1]+=fy;f[i,2]+=fz
        f[j,0]-=fx;f[j,1]-=fy;f[j,2]-=fz
    return f

class DiscretePDConcrete:
    def __init__(self, xyz, I, J, lengths, weights):
        self.xyz=xyz;self.I=I;self.J=J;self.lengths=lengths
        self.weights=weights;self.broken=np.zeros(len(I),dtype=np.bool_)
        self.prev_stretch=np.zeros(len(I),dtype=np.float64)
    def evaluate(self, displacement, dt, commit=False, override_rate=None):
        deformed=self.xyz+np.asarray(displacement,float)
        d=deformed[self.J]-deformed[self.I]
        cur=np.linalg.norm(d,axis=1)
        stretch=(cur-self.lengths)/self.lengths
        rate=(stretch-self.prev_stretch)/dt if override_rate is None else np.broadcast_to(override_rate,stretch.shape)
        scalar,failed=bond_scalar_force(stretch,rate,self.lengths,self.broken)
        amp=scalar*(V**2)*self.weights
        direction=d/cur[:,None]
        F=accumulate_unique_pairs(self.I,self.J,direction,amp,len(self.xyz))
        if commit:
            self.broken |= failed
            self.prev_stretch[:]=stretch
        return F,stretch,failed
    def damage(self):
        # Per-node unweighted broken-bond fraction; different from weighted material-point damage in thesis.
        return np.bincount(np.r_[self.I,self.J],weights=np.r_[self.broken,self.broken],minlength=len(self.xyz))/counts

model=DiscretePDConcrete(xyz,I,J,LEN,PAIR_WEIGHT)
print('Ready:',len(model.I),'bonds; Numba:',HAVE_NUMBA)
